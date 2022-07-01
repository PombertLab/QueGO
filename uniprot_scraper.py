#!/usr/bin/python

name = "uniprot_scraper.py"
version = "1.0.0"
updated = "2022-07-01"

usage = f"""\n
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	Navigates UniProt and aquires the metadata of desired proteins as well as 
		downloads the corresponding 3D structures and FASTA files

COMMAND		{name} \\
		  -v \\
		  -s \\
		  -ds \\
		  -df \\
		  -m X-ray \\
		  -c '(reviewed:yes)AND(organism:"Drosophila melanogaster (Fruit fly) [7227]")AND(structure_3d:true)'


OPTIONS
-g (--go_keyword)		Gene ontoloy keyword to search
-v (--verified_only)		Get downloads for genes that have been verified by UniProt
-s (--structures)		Find accessions with 3D structures
-ds (--dwnld_strts)		Download 3D structures
-df (--dwnld_fasta)		Download FASTA files
-m (--method)			Method used to obtain structure [Default = All] (i.e., X-ray, NMR)
-o (--outdir)			Output directory for downloading files [Default = ./UNIPROT_SCRAP_RESULTS]
-c (--custom)			Custom UniProt search

"""

from sys import exit,argv

# if (len(argv) == 1):
#     print(f"{usage}")
#     exit()

import re
import os
import argparse
from os import system, mkdir, path
from time import strftime, localtime, time


pipeline_location = os.path.dirname(argv[0])

## Setup GetOptions
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-g","--go_keyword")
parser.add_argument("-v","--verified_only",action='store_true')
parser.add_argument("-s","--structures",action='store_true')
parser.add_argument("-ds","--download_structures",action='store_true')
parser.add_argument("-df","--download_fasta",action='store_true')
parser.add_argument("-f","--fasta",action='store_true')
parser.add_argument("-m","--method",nargs='+')
parser.add_argument("-o","--outdir",default="./UNIPROT_SCRAP_RESULTS")
parser.add_argument("-c","--custom")

args = parser.parse_args()

## Acquire argparse options
go_keyword = args.go_keyword
reviewed = args.verified_only
structures = args.structures
download_structures = args.download_structures
download_fasta = args.download_fasta
methods = args.method
outdir = args.outdir
custom = args.custom
fastadir = outdir + "/FASTA"
pdbdir = outdir + "/PDBs"

LOG = open(f"{outdir}/search.log", "w")

LOG.write(f"{argv[0]}\n")

## Prepare keywords for searching
keywords = []

if(go_keyword):
	keywords.append(f"(go:{go_keyword})")
if(structures):
	keywords.append("(structure_3d:true)")
if(reviewed):
	keywords.append("(reviewed:true)")

## If a custom search phrase is provided, overwrite keywords
if(custom):
	keywords = custom

if not (os.path.exists(outdir)):
	os.mkdir(outdir)
	os.mkdir(fastadir)

if not (os.path.exists(fastadir)):
	os.mkdir(fastadir)

if not (os.path.exists(pdbdir)):
	os.mkdir(pdbdir)

## Loading all the necessary packages for web scraping
from selenium import webdriver 
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options

start = strftime("%Y-%m-%d, %H:%M:%S",localtime())

## Creating web scraper
options = Options()
options.headless = True
driver = webdriver.Firefox(options=options,service_log_path=os.path.devnull)

url = "https://www.uniprot.org/uniprotkb?query="

url = url + keywords[0]

for ops in keywords[1:]:
	url = url + "AND" + ops

# LOG.write(f"Source from {url}\n")
# LOG.write(f"Keywords\n\t{keywords}\n")

print(f"\nConnecting to {url}")

## Connect to UniProt
driver.get(url)

####################################################
##### All waiting times are set to 10 mins max #####
####################################################

## Remove annoying display preference option provided by new UniProt site if it exists
try:
	element = WebDriverWait(driver,600).until(EC.presence_of_element_located((By.CLASS_NAME,"u_JGw")))
	for item in element.find_elements_by_css_selector("span"):
		if item.text == "Table":
			item.click()
			break
	element.find_element_by_css_selector("button").click()

except:
	next

section = driver.find_element_by_class_name("sidebar-layout__content")
divs = section.find_elements_by_css_selector("div")
buttons = divs[[i.get_attribute("class") for i in divs].index("button-group GUgbz")].find_elements_by_css_selector("button")
buttons[[i.text for i in buttons].index("Download")].click()

## Waiting for sidebar to load in

download_panel = ""

while True:
	try:
		download_panel = driver.find_element_by_css_selector("aside")
		break
	except:
		next

panel_content = download_panel.find_element_by_class_name("sliding-panel__content")
fields = panel_content.find_elements_by_css_selector("fieldset")

options = fields[0].find_elements_by_css_selector("option")
options[[i.text for i in options].index("List")].click()

options = fields[1].find_elements_by_css_selector("label")
options[[i.text for i in options].index("No")].click()

sections = panel_content.find_elements_by_css_selector("section")
download = sections[[i.get_attribute("class") for i in sections].index("button-group sliding-panel__button-row rUH91 dx6Wo")].find_element_by_css_selector("a").get_attribute("href")

## Downloading list of accession that match the searched parameters
system(f"wget \"{download}\" -O {outdir}/accessions.list 1>/dev/null 2>error.log")


## Parse list of accessions
accession_numbers = []

FILE = open(f"{outdir}/accessions.list","r")
for accession in FILE:
	accession = accession.replace("\n","")
	accession_numbers.append([accession,f"https://www.uniprot.org/uniprotkb/{accession}/entry"])
FILE.close()

## Begin trolling accession pages for metadata and links

metadata = {}
Downloads = {}

for accession,page in accession_numbers:

	## Store ["METADATA",[STRUCTURE_LINKS]]
	info = ["",[]]

	## Jump to Names & Taxonomy section to get metadata
	driver.get(f"{page}#names_and_taxonomy")

	names_info = driver.find_element_by_id("names_and_taxonomy")
	name_content = ""

	## Things load slow, so just keep trying

	while True:
		try:
			name_content = names_info.find_element_by_class_name("card__content")
			break
		
		except:
			next


	list_content = name_content.find_elements_by_css_selector("div.card__content > ul")

	protein_index = [i.get_attribute("data-article-id") for i in name_content.find_elements_by_css_selector("h3")].index("protein_names")
	name_index = [i.get_attribute("data-article-id") for i in name_content.find_elements_by_css_selector("h3")].index("organism-name")

	prot_name = list_content[protein_index].find_element_by_css_selector("strong").text

	for alt in list_content[protein_index].find_element_by_css_selector("strong").find_elements_by_css_selector("button"):
		prot_name = prot_name.replace(alt.text,"")

	prot_name = prot_name.replace("\n","")

	org_name = list_content[name_index].find_element_by_css_selector("a").text
	
	info[0] = org_name
	metadata[accession] = [org_name,prot_name,[],[]]

	print(f"Visiting page for {accession}\t{prot_name}")

	results = ""

	while True:

		time.sleep(1)

		try:

			driver.get(f"{page}#structure")
			structure_card = driver.find_element_by_id("structure")
			structure_display = structure_card.find_element_by_class_name("card__content")
			table_structure = structure_display.find_element_by_css_selector("protvista-uniprot-structure")
			table = table_structure.find_element_by_css_selector("table > tbody")
			results  = table.find_elements_by_css_selector("tr")
			break

		except:
			next

	for result in results:

		method = result.find_element_by_css_selector("td:nth-child(4)").text
		if method in methods:
			download_link = result.find_element_by_css_selector("td:nth-child(8) > a").get_attribute("href")
			info[1].append(download_link)

	Downloads[accession] = info

## Close the scraper, we are done surfing
driver.close()

downloaded = []
sources = ""

SOURCES = open(f"{outdir}/download.log", "w")
METADATA = open(f"{outdir}/metadata.log","w")

## Iterate over all UniProt accession found

print(f"Downloading data for {len(Downloads.keys())} accession!")

for accession in Downloads:

	## Split the data into workable chunks
	org_name,structure_links = Downloads[accession]

	print(f"Getting data for {accession}")

	if download_fasta:

		SOURCES.write(f"{accession}\t{org_name}\n")
		SOURCES.write("FASTA:\n")
		METADATA.write(f"{accession}\t{metadata[accession][0]}\t{metadata[accession][1]}\n")
		METADATA.write(f"\tFASTA\n")

		link = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

		print(f"Acquiring FASTA file {link}")

		fasta_name = f"{accession}.fasta"
		
		SOURCES.write(f"\t{link}\n")
		METADATA.write(f"\t\t{fasta_name}\n")

		## Download the FASTA file if not already
		if not path.exists(f"{fastadir}/{fasta_name}"):
			system(f"wget {link} -O {fastadir}/{fasta_name} 1>/dev/null 2>>{outdir}/download.errors")

	if download_structures:
		SOURCES.write("STRUCTURE:\n")
		METADATA.write(f"\tSTRUCTURE\n")

		for link in structure_links:
			
			## PDB prefix
			pdb_prefix = ""

			## Determine PDB source\
			db = ""
			
			### AlphaFold
			if(re.search("AF-(\w+)-\w\d-model",link)):
				pdb_prefix = re.search("AF-(\w+)-\w\d-model",link).groups(0)[0]
				db = "AF"
			
			### RCSB
			if(re.search("pdb(\w+)\.ent",link)):
				pdb_prefix = re.search("pdb(\w+)\.ent",link).groups(0)[0]
				db = "RCSB"

			SOURCES.write(f"\t{link}\n")

			## Download the RCSB PDB file if it hasn't been downloaded/split
			if not path.exists(f"{pdbdir}/{pdb_prefix}.pdb") and not path.isdir(f"{pdbdir}/{pdb_prefix}"):
				system(f"wget {link} -O {pdbdir}/{pdb_prefix}.pdb 1>/dev/null 2>>{outdir}/download.errors")

			## Split RCSB PDB file by chain and perform BLAST search to identify wanted chains
			if(db == "RCSB"):

				## Split the PDB file into different chains if not yet done
				if not path.exists(f"{pdbdir}/{pdb_prefix}") and path.exists(f"{pdbdir}/{pdb_prefix}.pdb"):
					system(f"{pipeline_location}/split_PDB.pl -p {pdbdir}/{pdb_prefix}.pdb -o {pdbdir}/{pdb_prefix}")

				## Extract the FASTA sequences from PDB files
				if not path.exists(f"{pdbdir}/{pdb_prefix}/FASTA/ALL.fasta") and any(File.endswith(".pdb") for File in os.listdir(f"{pdbdir}/{pdb_prefix}")):
					
					## Split the PDBS
					system(f"{pipeline_location}/extract_pdb_sequence.pl -p {pdbdir}/{pdb_prefix}/*.pdb -o {pdbdir}/{pdb_prefix}/FASTA")

					## Concatenate the FASTA files
					system(f"cat {pdbdir}/{pdb_prefix}/FASTA/*.fasta > {pdbdir}/{pdb_prefix}/FASTA/ALL.fasta")

				if path.exists(f"{pdbdir}/{pdb_prefix}/FASTA/ALL.fasta"):
					## If the sequences have been extracted, perform a BLAST search
					if path.getsize(f"{pdbdir}/{pdb_prefix}/FASTA/ALL.fasta") > 0:

						## Make BLAST db
						system(f"diamond makedb --in {pdbdir}/{pdb_prefix}/FASTA/ALL.fasta --db {pdbdir}/{pdb_prefix}/FASTA/3Database &>/dev/null")

						## Perform BLAST search
						system(f"diamond blastp -q {fastadir}/{fasta_name} --db {pdbdir}/{pdb_prefix}/FASTA/3Database --out {pdbdir}/{pdb_prefix}/FASTA/results.diamond.6 &>/dev/null")

						## Parse BLAST file for best results
						with open(f"{pdbdir}/{pdb_prefix}/FASTA/results.diamond.6","r") as BLAST:

							for line in BLAST:
								
								file = line.split("\t")[1]
								
								## Move the best PDB chain from its parent into home
								system(f"cp {pdbdir}/{pdb_prefix}/{file}.pdb {pdbdir}/{file}.pdb")

								## Add metadata to log
								if f"{file}.pdb" not in metadata[accession][3]:
									metadata[accession][3].append(f"{file}.pdb")
									METADATA.write(f"\t\t{file}.pdb\n")
								break
			
			elif f"{file}.pdb" not in metadata[accession][3]:
				metadata[accession][3].append(f"{file}.pdb")
				METADATA.write(f"\t\t{file}.pdb\n")

if download_structures:
	## Remove unneeded PDBs to prevent unwanted hits being returned
	for item in os.listdir(f"{pdbdir}/"):
		if (os.path.isdir(f"{pdbdir}/{item}")):
			system(f"rm -r {pdbdir}/{item}")

stop = time.strftime("%Y-%m-%d, %H:%M:%S",time.localtime())

LOG.write(f"Start\t{start}\nStop\t{stop}\n")

LOG.close()
METADATA.close()
SOURCES.close()

## Create an archive for later use
system(f"tar -czf {outdir}.tar.gz {outdir}")