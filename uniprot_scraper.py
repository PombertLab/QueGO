#!/usr/bin/python

name = "uniprot_scraper.py"
version = "0.5"
updated = "2022-04-25"

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
		  -c 'toll database:(type:pdb) AND reviewed:yes AND organism:"Drosophila melanogaster (Fruit fly) [7227]"'


OPTIONS
-g (--go_annotation)		Gene ontoloy keyword to search
-v (--verified_only)		Get downloads for genes that have been verified by UniProt
-s (--structures)		Find accessions with 3D structures
-ds (--dwnld_strts)		Download 3D structures
-df (--dwnld_fasta)		Download FASTA files
-m (--method)			Method used to obtain structure [Default = All] (i.e., X-ray, NMR)
-c (--custom)			Custom UniProt search
-o (--outdir)			Output directory for downloading files [Default = ./UNIPROT_SCRAP_RESULTS]

"""

from sys import exit,argv

if (len(argv) == 1):
    print(f"{usage}")
    exit()

import re
import os
import time
import argparse
from os import system, mkdir, path


pipeline_location = os.path.dirname(argv[0])

## Setup GetOptions
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-g","--go_annotation")
parser.add_argument("-v","--verified_only",action='store_true')
parser.add_argument("-s","--structures",action='store_true')
parser.add_argument("-ds","--download_structures",action='store_true')
parser.add_argument("-df","--download_fasta",action='store_true')
parser.add_argument("-f","--fasta",action='store_true')
parser.add_argument("-m","--method")
parser.add_argument("-c","--custom")
parser.add_argument("-o","--outdir",default="./UNIPROT_SCRAP_RESULTS")

args = parser.parse_args()

## Acquire argparse options
go_annotation = args.go_annotation
reviewed = args.verified_only
structures = args.structures
download_structures = args.download_structures
download_fasta = args.download_fasta
method = args.method
custom = args.custom
outdir = args.outdir
fastadir = outdir + "/FASTA"
pdbdir = outdir + "/PDBs"


LOG = open(f"{outdir}/search.log", "w")

LOG.write(f"{argv[0]}\n")

## Prepare key words for searching
keywords = ""

if(go_annotation):
	keywords = f"goa:({go_annotation})"
if(reviewed):
	keywords = f"{keywords} AND reviewed:yes"
if(method):
	keywords = f"{keywords} AND method:({method})"
if(structures):
	keywords = f"{keywords} AND database:(type:pdb)"

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

start = time.strftime("%Y-%m-%d, %H:%M:%S",time.localtime())

## Creating web scraper
options = Options()
options.headless = True
driver = webdriver.Firefox(options=options,service_log_path=os.path.devnull)

url = "https://www.uniprot.org/"

LOG.write(f"Source from {url}\n")
LOG.write(f"Keywords\n\t{keywords}\n")

print(f"\nConnecting to {url}")

## Connect to UniProt
driver.get(url)

####################################################
##### All waiting times are set to 10 mins max #####
####################################################

## Wait for page load
WebDriverWait(driver,600).until(EC.presence_of_element_located((By.ID, "query")))
## Find search bar
search_bar = driver.find_element_by_id("query")
## Enter query into search bar
search_bar.send_keys(keywords)
## Find the search button
search_button = driver.find_element_by_id("search-button")
## Click the search button
search_button.click()

print(f"Searching for keywords '{keywords}'")

## Wait for dropdown (tells us the results are in)
WebDriverWait(driver,600).until(EC.presence_of_element_located((By.CLASS_NAME,"limitdropDown")))
## Find the number of results drop down box
item_limit = driver.find_element_by_class_name("limitdropDown")
## Get the options
options = item_limit.find_elements_by_css_selector("option")
## Select the max result for SPEED
options[4].click()

print("Aquiring matches")

LOG.write("Accession Found\n")

accession_numbers = []

## Visit next page of results if there is one
loop = True
while loop == True:
	## Wait for results to show up
	WebDriverWait(driver,600).until(EC.presence_of_element_located((By.ID, "results")))
	## HTML hunting for accession numbers and corresponding links
	content_table = driver.find_element_by_id("results")
	table = content_table.find_element_by_css_selector("tbody")

	results = table.find_elements_by_css_selector("tr")

	## For each accession number, get its corresponding link
	for i in results:
		
		accession = i.get_attribute("id")
		
		print(f"Found Query\t{accession}")
		
		LOG.write(f"\t{accession}\n")

		href = i.find_element_by_class_name("entryID > a").get_attribute("href")
		accession_numbers.append([accession,href])
	
	## Try to click on next page else, start accession number surfing
	try:
		nextPage = driver.find_element_by_class_name("nextPageLink").get_attribute("href")
		driver.get(nextPage)
		print("Moving to next page of results")
	except:
		loop = False

print()

metadata = {}
Downloads = {}

for pages in accession_numbers:

	## Store ["METADAT",[FASTA_LINKS],[STRUCTURE_LINKS]]
	info = ["",[],[]]

	## Get accession number and go to its page
	accession = pages[0]
	driver.get(pages[1])

	## Wait for content to load
	WebDriverWait(driver,600).until(EC.presence_of_element_located((By.ID,"content-organism")))
	## HTML hunting for the name of the organism the query belongs to
	organism_title = driver.find_element_by_id("content-organism")
	organism_name = organism_title.find_element_by_tag_name("em").text
	info[0] = organism_name

	## Wait for content to load
	WebDriverWait(driver,600).until(EC.presence_of_element_located((By.ID,"content-protein")))
	## HTML hunting for the protein name
	protein_title = driver.find_element_by_id("content-protein")
	protein_name = protein_title.find_element_by_tag_name("h1").text

	print(f"Visiting page for {protein_name}")

	metadata[accession] = [organism_name,protein_name,[],[]]

	## If downloading FASTA files, go HTML hunting
	if(download_fasta):

		print(f"\tSearching for FASTA links for {accession} found in {organism_name}")

		try:

			## Wait to load content
			WebDriverWait(driver,600).until(EC.presence_of_element_located((By.ID,"sequences-section")))
			sequences_section = driver.find_element_by_id("sequences-section")
			download_button = sequences_section.find_elements_by_tag_name("a")
			
			## Search all the HTML for the FASTA link
			for a in download_button:
				## Get the first fasta link, and break out
				if(a.get_attribute("class") == "tooltipped icon icon-functional button inlineDisplayThis"):
					fasta_link = a.get_attribute("href")
					info[1].append(fasta_link)
					print(f"\t\tFASTA link\n\t\t\t{fasta_link}")
					break

		except:

			LOG.write(f"Unable to acquire fasta download link for {accession}\n")
			continue

	## If downloading PDB files, go HTML hunting
	if(download_structures):

		print(f"\tSearching for 3D structure links for {accession} found in {protein_name}")

		try:

			## Wait for content to load
			WebDriverWait(driver,600).until(EC.presence_of_element_located((By.TAG_NAME, "protvista-datatable")))
			## HTML hunting for PDB file links
			datatable = driver.find_element_by_tag_name("protvista-datatable")
			table = datatable.find_element_by_tag_name("table")
			tbody = datatable.find_element_by_tag_name("tbody")
			rows = tbody.find_elements_by_tag_name("tr")

			print("\t\tPDB Links")

			## Iterate over all the HTML garbage for the download links
			for row in rows:
				## Get the PDB link
				item = row.find_elements_by_tag_name("td")[8]
				download_link = item.find_element_by_tag_name("a").get_attribute("href")
				info[2].append(download_link)
				print(f"\t\t\t{download_link}\n")

		except:
			LOG.write(f"Unable to acquire structure download links for {accession}\n")
			continue

	if not download_fasta and not download_structures:
		print("Acquiring metadata only!")
	
	Downloads[accession] = info

## Close the scraper, we are done surfing
driver.close()

downloaded = []
sources = ""

SOURCES = open(f"{outdir}/download.log", "w")
METADATA = open(f"{outdir}/metadata.log","w")
## Iterate over all UniProt accession found
### For each accession there can be a FASTA file

print(f"Downloading data for {len(Downloads.keys())} accession!")

for accession in Downloads:

	## Split the data into workable chunks
	org_name,fasta_link,structure_links = Downloads[accession]

	SOURCES.write(f"{accession}\t{org_name}\n")
	SOURCES.write("FASTA:\n")
	METADATA.write(f"{accession}\t{metadata[accession][0]}\t{metadata[accession][1]}\n")
	METADATA.write(f"\tFASTA\n")

	print(f"Getting data for {accession}")
		
	## Iterate over all FASTA links (This is redundant because currently hardcoded to get (1) FASTA link)
	for link in fasta_link :

		print(f"Acquiring FASTA file {link}")

		fasta_name = re.search("\/(\w+\.fasta)$",link).groups(0)[0]
		
		## Download the FASTA file if not already
		if not path.exists(f"{fastadir}/{fasta_name}"):
			system(f"wget {link} -O {fastadir}/{fasta_name} 1>/dev/null 2>>{outdir}/download.errors")

		SOURCES.write(f"\t{link}\n")
		METADATA.write(f"\t\t{fasta_name}\n")
	
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
	
## Remove unneeded PDBs to prevent unwanted hits being returned
for item in os.listdir(f"{pdbdir}/"):
	if (os.path.isdir(f"{pdbdir}/{item}")):
		system(f"rm -r {pdbdir}/{item}")

stop = time.strftime("%Y-%m-%d, %H:%M:%S",time.localtime())

LOG.write(f"Start\t{start}\nStop\t{stop}\n")

with open(f"{outdir}/metadata.log","w") as METADATA:
	for key in metadata.keys():
		METADATA.writelines(f"{key}\t{metadata[key][0]}\t{metadata[key][1]}\n")
		METADATA.writelines(f"\tFASTA\n")
		for fasta in metadata[key][2]:
			METADATA.writelines(f"\t\t{fasta}\n")
		METADATA.writelines(f"\tSTRUCTURE\n")
		for struct in metadata[key][3]:
			METADATA.writelines(f"\t\t{struct}\n")