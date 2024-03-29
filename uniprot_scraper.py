#!/usr/bin/python

name = "uniprot_scraper.py"
version = "1.8.0"
updated = "2022-12-04"

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
		  -c '(go:toll)AND(reviewed:true)AND(organism_id:7227)'


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

def die(string):
	exit(string)

import struct
from sys import exit,argv

if (len(argv) == 1):
	die(f"{usage}")

import re
import argparse
from os import system, mkdir, path, listdir
from time import sleep
from datetime import datetime

start_time = datetime.today()

pipeline_location = path.dirname(argv[0])

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

## Setup working directory
if not (path.exists(outdir)):
	mkdir(outdir)
	mkdir(fastadir)

if not (path.exists(fastadir)):
	mkdir(fastadir)

if not (path.exists(pdbdir)):
	mkdir(pdbdir)

## Loading all the necessary packages for web scraping
from selenium import webdriver 
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options

###################################################################################################
## Begin logging
###################################################################################################

## OPS stores basic operational information
OPS = open(f"{outdir}/uniprot_scraper.log","w")
OPS.write(f">NAME\n  {name}\n\n")
OPS.write(f">VERSION\n  {version}\n\n")

previous_downloads = {}
accession = ""
data_key = ""
if path.isfile(f"{outdir}/metadata.log"):
	META = open(f"{outdir}/metadata.log","r")
	for line in META:
		line = line.replace("\n","")
		if line:
			if line[0] == ">":
				accession = line.replace(">","",1)
			elif line[0:2] == "\t\t":
				if accession in previous_downloads.keys():
					if data_key in previous_downloads[accession].keys():
						previous_downloads[accession][data_key] = previous_downloads[accession][data_key] + "\n\t\t" + line.replace("\t\t","",1)
					else:
						previous_downloads[accession][data_key] = line.replace("\t\t","",1)
				else:
					previous_downloads[accession] = {}
					previous_downloads[accession][data_key] = line.replace("\t\t","",1)
			elif line[0] == "\t":
				data_key = line.replace("\t","")
	META.close()


## META stores the accession numbers, the organism names, the FASTA and 3D Structure links, and 3D structure features
META = open(f"{outdir}/metadata.log","w")

###################################################################################################
## Setting up webscaper object
###################################################################################################

## Create options for the scraper (really only used to make it headless)
options = Options()
options.headless = True

## Create the scraper object
driver = webdriver.Firefox(options=options,service_log_path=path.devnull)

## Setting up the results url using keywords
url = "https://www.uniprot.org/uniprotkb?query="

## Prepare keywords for searching
keywords = []

if(go_keyword):
	keywords.append(f"(go:{go_keyword})")
if(structures):
	keywords.append("(structure_3d:true)")
if(reviewed):
	keywords.append("(reviewed:true)")

## If a custom search phrase is provided, overwrite keywords

OPS.write(f">KEYWORDS\n")

if(custom):
	url = url + custom
	OPS.write(f"  {custom}\n")
else:
	url = url + keywords[0]
	OPS.write(f"  {keywords[0]}\n")
	for ops in keywords[1:]:
		url = url + "AND" + ops
		OPS.write(f"  {ops}\n")
OPS.write(f"\n")

OPS.write(f">DOWNLOAD_FASTA\n")
if(download_fasta):
	OPS.write(f"  TRUE\n")
else:
	OPS.write(f"  FALSE\n")
OPS.write(f"\n")

OPS.write(f">DOWNLOAD_STRUCTURE\n")
if(download_structures):
	OPS.write(f"  TRUE\n\n")
	OPS.write(f"  >METHODS\n")
	if(methods):
		for method in methods:
			OPS.write(f"    {method}\n")
	else:
		OPS.write(f"    ALL\n")
else:
	OPS.write(f"  FALSE\n")
OPS.write(f"\n")

OPS.write(f">SEARCH_URL\n  {url}\n\n")
OPS.write(f">LAUNCHED\n  {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

print(f"\nConnecting to {url}")

###################################################################################################
## Acquire accessions list for given keywords
###################################################################################################

### When the following structure is seen, we are waiting for things to load in
## try:
##		some command
##		break
## except Exception:
##		next

try:
	## Connect to UniProt results page
	driver.get(url)

	## Remove annoying display preference option provided by new UniProt site if it exists
	try:
		element = WebDriverWait(driver,600).until(EC.presence_of_element_located((By.CLASS_NAME,"u_JGw")))
		for item in element.find_elements_by_css_selector("span"):
			if item.text == "Table":
				item.click()
				break
		element.find_element_by_css_selector("button").click()
	except Exception:
		next

	## Open the download options bar
	section = driver.find_element_by_class_name("sidebar-layout__content")
	divs = section.find_elements_by_css_selector("div")
	buttons = divs[[i.get_attribute("class") for i in divs].index("button-group GUgbz")].find_elements_by_css_selector("button")
	buttons[[i.text for i in buttons].index("Download")].click()

	## Waiting for sidebar to load in
	download_panel = ""
	while True:
		sleep(1)
		try:
			download_panel = driver.find_element_by_css_selector("aside")
			break
		except Exception:
			next

	panel_content = download_panel.find_element_by_class_name("sliding-panel__content")
	fields = panel_content.find_elements_by_css_selector("fieldset")

	## Changing download to accession only List
	while True:
		sleep(1)
		try:
			options = fields[0].find_elements_by_css_selector("option")
			options[[i.text for i in options].index("List")].click()
			break
		except Exception:
			next

	## Changing download to decompressed type (not a large download, doesn't really matter)
	while True:
		sleep(1)
		try:
			options = fields[1].find_elements_by_css_selector("label")
			options[[i.text for i in options].index("No")].click()
			break
		except Exception:
			next

	sections = panel_content.find_elements_by_css_selector("section")
	download = sections[[i.get_attribute("class") for i in sections].index("button-group sliding-panel__button-row rUH91 dx6Wo")].find_element_by_css_selector("a").get_attribute("href")

	## Downloading list of accession that match the searched parameters
	system(f"wget \"{download}\" -O {outdir}/accessions.list 1>/dev/null 2>>{outdir}/download.error")

	###################################################################################################
	## Scraping accession pages for metadata
	###################################################################################################

	## Parse list of accessions
	accession_numbers = []

	FILE = open(f"{outdir}/accessions.list","r")
	for accession in FILE:
		accession = accession.replace("\n","")
		accession_numbers.append([accession,f"https://www.uniprot.org/uniprotkb/{accession}/entry"])
	FILE.close()

	## Begin trolling accession pages for metadata and links

	print("\nAcquiring metadata:\n")

	scrap_results = {}
	for count,(accession,page) in enumerate(sorted(accession_numbers,key=lambda x: x[0])):

		print(f"[{count+1:0>{len(str(len(accession_numbers)))}}/{len(accession_numbers)}]\t{accession}")

		if accession in previous_downloads.keys():
			print(f"\tData previously acquired... Skipping...")
			META.write(f">{accession}\n")
			META.write(f"\tPROTEIN_NAME\n")
			META.write(f"\t\t{previous_downloads[accession]['PROTEIN_NAME']}\n")
			META.write(f"\tORGANISM_NAME\n")
			META.write(f"\t\t{previous_downloads[accession]['ORGANISM_NAME']}\n")
			META.write(f"\tFEATURES\n")
			META.write(f"\t\t{previous_downloads[accession]['FEATURES']}\n")
			META.write(f"\tSTRUCTURES\n")
			
			if 'STRUCTURES' in previous_downloads[accession].keys():
				META.write(f"\t\t{previous_downloads[accession]['STRUCTURES']}\n")
			else:
				META.write(f"\n")
			
			continue

		## acession = [FASTA link,features,[structures]]
		scrap_results[accession] = [f"https://rest.uniprot.org/uniprotkb/{accession}.fasta","",[]]

		##################################################
		## Name and Taxonomy Metadata
		##################################################
		
		## Jump to Names & Taxonomy section to get metadata
		driver.get(f"{page}#names_and_taxonomy")

		names_info = driver.find_element_by_id("names_and_taxonomy")
		name_content = ""

		## Wait for Names & Taxonomy content to load
		while True:
			sleep(1)
			try:
				name_content = names_info.find_element_by_class_name("card__content")
				break
			
			except Exception:
				next


		header_content = name_content.find_elements_by_css_selector("h3")
		list_content = name_content.find_elements_by_class_name("info-list")

		## Index of protein info in list_content
		protein_index = [i.text for i in header_content].index("Protein names")
		print(f"Protein Index: {protein_index}")
		## Index of name info in list_content
		name_index = [i.text for i in header_content].index("Organism names")
		print(f"Name Index: {name_index}")

		## Get protein name
		prot_name = list_content[protein_index].find_element_by_css_selector("strong").text
		prot_name = prot_name.replace("\n","")
		## Remove "curate" and "publication" text from the protein name as it is not relevant!
		for alt in list_content[protein_index].find_elements_by_css_selector("button"):
			prot_name = prot_name.replace(alt.text,"")

		## Get organism name
		org_name = list_content[name_index].find_element_by_css_selector("div > div.decorated-list-item__content").text.replace("\n","")
		## Remove "curate" and "publication" text from the protein name as it is not relevant!
		for alt in list_content[name_index].find_elements_by_css_selector("button"):
			org_name = org_name.replace(alt.text,"")

		META.write(f">{accession}\n")
		META.write(f"\tPROTEIN_NAME\n\t\t{prot_name}\n")
		META.write(f"\tORGANISM_NAME\n\t\t{org_name}\n")
		print(f"\t{accession} ({prot_name})\n")
		if (not path.exists(f"{fastadir}/{accession}.fasta")):
			print(f"\t\tDownloading FASTA file for {accession}")
			system(f"wget https://rest.uniprot.org/uniprotkb/{accession}.fasta -O {fastadir}/{accession}.fasta 1>/dev/null 2>>{outdir}/download.error")
		else:
			print(f"\t\tFASTA previously downloaded for {accession}...Skipping...")

		##################################################
		## Structure Metadata
		##################################################

		## Jump to Structure section to get metadata
		driver.get(f"{page}#structure")

		## Get structure results
		results = False
		while True:
			sleep(1)
			try:
				structure_card = driver.find_element_by_id("structure")
				structure_display = structure_card.find_element_by_class_name("card__content")
				table_structure = structure_display.find_element_by_css_selector("protvista-uniprot-structure")
				table = table_structure.find_element_by_css_selector("table > tbody")
				results  = table.find_elements_by_css_selector("tr")
				break
			except Exception:
				next

			try:
				structure_card = driver.find_element_by_id("structure")
				structure_display = structure_card.find_element_by_class_name("card__content")
				table_structure = structure_display.find_element_by_css_selector("protvista-uniprot-structure")
				table_structure.find_element_by_class_name("protvista-no-results")
				break
			except Exception:
				next
		
		## Get feature results
		attributes = False
		while True:
			sleep(1)
			try:
				structure_card = driver.find_element_by_id("structure")
				structure_display = structure_card.find_element_by_class_name("card__content")
				headers = [i.text for i in structure_display.find_elements_by_css_selector("h3")]
				if "Features" in headers:
					attributes = True
				break
			except Exception:
				next

		META.write(f"\tFEATURES\n")
		
		struct_atts = {}
		
		if attributes:

			feature_table = structure_display.find_elements_by_css_selector("protvista-datatable")[-1]
			table_object = feature_table.find_element_by_css_selector("table > tbody")
			feature_results = table_object.find_elements_by_css_selector("tr")
			
			for result in feature_results:
				if "hidden" not in result.get_attribute("class"):
					feat = result.find_element_by_css_selector("td:nth-child(2)").text
					if feat not in struct_atts.keys():
						struct_atts[feat] = 0
					struct_atts[feat] += 1
			
			for key in sorted(struct_atts.keys()):
				META.write(f"\t\t{key.upper()}:{struct_atts[key]}\n")
			
			scrap_results[accession][1] = struct_atts

		else:
			
			META.write(f"\t\tNone Available\n")

		META.write(f"\tSTRUCTURES\n")

		structure_data = []

		if results:

			for result in results:

				## Get the structure type
				method = result.find_element_by_css_selector("td:nth-child(4)").text

				## Set the chain for non-predicted structures, and set the path to the RCSB PDB rather than the PDBe (don't know if it makes a difference other than who gets credit)
				if method != "Predicted":
					pdb = result.find_element_by_css_selector("td:nth-child(3)").text
					chain = result.find_element_by_css_selector("td:nth-child(6)").text

					if "/" in chain:
						chain = chain.split("/")[0]

					download_link = f"https://files.rcsb.org/download/{pdb}.pdb"
				else:
					pdb = accession
					chain = "-"
					download_link = result.find_element_by_css_selector("td:nth-child(9) > a").get_attribute("href")


				## Get the PDB
				structure_data.append([pdb,chain,method,download_link])
				if method in methods:
					META.write(f"\t\t{pdb}\t{chain}\t{method}\t{download_link}\n")

			def get_pdb(struct_link,pdb_code,method,chain):

				if ((not path.isdir(f"{pdbdir}/{pdb_code}")) and (not path.isfile(f"{pdbdir}/{pdb_code}.pdb.gz"))):
					system(f"wget {struct_link} -O {pdbdir}/{pdb_code}.pdb 1>/dev/null 2>>{outdir}/download.error")
					print(f"\t\tDownloading {struct_link}")
				else:
					print(f"\t\tSkipping {pdb_code}, already downloaded")

				if method != "Predicted":
					if (not path.isdir(f"{pdbdir}/{pdb_code}")):
						system(f"""
							{pipeline_location}/split_PDB.pl \\
							-p {pdbdir}/{pdb_code}.pdb \\
							-o {pdbdir}/{pdb_code}
						""")
					system(f"""
						if [ -f {pdbdir}/{pdb_code}/{pdb_code}_{chain}.pdb ]
						then
							if ! [ -f {pdbdir}/{pdb_code}_{chain}.pdb.gz ]
							then
								mv {pdbdir}/{pdb_code}/{pdb_code}_{chain}.pdb {pdbdir}
								gzip {pdbdir}/{pdb_code}_{chain}.pdb
							fi
						fi
					""")
				else:
					if not path.isfile(f"{pdbdir}/{pdb_code}.pdb.gz"):
						system(f"""
							gzip {pdbdir}/{pdb_code}.pdb
						""")
			
			for set in structure_data:

				pdb, chain, method, download_link = set

				if methods:
					if method in methods:
						get_pdb(download_link,pdb,method,chain)
				else:
					get_pdb(download_link,pdb,method,chain)

			print()
			
			scrap_results[accession][2] = structure_data

		else:

			META.write(f"\t\tNONE\n")

	print()

	## Close the scraper, we are done surfing
	driver.close()

	if download_structures:
		## Remove unneeded PDBs to prevent unwanted hits being returned
		for item in listdir(f"{pdbdir}/"):
			if (path.isdir(f"{pdbdir}/{item}")):
				system(f"rm -r {pdbdir}/{item}")

	stop_time = datetime.today()
	OPS.write(f">COMPLETED\n  {stop_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

	hours, rest = divmod((stop_time-start_time).total_seconds(),60*60)
	mins, secs = divmod(rest,60)

	OPS.write(f">RUNTIME\n  {int(hours)}:{int(mins)}:{round(secs,2)}\n")

	OPS.close()
	META.close()

except Exception:
	driver.close()
	exit("It appears that UniProt has changed its HTML... Please inform QueGO developers to resolve error")