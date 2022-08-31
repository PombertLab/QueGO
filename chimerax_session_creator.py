#!/usr/bin/python
from chimerax.core.commands import run
from sys import argv
import argparse
import re
import os

name = 'chimerax_session_creator.py'
version = '0.4.0'
updated = '2022-08-31'

usage = f'''
NAME		{name}
VERSION		{version}
UPDATED		{updated}

SYNOPSIS	This script is used to align a reference .pdb to a predicted .pdb,
		changes the predicted .pdb color, hides all atoms, shows only matched
		chains, and saves the result as a ChimeraX session, .cxs file. This version
		is tested and functional as of ChimeraX 1.3.1.

COMMAND		{name} \\
			-p ...preference \\
			-r ...reference

OPTIONS

-p (--provided_pdb)		Predicted .pdb file
-m (--match_pdb)		RCSB .pdb file
-m (--match)	RCSB match name
-c (--chain)	RCSB matched chain
-o (--outdir)	Output directory for .cxs files [Default: ./3D_Visualizations]
'''

if len(argv) < 2:
	print(f"\n\n{usage}\n\n")
	exit()

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument('-p','--provided_pdb',type=str,required=True)
parser.add_argument('-m','--match_pdb',type=str,required=True)
parser.add_argument('-o','--outdir',type=str,default="./3D_Visualizations")
parser.add_argument('--nogui')

args = parser.parse_args()
provided_pdb = args.provided_pdb
match_pdb = args.match_pdb
if(args.outdir):
	outdir = args.outdir

provided_pdb_name = (os.path.splitext(os.path.basename(provided_pdb))[0]).replace(".temp","")
match_pdb_name = (os.path.splitext(os.path.basename(match_pdb))[0]).replace(".temp","")

## Load pdb files
model_provided_pdb = run(session,f"open {provided_pdb}")[0]
model_provided_pdb_name = (model_provided_pdb.id_string)

model_rcsb = run(session,f"open {match_pdb}")[0]
model_rcsb_name = (model_rcsb.id_string)


## Prepare file for display by hiding everything
run(session,"hide atoms")
run(session,"hide ribbons")

match = run(session,f"match #{model_provided_pdb_name} to #{model_rcsb_name}")

## Color reference structure a diferrent color
run(session,f"color #{model_rcsb_name} #00FFFF ribbons")

## Show only matching chains
run(session,f"show #{model_provided_pdb_name} ribbons")
run(session,f"show #{model_rcsb_name} ribbons")

## Orient the chain to view
run(session,"view")

## Save match as a new file
run(session,f"save {outdir}/{provided_pdb_name}_{match_pdb_name}.cxs format session")

quit()