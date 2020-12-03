#!/usr/bin/python

# This works identically with either Python2 or Python3

# sends output to stdout.

# This code was originally written to parse PanGIA output, and now
# appears to be working correctly for parsing DIAMOND output after some minor changes and additions 101717 Geoffrey House.

# It looks up the full taxonomic classification for each taxon ID in the input file (using module taxonomy.py) and returns
# that, which is then appended (tab-delimited) after the number of reads assigned to each taxonomic group (NOT rollup counts).
# Any taxonomic group with 0 assigned reads (regardless of rollup) are dropped from the tab_tree output.

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = ""
__date__      = ""
__copyright__ = ""

import argparse as ap
import sys
import os
from re import search, findall
import taxonomy as t

def parse_params():
	p = ap.ArgumentParser(prog='convert_diamond2tabTree.py', description="""Convert DIAMOND output (list file) to tab_tree format.""")

	p.add_argument('-i', '--inputFile', type = str, help="""Input taxonomy list file to convert to tab_tree format.""")

	p.add_argument('-dp', '--dbPath',
			metavar='[PATH]', type=str, default=None,
					help="""Path of databases. (i.e. parent dir to nodes.dmp and names.dmp like '/panfs/biopan01/refdb/usrdb/diamond_db/taxonomy'). 
					If dbPath isn't specified but a path is provided in "--database" option, this path of database will also be used in dbPath. 
					Otherwise, the program will search "database/" in program directory.
					[default: database/]""")

	args_parsed = p.parse_args()

	if not args_parsed.dbPath:
		if args_parsed.database and "/" in args_parsed.database:
			db_dir = search( '^(.*?)[^\/]+$', args_parsed.database )
			args_parsed.dbPath = db_dir.group(1)
		else:
			bin_dir = os.path.dirname(os.path.realpath(__file__))
			args_parsed.dbPath = bin_dir + "/database"

	return args_parsed

if __name__ == '__main__':
	argvs = parse_params()
	t.loadTaxonomy( argvs.dbPath )
	
	# Originally read input from stdin
	#for line in sys.stdin:

	with open(argvs.inputFile, 'r') as inFile:	

		for line in inFile:
			# skip header
			if line.startswith("LEVEL"):
				continue
			temp = line.split('\t')
		
			'''I don't want this, I need it to be at all taxon. levels'''
			# only process species level
			#if temp[0] != "species":
				#continue
	

			# New - skip entries (lines) that have 0 assigned reads directly to them
			# (regardless of whether they have any rollup reads assigned to them or not).
			# This makes the output tab_tree consistent with the formatting from kraken-mini.

			if int(temp[3]) == 0:
				continue	
	
			'''
			# only process unfiltered
			if not argvs.displayAll and temp[15]:
				continue
			# display pathogen only
			if not argvs.pathogenOnly and not temp[10]:
				continue
			'''
        
			# Original
			#lineage = t.taxid2lineage(temp[2])
        
			'''To parse the correct tax ID column from the diamond list output'''
			lineage = t.taxid2lineage(temp[4])
        
			# This puts the number of assigned reads for each taxonomic group (NOT the rollup) as the 
			# first column, and then the full taxonomy of the organism as a tab delimited list that follows.
			# This splits the elements of lineage on '|' and then joins them using '\t' instead.
			print( "%s\t%s" %
				( temp[3],
				  '\t'.join( lineage.split('|') )
			))
