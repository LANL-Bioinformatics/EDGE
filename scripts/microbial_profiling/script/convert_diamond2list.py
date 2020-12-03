#!/usr/bin/env python

# This is Paul's original script called 'diamond_class_reporter.py'. The name has been changed
# in EDGE to match the naming for the other tools.
		
# Modified to work with either Python3 or Python2 by converting the counts to floats so
# the relative abundance is calculated correctly. This is working correctly 101617 Geoffrey House
	
__author__  = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__ = "0.1"
__date__    = "2017/08/26"
__license__ = "GPLv3"

import argparse as ap, textwrap as tw
import sys, os, time, subprocess
import taxonomy as t

def parse_params(ver):
	p = ap.ArgumentParser(prog='diamond_class_report.py', description="""Diamond classification reporter %s""" % ver)

	p.add_argument('-i', '--input', 
			metavar='[FILE]', type=ap.FileType('r'), default=sys.stdin,
			        help="Specify a file of Diamond read classification results. [default: STDIN]")

	p.add_argument('-tp', '--taxaPath',
			metavar='[PATH]', type=str, default='./taxonomy',
	                help="""Path of taxonomy information (taxonomy.tsv or names.dmp/nodes.dmp). [default: ./taxonomy]""")

	p.add_argument( '-o','--output', metavar='[FILE]', type=ap.FileType('w'), default=sys.stdout,
					help="Output file [default: STDOUT]")
	
	args_parsed = p.parse_args()
	return args_parsed

"""
The function "parsing" will take Diamond classification result line by line
"""
def parsing():
	cnt = 0

	# parsing results
	for line in argvs.input:
		cnt += 1
		# Example line: D85TC5M1:229:C2JUAACXX:1:1101:2077:2173 32008   7.5e-05
		# Only process complete lines
		temp = line.split('\t')

		(readid, taxid, evalue) = ("","","")
		if len(temp) == 3:
			(readid, taxid, evalue) = temp
		else:
			continue

		# reference doesn't have a mapped taxid or taxid doesn't have tax info
		if taxid == '0': continue

		lineage = t.taxid2lineageDICT( taxid, 1, 1 )
		assigned_rank = False

		# Directory "res_rollup" maintains a directory of "ASGN"|"ROLL"|"NAME" for all major ranks and mapped taxonomy id.
		for rank in sorted( major_ranks, key=major_ranks.__getitem__, reverse=True ):
			if not rank in lineage:
				continue
			else:
				tid = lineage[rank]['taxid']
				name = lineage[rank]['name']
			
				# count the read to the latest rank
				if not assigned_rank:
					assigned_rank = True
					if tid in res_rollup[rank]:
						res_rollup[rank][tid]["ASGN"] += 1
					else:
						res_rollup[rank][tid]["ASGN"] = 1
						res_rollup[rank][tid]["ROLL"] = 0
						res_rollup[rank][tid]["NAME"] = name

				# count for rolling up taxonomies
				if tid in res_rollup[rank]:
					res_rollup[rank][tid]["ROLL"] += 1
				else:
					res_rollup[rank][tid]["ROLL"] = 1
					res_rollup[rank][tid]["ASGN"] = 0
					res_rollup[rank][tid]["NAME"] = name

		if ( cnt % 1000 == 0):
			sys.stderr.write( "[INFO] Processing %s read classifications...\r"%cnt )
	
	return cnt

def write_report(f, tol_read_count):
	# print header
	f.write( "LEVEL\tNAME\tROLLUP\tASSIGN\tTAXID\tREL_ABU\n" )
	# print report
	for rank in sorted( major_ranks, key=major_ranks.__getitem__ ):
		# sort by rollup read count
		for (tid,v) in sorted( res_rollup[rank].items(), key=lambda r: r[1]["ROLL"], reverse=True):
			f.write( "%s\t%s\t%s\t%s\t%s\t%.4f\n" % (
				rank,
				res_rollup[rank][tid]["NAME"],
				res_rollup[rank][tid]["ROLL"],
				res_rollup[rank][tid]["ASGN"],
				tid,
				float(res_rollup[rank][tid]["ROLL"])/float(tol_read_count)
			))

if __name__ == '__main__':
	argvs = parse_params( __version__ )
	#load taxonomy
	sys.stderr.write( "[INFO] Loading taxonomy...\n" )
	t.loadTaxonomy( argvs.taxaPath )
	sys.stderr.write( "[INFO] Done.\n" )
	
	major_ranks    = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}
	res_rollup     = t._autoVivification()
	tol_read_count = 0

	sys.stderr.write( "[INFO] Start processing read classifications...\n" )
	tol_read_count = parsing()
	sys.stderr.write( "[INFO] Done. %s read classifications are processed.\n"%tol_read_count )
	write_report( argvs.output, tol_read_count )
	sys.stderr.write( "[INFO] Done writing report.\n" )
