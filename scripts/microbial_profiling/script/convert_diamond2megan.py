#!/usr/bin/python

# This is a very basic script to convert DIAMOND's list file to a file for megan (if necessary in EDGE),
# which has two columns (with no header):
# 	1) Taxon name
#	2) number of assigned reads (NOT rollup) for each taxa. Taxa with 0 assigned reads
#	(regardless of rollup) are dropped from the output

# This script prints its output to stdout for easy redirection to a file.

__author__    = "Geoffrey House, Bioscience Division, Los Alamos National Laboratory"
__version__   = ""
__date__      = ""
__copyright__ = "" 


import argparse

parser = argparse.ArgumentParser(description = ("Converts DIAMOND list file to a file for megan."))

parser.add_argument('-i', '--inputFile', type = str, help = ("Input DIAMOND taxonomy list file to convert."))

args = parser.parse_args()

inFileName = args.inputFile

with open(inFileName, 'r') as inFile:
	for line in inFile:
		
		splitLine = line.strip().split('\t')

		# Continue for the header
		if splitLine[0] == "LEVEL":
			continue
		
		# Continue if the number of assigned reads is 0
		if int(splitLine[3]) == 0:
			continue

		# Otherwise process the line to be in krona format and output to stdout.	
		print(splitLine[1] + '\t' + splitLine[3]) 
		
		
