#!/usr/bin/env python
import sys
import re
import ConfigParser
import json
from Bio import SeqIO
import argparse


def getGFFInfo (rgiHitList, gffTableList, outGFFList, outCoordsList, idList):
	(orfID, contig, start, stop, orientation, cutOff, passEvalue, bestHitBitscore, bestHitAROName, bestHitARO, bestHitCategories, identities, aros, aroNames, modelType, snps, aroCategories, bitscores, predictedProtein, cardProtein, label, hitID) = rgiHitList

	geneRGI = orfID.split(" ")[0]
	for gffEntry in gffTableList:
		(gffSeqName, gffSource, gffFeature, gffStart, gffStop, gffScore, gffStrand, gffFrame, gffAttributes) = gffEntry
		gffAttributeList = gffAttributes.split(';')
		newAttributeList = []
		if geneRGI in gffAttributeList[0]:
			#print gffAttributeList[0]
			#print geneRGI, gffAttributeList
			for gffAttribute in gffAttributeList:
				if "ID=" in gffAttribute:
					newID = "ID=" + bestHitAROName
					if newID in idList:
						newID += ":" + str(idList.count(newID))
					newAttributeList.append(newID)
					idList.append("ID=" + bestHitAROName)
				elif "gene=" in gffAttribute:
					newAttributeList.append(gffAttribute + ":" + bestHitAROName)
				elif "inference=" in gffAttribute:
					newAttributeList.append("inference=RGI")
				elif "locus_tag=" in gffAttribute:
					newAttributeList.append(gffAttribute)
				elif "product=" in gffAttribute:
					newAttributeList.append(gffAttribute)
			gffFeature = "AR"
			gffSource = "RGI"
			newAttributeList.append("rgi_best_ARO=" + bestHitARO)
			newAttributeList.append("rgi_CARD_link=https://card.mcmaster.ca/aro/" + bestHitARO)
			newAttributeList.append("rgi_best_categories=" + bestHitCategories)
			newAttributeList.append("rgi_cutoff_type=" + cutOff)
			newAttributeList.append("rgi_model_type=" + modelType)
			newAttributeList.append("rgi_AROs=" + aros)
			newAttributeList.append("rgi_ARO_names=" + aroNames)
			newAttributeList.append("rgi_SNPs=" + snps)
			newAttributeList.append("rgi_ARO_categories=" + aroCategories)
			newAttributeList.append("rgi_best_bitscore=" + bestHitBitscore)
			newAttributeList.append("Card_rotein_sequence=" + cardProtein)
			newAttributeList.append("rgi_protein_sequence=" + predictedProtein)
			gffAttributes = ";".join(newAttributeList)
			gffRGIEntry = [gffSeqName, gffSource, gffFeature, gffStart, gffStop, gffScore, gffStrand, gffFrame,
			               gffAttributes]
			outGFFList.append(gffRGIEntry)
			outCoordsList.append([gffSeqName, bestHitAROName, gffStart, gffStop, gffStrand])

	return (outGFFList, outCoordsList, idList)

def processGFFFile (inGFFFileName):
	try:
		inGFFFile = open(inGFFFileName, 'r')
	except IOError as e:
		exit("IO Error: Could not locate " + e.filename + ".")
	# List of split lines in GFF file
	gffTableList = []
	for line in inGFFFile:
		if line.startswith("##FASTA"):
			break
		elif line.startswith("#"):
			continue
		line = line.strip()
		splitLine = line.split("\t")
		gffTableList.append(splitLine)
	inGFFFile.close()
	return gffTableList

def processRGITable(inRGITableFileName):
	try:
		inRGITableFile = open(inRGITableFileName, 'r')
	except IOError as e:
		exit("IO Error: Could not locate " + e.filename + ".")
	rgiTableList = []
	header = inRGITableFile.readline()
	for line in inRGITableFile:
		line = line.strip('\n')
		splitLine = line.split("\t")
		rgiTableList.append(splitLine)
	inRGITableFile.close()
	# If no results exit quietly
	if len(rgiTableList) == 0:
		print "No AR Results"
		exit(0)
	return rgiTableList


def printListToFile (outList, prefix):
	try:
		print "Writing to file:"+prefix
		outFile = open(prefix, 'w')
	except IOError as e:
		exit("ERROR: Could not locate " + e.filename + ". Please make sure it exists.")

	for line in outList:
		outFile.write("\t".join(line)+"\n")
	outFile.close()

def printJsonToFile (jsonDict, prefix):
	print "Writing to File"+prefix
	try:
		outFile = open(prefix, 'w')
	except IOError as e:
		exit("ERROR: Could not locate " + e.filename + ". Please make sure it exists.")

	json.dump(jsonDict, outFile, indent=4, sort_keys=True)

	outFile.close()

#------------------------------------ Main ------------------------------

def main(args):
	if args.inGFF and args.inRGITable:
		gffTableList = processGFFFile(args.inGFF)
		rgiTableList = processRGITable(args.inRGITable)
		outGFFList = []
		outCoordsList = []
		idList = []
	else:
		exit("Please provide a gff file and rgi table file")

	prefix = args.prefix
	categoryDict = {}

	print "Processing RGI results"

	for rgiHitList in rgiTableList:
		# print rgiHitList[0]

		(orfID, contig, start, stop, orientation, cutOff, passEvalue, bestHitBitscore, bestHitAROName, bestHitARO, bestHitCategories, identities, aros, aroNames, modelType, snps, aroCategories, bitscores, predictedProtein, cardProtein, label, hitID) = rgiHitList

		(outGFFList, outCoordsList, idList) = getGFFInfo(rgiHitList, gffTableList, outGFFList, outCoordsList, idList)

		arData = {'orfID': orfID,
		          'bestHitName': bestHitAROName,
		          'bestHitARO': bestHitARO,
		          'bestHitCategories': bestHitCategories,
		          'snps': snps,
		          'otherHits': aroNames}

		for category in bestHitCategories.split(","):
			category = category.strip()
			if category == "":
				category = "No_category"
			if category in categoryDict:
				categoryDict[category].append(arData)
			else:
				categoryDict[category] = [arData]

	print "Done processing RGI results"

	printListToFile(outGFFList, prefix + ".gff")
	printListToFile(outCoordsList, prefix + "_coords.txt")
	printJsonToFile(categoryDict, prefix + "_table.json")
	print "Done"


def run():
	parser = argparse.ArgumentParser(description='Process Antibiotic Resistance Results for Resistance Gene Identifier (RGI).  To be used with ORF analysis only')
	parser.add_argument('-i', '--inRGITable', dest='inRGITable', help='RGI result table (after converting the json to table)')
	parser.add_argument('-g', '--gff', dest='inGFF', help='Project annotation GFF file, used with ORF analysis only')
	parser.add_argument('-p', '--prefix', dest='prefix', help='Prefix used for output')
	args = parser.parse_args()
	main(args)

if __name__ == '__main__':
	run()