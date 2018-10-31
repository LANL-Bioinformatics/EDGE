#!/usr/bin/env python
import sys
import os
import csv
import re
import ConfigParser
import json
from Bio import SeqIO
import argparse


#---------------------- End RGI filepaths.py ---------------------------

# Code from conda RGI convertJsonToTSV.py Version 3.1.1
# Edited by Logan Voegtly

#output the information particular field from alignment.Title by splicing it by '|'
def findnthbar(bunchstr, start):

	barc = 0
	over = start+1
	temp = ""

	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if barc == start:
			if eachc == '|':
				pass
			else:
				temp += eachc
		if barc == over:
			break

	return temp


#output the information particular field from alignment.Title by splicing it by '#'
def findnthbar2(bunchstr, n):
	arr = bunchstr.split("#")
	if n < len(arr):
		# gene id
		if n == 1 and arr[n]:
			return int(arr[n])
		elif n == 2:
			return int(arr[n])
		elif n == 3:
			if int(arr[n]) == 1:
				# positive
				return "+"
			else:
				# neg
			    return "-"
		else:
			return arr[n]
	else:
		return ""


def findORFfrom (bunchstr):
	barc = 0
	start = 6
	temp = ""
	allout = False

	for eachc in bunchstr:
		if eachc == '|':
			barc += 1
		if allout or barc == start:
			allout = True
			temp += eachc

	return temp[1:]


def convert(input):
	if isinstance(input, dict):
		return dict((convert(key), convert(value)) for key, value in input.iteritems())
	elif isinstance(input, list):
		return [convert(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input


def checkKeyExisted(key, my_dict):
	try:
		nonNone = my_dict[key] is not None
	except KeyError:
		nonNone = False
	return nonNone


# Logan Voegtly modified function
# Function originally printCSV

def processJson(resultfile, orf):
	try:
		with open(resultfile, 'r') as f:
			data = json.load(f)

	except ValueError:
		print>> sys.stderr, "convertJsonToTSV expects a file contains a VALID JSON string."
		exit()

	rgiTable = []
	rgiTable.append(
		["ORF_ID", "CONTIG", "START", "STOP", "ORIENTATION", "CUT_OFF", "PASS_BitScore", "Best_Hit_BitScore", "Best_Hit_evalue",
		 "Best_Hit_ARO_Name", "Best_Hit_ARO", "Best_Hit_Categories", "Best_Identities", "ARO", "ARO_name",
		 "Model_type", "SNP", "ARO_category", "bit_score", "Predicted_Protein", "CARD_Protein_Sequence", "LABEL", "ID"])
	for item in data:
		bestevalue = False
		maxBitScore = False
		startCompare = False
		bestARO = 0
		bestAROCat = 0
		bestAROName = 0
		AROlist = []
		AROnameList = []
		bitScoreList = []
		AROcatList = []
		snpList = []
		cutoffList = []
		typeList = []
		evalueList = []
		identityList = []
		SequenceFromBroadStreet = ""
		predictedProtein = ""
		geneID = ""
		hitID = ""

		if item not in ["_metadata", "data_type"]:
			geneID = item

			for it in data[item]:
				cgList = []
				AMRfamList = []
				if checkKeyExisted("ARO_category", data[item][it]):
					for aroctkey in data[item][it]["ARO_category"]:
						if(data[item][it]["ARO_category"][aroctkey]["category_aro_class_name"] == "AMR Gene Family"):                                                                   
						    AMRfamList.append(str(data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii','replace')))     
						else:
							cgList.append(str(
							data[item][it]["ARO_category"][aroctkey]["category_aro_name"].encode('ascii',
							                                                                     'replace')))                                                                
						                                                                    		

				if data[item][it]["model_type_id"] == 40293:
					temp = data[item][it]["snp"]["original"] + str(data[item][it]["snp"]["position"]) + \
					       data[item][it]["snp"]["change"]
					snpList.append(convert(temp))
				elif data[item][it]["model_type_id"] == 40292:
					snpList.append("n/a")
				
				AROlist.append(convert(data[item][it]["ARO_accession"]))
				AROnameList.append(convert(data[item][it]["ARO_name"]))
				bitScoreList.append(data[item][it]["bit_score"])
				pass_bitscore = str(data[item][it]["pass_bitscore"]).split("|")[0]
				AROcatList.append(cgList)
				typeList.append(convert(data[item][it]["model_type"]))
				cutoffList.append(convert(data[item][it]["type_match"]))
				idenPercent = float(data[item][it]["max_identities"]) / len(data[item][it]["query"])
				'''print>>sys.stderr, data[item][it]["max_identities"]
				print>>sys.stderr, len(data[item][it]["query"])
				print (str(269/289) + "haha")
				print>>sys.stderr, float(data[item][it]["max_identities"] % len(data[item][it]["query"]))'''
				identityList.append(idenPercent)

				if startCompare:
					if maxBitScore < data[item][it]["bit_score"]:
						bestevalue = data[item][it]["evalue"]
						maxBitScore = data[item][it]["bit_score"]
						bestARO = data[item][it]["ARO_accession"]
						bestAROName = data[item][it]["ARO_name"]
						bestAROCat = AMRfamList
						SequenceFromBroadStreet = data[item][it]["sequence_from_broadstreet"]
						if "orf_prot_sequence" in data[item][it]:
							predictedProtein = data[item][it]["orf_prot_sequence"]
						if "hsp_num:" in it:
							hitID = it
				else:
					startCompare = True
					maxBitScore = data[item][it]["bit_score"]
					bestevalue = data[item][it]["evalue"]
					bestAROName = data[item][it]["ARO_name"]
					bestARO = data[item][it]["ARO_accession"]
					bestAROCat = AMRfamList
					SequenceFromBroadStreet = data[item][it]["sequence_from_broadstreet"]
					if "orf_prot_sequence" in data[item][it]:
						predictedProtein = data[item][it]["orf_prot_sequence"]
					if "hsp_num:" in it:
						hitID = it

		clist = set(cutoffList)
		tl = set(typeList)
		arocatset = set(AROnameList)

		if set(snpList) == set(['n/a']):
			snpList = 'n/a'
		else:
			snpList = ', '.join(snpList)

		from itertools import chain
		AROcatList = list(chain.from_iterable(AROcatList))
		AROcatalphaSet = set(AROcatList)
		AROsortedList = sorted(list(AROcatalphaSet))
		
		if typeList:
			if orf == "1":
				# for protein RGI runs where there's no | or seq_start/stop/strand
				if findnthbar(item, 4) == "":
					rgiTable.append([item, "", "", "", "", ', '.join(list(clist)), pass_bitscore, maxBitScore, bestevalue, bestAROName, bestARO, ','.join(bestAROCat),
					                 max(identityList), ', '.join(map(lambda x: "ARO:" + x, AROlist)),
					                 ', '.join(list(arocatset)), ', '.join(list(tl)), snpList,
					                 ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)), predictedProtein,
					                 SequenceFromBroadStreet, geneID, hitID])
				else:
					rgiTable.append([findnthbar(item, 0), findORFfrom(item), int(findnthbar(item, 4)) - 1,
					                 int(findnthbar(item, 5)) - 1, findnthbar(item, 3), ', '.join(list(clist)),
					                 pass_bitscore, maxBitScore, bestevalue, bestAROName, bestARO, ','.join(bestAROCat), max(identityList),
					                 ', '.join(map(lambda x: "ARO:" + x, AROlist)), ', '.join(list(arocatset)),
					                 ', '.join(list(tl)), snpList, ', '.join(AROsortedList),
					                 ', '.join(map(str, bitScoreList)), predictedProtein, SequenceFromBroadStreet,
					                 geneID, hitID])
			else:
				if findnthbar2(item, 1) == "":
					rgiTable.append([item, "", "", "", "", ', '.join(list(clist)), pass_bitscore, maxBitScore, bestevalue, bestAROName, bestARO, ','.join(bestAROCat),
					                 max(identityList), ', '.join(map(lambda x: "ARO:" + x, AROlist)),
					                 ', '.join(list(arocatset)), ', '.join(list(tl)), snpList,
					                 ', '.join(AROsortedList), ', '.join(map(str, bitScoreList)), predictedProtein,
					                 SequenceFromBroadStreet, geneID, hitID])
				else:
					rgiTable.append([findnthbar2(item, 0),
					                 findnthbar2(item, 4).strip(" "),
					                 int(findnthbar2(item, 1)) - 1,
					                 int(findnthbar2(item, 2)) - 1,
					                 findnthbar2(item, 3),
					                 ', '.join(list(clist)), pass_evalue, maxBitScore, pass_bitscore, bestAROName, bestARO, ','.join(bestAROCat), max(identityList),
					                 ', '.join(map(lambda x: "ARO:" + x, AROlist)), ', '.join(list(arocatset)),
					                 ', '.join(list(tl)), snpList, ', '.join(AROsortedList),
					                 ', '.join(map(str, bitScoreList)), predictedProtein, SequenceFromBroadStreet,
					                 geneID, hitID])

	return rgiTable
#------------------------------- End RGI convertJsonToTSV.py --------------------------


def getGFFInfo (rgiHitList, gffTableList, outGFFList, outCoordsList, idList):
	(orfID, contig, start, stop, orientation, cutOff, PASS_BitScore, bestHitBitScore, bestHitEvalue, bestHitAROName,
	 bestHitARO, bestHitCategories, identities, aros, aroNames, modelType, snps, aroCategories, bitscores,
	 predictedProtein, cardProtein, label, hitID) = rgiHitList

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
			newAttributeList.append("rgi_pass_bitscore=" + PASS_BitScore)
			newAttributeList.append("rgi_ARO_categories=" + aroCategories)
			newAttributeList.append("rgi_best_bitscore=" + bestHitBitScore)
			newAttributeList.append("rgi_best_evalue=" + bestHitEvalue)
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


def convertListOfLists(inListOfLists):
	newListOfLists = []
	for inList in inListOfLists:
		listString = '\t'.join(str(v) for v in inList)
		newListOfLists.append(listString.split('\t'))
	return newListOfLists
#------------------------------------ Main ------------------------------

def main(args):
	if args.inGFF and args.inRGIJson:
		gffTableList = processGFFFile(args.inGFF)
		print "Processing Json file"
		rgiTableList = processJson(args.inRGIJson, 0)
		print "Done processing Json file"
		rgiTableHeader = rgiTableList.pop(0)
		rgiTableList = convertListOfLists(rgiTableList)
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

		(orfID, contig, start, stop, orientation, cutOff, passBitScore, bestHitBitScore, bestHitEvalue, bestHitAROName, bestHitARO, bestHitCategories, identities, aros, aroNames, modelType, snps, aroCategories, bitscores, predictedProtein, cardProtein, label, hitID) = rgiHitList

		(outGFFList, outCoordsList, idList) = getGFFInfo(rgiHitList, gffTableList, outGFFList, outCoordsList, idList)

		arData = {'orfID': orfID,
		          'bestHitName': bestHitAROName,
		          'bestHitARO': bestHitARO,
		          'bestHitCategories': bestHitCategories,
		          'bestHitBitScore': bestHitBitScore,
		          'bestHitEvalue': bestHitEvalue,
		          'cutOff': cutOff,
		          'passBitScore': passBitScore,
		          'aroCategory': aroCategories,
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
	rgiTableList.insert(0, rgiTableHeader)
	printListToFile(rgiTableList, prefix + ".txt")
	printJsonToFile(categoryDict, prefix + "_table.json")
	print "Done"


def run():
	parser = argparse.ArgumentParser(description='Process Antibiotic Resistance Results for Resistance Gene Identifier (RGI).  To be used with ORF analysis only')
	parser.add_argument('-i', '--inRGIJson', dest='inRGIJson', help='RGI result json file')
	parser.add_argument('-g', '--gff', dest='inGFF', help='Project annotation GFF file, used with ORF analysis only')
	parser.add_argument('-p', '--prefix', dest='prefix', help='Prefix used for output')
	args = parser.parse_args()
	main(args)


if __name__ == '__main__':
	run()