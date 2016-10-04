#!/usr/bin/env python
import sys
import argparse
import re
import mysql.connector as mysql
from mysql.connector import errorcode
import ConfigParser
import json
import subprocess
from Bio import SeqIO

#-------------------------- Functions -----------------------------------

def connectMysql(mysqlConfigName):
	mysqlConfig = ConfigParser.ConfigParser()
	try:
		mysqlConfig.read(mysqlConfigName)
	except IOError as e:
		exit("IO Error: Could not open " + e.filename + ".")
	# Read in mysql credentials
	mysqlUsername = mysqlConfig.get('LOGIN', 'username')
	mysqlPassword = mysqlConfig.get('LOGIN', 'password')
	mysqlHost = mysqlConfig.get('LOGIN', "host")
	mysqlPort = mysqlConfig.get('LOGIN', 'port')
	mysqlDatabase = mysqlConfig.get('LOGIN', 'database')

	# Create connection into database
	try:
		mysqlCnx = mysql.connect(user=mysqlUsername, password=mysqlPassword, host=mysqlHost, port=mysqlPort,
		                         database=mysqlDatabase)
	except mysql.Error as err:
		if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
			print("Something is wrong with your user name or password")
			exit()
		elif err.errno == errorcode.ER_BAD_DB_ERROR:
			print("Database does not exist")
			exit()
		else:
			print(err)
			exit()
	return mysqlCnx


def getMysqlResults (mysqlCursor, sqlStatement, inList):
	try:
		mysqlCursor.execute(sqlStatement, inList)
		mysqlResult = mysqlCursor.fetchall()
	except mysql.Error as err:
		exit(err)
	return mysqlResult


def processResultFile (inSBResultsFileName):
	try:
		inSBResultFile = open(inSBResultsFileName, 'r')
	except IOError as e:
		exit("IO Error: Could not open " + e.filename + ".")

	# Put results GI and score into dict
	noSBResultCount = 0
	totalVFsCount = 0
	sbResultCount = 0
	sbResultDict = {}
	# Remove header line
	inSBResultFile.readline()
	for line in inSBResultFile:
		line = line.strip()
		splitLine = line.split("\t")
		score = float(splitLine[1])
		gi = splitLine[0].split("|")[3]
		totalVFsCount += 1
		if score > 0:
			sbResultCount += 1
			sbResultDict[gi] = str(score)
		else:
			noSBResultCount += 1
	if len(sbResultDict) == 0:
		print "No Virulence Results"
		sys.exit(0)

	inSBResultFile.close()
	return sbResultDict


def processSBHitsFile (inSBHitsFileName):
	try:
		inSBHitsFile = open(inSBHitsFileName, 'r')
	except IOError as e:
		exit("IO Error: Could not locate " + e.filename + ".")

	# SBhits File from ShortBRED output
	sbHitsList = []
	for line in inSBHitsFile:
		line = line.strip()
		splitLine = line.split("\t")
		sbHitsList.append(splitLine)
	inSBHitsFile.close()

	# Dictionary of marker GI and set of genes
	sbGeneSetDict = {}
	# Loop through the sbHits
	for sbHit in sbHitsList:
		# Gets GI for the SB Hit
		sbMarkerGI = sbHit[0].split("|")[3]
		# gets the called gene from edge like (OBrien_E_coli_O157s_933Sr_S21_00428)
		sbGene = sbHit[1]
		try:
			# Add to the gene to the set within the dictionary with gi as the key
			# There may be cases where the same gi will map to multiple genes
			sbGeneSetDict[sbMarkerGI].add(sbGene)
		except KeyError:
			# If the gi key does not exist
			sbGeneSetDict[sbMarkerGI] = set()
			sbGeneSetDict[sbMarkerGI].add(sbGene)
	return sbGeneSetDict


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


def generateKnownClassDict(mysqlResult):
	# Generate dict of known (not Unknown) classes and VF's
	knownClassDict = {}
	for resultList in mysqlResult:
		currentClass = resultList[0]
		currentVF = resultList[1]
		currentVFNumber = resultList[4]
		if currentClass != "Unknown":
			try:
				if currentClass in knownClassDict:
					if currentVF in knownClassDict[currentClass]:
						knownClassDict[currentClass][currentVF].add(currentVFNumber)
					else:
						knownClassDict[currentClass][currentVF] = set()
						knownClassDict[currentClass][currentVF].add(currentVFNumber)
				else:
					knownClassDict[currentClass] = {}
					knownClassDict[currentClass][currentVF] = set()
					knownClassDict[currentClass][currentVF].add(currentVFNumber)
			except KeyError:
				pass
	return knownClassDict


def generateGeneList (sbGeneSetDict, gffTableList):
	# Dict with gi as key and set of gff entries
	gffGeneListDict = {}
	# Loop through the dictionary with sets of genes as the value and GIs as key
	for sbMarkerGI, sbGeneSet in sbGeneSetDict.iteritems():
		# Loop through the gene set for the current gi
		for sbGene in sbGeneSet:
			# Loop through the GFF table
			for gffEntry in gffTableList:
				# See if the gene is in the GFF talbe
				if sbGene in gffEntry[8]:
					# Add to gffGeneListDict
					try:
						gffGeneListDict[sbMarkerGI].append(gffEntry)
					except KeyError:
						gffGeneListDict[sbMarkerGI] = [gffEntry]
	return gffGeneListDict


def generateKronaPlot(outResultList, prefix):
	outKronaListFile = open(prefix + "_krona_list.txt", "w")
	listCount = 0
	for result in outResultList:
		# print result, len(result)
		(vfclass, vf, vfGene, generaRepresented, vfNumber, resultGI, vfdbInfo, newGeneAccession, genomeSource, taxid,
		 ncbiComment, resultScore) = result
		kronaResult = [vfNumber, taxid, resultScore]
		if listCount == 0:
			outKronaListFile.write("#" + "\t".join(kronaResult) + "\n")
		else:
			outKronaListFile.write("\t".join(kronaResult) + "\n")
		listCount += 1

	outKronaListFile.close()

	doKrona = True
	if doKrona:
		print "Generating Krona Plot"
		runKronaString = "ktImportTaxonomy -o " + prefix + ".krona.html " + outKronaListFile.name
		subprocess.call(runKronaString.split(' '))


def getGFFInfo (resultList, gffGeneListDict, outGFFList, outCoordsList, idList):
	(vfclass, vf, vfGene, generaRepresented, vfNumber, resultGI, vfdbInfo, newGeneAccession, genomeSource, taxid,
	 ncbiComment, resultScore) = resultList
	for gi, gffEntryList in gffGeneListDict.iteritems():
		if gi == resultGI:
			for gffEntry in gffEntryList:
				(gffSeqName, gffSource, gffFeature, gffStart, gffStop, gffScore, gffStrand, gffFrame, gffAttributes) = gffEntry
				gffAttributeList = gffAttributes.split(";")
				newAttributeList = []
				for gffAttribute in gffAttributeList:
					if "ID=" in gffAttribute:
						newID = "ID=" + vfGene
						if newID in idList:
							# newID += ":"+str(idList.count(newID))
							newID += ":" + vfNumber
						newAttributeList.append(newID)
						idList.append("ID=" + vfGene)
					elif "gene=" in gffAttribute:
						newAttributeList.append(gffAttribute + ":" + vfGene)
					elif "inference=" in gffAttribute:
						newAttributeList.append("inference=ShortBRED")
					elif "locus_tag=" in gffAttribute:
						newAttributeList.append(gffAttribute)
					elif "product=" in gffAttribute:
						newAttributeList.append(gffAttribute)
				gffFeature = "Virulence"
				gffSource = "ShortBRED"
				newAttributeList.append("vfclass=" + vfclass)
				newAttributeList.append("vf=" + vf)
				newAttributeList.append("accession=" + newGeneAccession)
				newAttributeList.append("typically_found_on=" + genomeSource)
				newAttributeList.append("found_in=" + generaRepresented)
				newAttributeList.append("vfinfo=" + vfdbInfo)
				newAttributeList.append("ncbi_comment=" + ncbiComment)
				newAttributeList.append("vf_number=" + vfNumber)
				gffAttributes = ";".join(newAttributeList)
				gffRGIEntry = [gffSeqName, gffSource, gffFeature, gffStart, gffStop, gffScore, gffStrand, gffFrame,
				               gffAttributes]
				outGFFList.append(gffRGIEntry)
				outCoordsList.append([gffSeqName, vfGene, gffStart, gffStop, gffStrand])
	return (outGFFList, outCoordsList, idList)


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

#------------------------------------ Main ----------------------------------------

def main (args):
	processingType = args.processingType
	prefix = args.prefix
	if processingType == 'orf':
		if args.inSBHits and args.inGFF:
			sbGeneSetDict = processSBHitsFile(args.inSBHits)
			gffTableList = processGFFFile(args.inGFF)
			gffGeneListDict = generateGeneList(sbGeneSetDict, gffTableList)
			outGFFList = []
			outCoordsList = []
			idList = []
		else:
			exit("Please provide SBHits file and GFF file when processing ORFs")
	elif processingType == 'reads':
			pass
	else:
		exit("Please provide a processing type reads or orf")

	if args.inSBResults:
		sbResultDict = processResultFile(args.inSBResults)
	else:
		exit("Please provide a ShortBRED results file")

	if args.mysqlConfigName:
		mysqlCnx = connectMysql(args.mysqlConfigName)
	else:
		exit("Please provide a mysql.ini config file")

	jsonDict = {}

	# Create cursor to query database
	mysqlCursor = mysqlCnx.cursor()

	outResultList = []
	outResultList.append(["VFClass", "VF", "VFGene", "VFGeneraRepresented", "VFNumber", "GI", "VFDBInfo",
	                      "NewGeneAccession", "GenomeSource", "TaxID", "NCBIComment", "SBScore"])

	getResultsFromDB = ("SELECT DISTINCT vftable.VFClass, vftable.VF, vftable.VFGene, "
	                    "vftable.VFGeneraRepresented, genes.VFNumber, genes.GI, genes.VFDBInfo, "
	                    "genes.NewGeneAccession, genes.GenomeSource, genes.taxID, genes.NCBIComment "
	                    "FROM vftable, genes, family, family_genes, family_vftable "
	                    "WHERE genes.gi in (" + ",".join(("%s",) * len(sbResultDict)) + ") "
                            "and genes.GeneID = family.ParentGeneID "
                            "and family_vftable.FV_FamilyID=family.FamilyID "
                            "and family_vftable.FV_VFTableID=vftable.VFTableID "
                            "order by vftable.VFGenus, vftable.VFClass, vftable.VFGene; ")
	giTuple = tuple(sbResultDict.keys())

	print "Executing MySQL command to get VF results:"
	print getResultsFromDB
	print "Using the following GI's:"
	print giTuple
	mysqlResult = getMysqlResults(mysqlCursor, getResultsFromDB, giTuple)
	print "Done executing MySQL command"
	print "Preprocessing results"
	knownClassDict = generateKnownClassDict(mysqlResult)
	print "Done preprocessing results"

	print "Processing results"
	for resultList in mysqlResult:
		# Appending score to the result
		resultGI = str(resultList[5])
		resultScore = str(sbResultDict[resultGI])
		resultString = "\t".join(str(v) for v in resultList) + "\t" + resultScore
		resultList = resultString.split("\t")
		(vfclass, vf, vfGene, generaRepresented, vfNumber, resultGI, vfdbInfo, newGeneAccession, genomeSource, taxid,ncbiComment, resultScore) = resultList
		# print resultList
		isNotDuplicate = True
		if vfclass == "Unknown":
			for knownClass, knownVFDict in knownClassDict.iteritems():
				for knownVF, knownVFNumberSet in knownVFDict.iteritems():
					if vfNumber in knownVFNumberSet:
						isNotDuplicate = False
					if knownClass in vf:
						vf = knownVF
						vfclass = knownClass
					else:
						vf = re.sub(" \(\w+\d+\)", "", vf)
						vf = vf.strip()
						vfclass = "Others"
		if isNotDuplicate:
			if processingType == "orf":
				(outGFFList, outCoordsList, idList) = getGFFInfo(resultList, gffGeneListDict, outGFFList, outCoordsList, idList)

			vfdata = {'vfgene': vfGene,
			          'foundin': generaRepresented,
			          'vfnumber': vfNumber,
			          'gi': resultGI,
			          'accession': newGeneAccession,
			          'source': genomeSource,
			          'taxid': taxid,
			          'vfdbinfo': vfdbInfo
			          }

			try:
				if vfclass in jsonDict:
					if vf in jsonDict[vfclass]:
						jsonDict[vfclass][vf].append(vfdata)
					else:
						jsonDict[vfclass][vf] = [vfdata]
				else:
					jsonDict[vfclass] = {vf: [vfdata]}
			except KeyError:
				exit(KeyError.message)

			outResultList.append(resultString.split("\t"))

	print "Done processing results"
	if processingType == 'orf':
		printListToFile (outGFFList, prefix+".gff")
		printListToFile (outCoordsList, prefix+"_coords.txt")
	printListToFile (outResultList, prefix+"_table.txt")
	printJsonToFile(jsonDict, prefix+"_table.json")
	generateKronaPlot(outResultList, prefix)
	print "Done"


def run():
	parser = argparse.ArgumentParser(description='Process Virulence ShortBRED Results')
	parser.add_argument('-t', '--processingType', dest='processingType', help='orf processing or reads processing')
	parser.add_argument('-i', '--sbResults', dest='inSBResults', help='ShortBRED result file')
	parser.add_argument('-m', '--mysqlConfig', dest='mysqlConfigName', help='MySQL configuration file')
	parser.add_argument('-s', '-sbHits', dest='inSBHits', help='ShortBRED SBHits file, used with ORF analysis only')
	parser.add_argument('-g', '--gff', dest='inGFF', help='Project annotation GFF file, used with ORF analysis only')
	parser.add_argument('-p', '--prefix', dest='prefix', help='Prefix used for output')
	args = parser.parse_args()
	main(args)


if __name__ == '__main__':
	run()
