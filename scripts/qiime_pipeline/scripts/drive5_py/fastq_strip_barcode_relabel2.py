#!/usr/bin/python
# Like fastq_strip_barcode_relabel except:
#	relabels with barcode label instead of barcode sequence.
# Expect seq = <barcode><primer><gene>
# Allow 2 mismatches with primer
# Allow 0 mismatches with barcode
# Strips primer & barcode, adds barcode label to seq label.

import sys
import fasta
import fastq
import primer

MAX_PRIMER_MISMATCHES = 2

FileName = sys.argv[1]
Primer = sys.argv[2]
BarcodeFileName = sys.argv[3]
LabelPrefix = sys.argv[4]

SeqCount = 0
OutCount = 0
BarcodeMismatchCount = 0
PrimerMismatchCount = 0

PL = len(Primer)

Barcodes = fasta.ReadSeqsDict(BarcodeFileName)

def MatchesPrimer(Seq, Primer):
	return primer.MatchPrefix(Seq, Primer)

def FindBarcode(Seq):
	global Barcodes
	for BarcodeLabel in Barcodes.keys():
		Barcode = Barcodes[BarcodeLabel]
#		print ""
#		print " BC", Barcode
#		print "Seq", Seq
		if Seq.startswith(Barcode):
#			print "+++++++++++++++++++++"
			return Barcode, BarcodeLabel
	return "", ""
	
def OnRec(Label, Seq, Qual):
	global PL, LabelPrefix, Barcode, SeqCount, OutCount, BarcodeMismatchCount, PrimerMismatchCount

	SeqCount += 1
	Barcode, BarcodeLabel = FindBarcode(Seq)
	if Barcode == "":
		BarcodeMismatchCount += 1
		return

	BarcodeLength = len(Barcode)
	Seq = Seq[BarcodeLength:]
	Qual = Qual[BarcodeLength:]

	Diffs = MatchesPrimer(Seq, Primer)
	if Diffs > MAX_PRIMER_MISMATCHES:
		PrimerMismatchCount += 1
		return

	OutCount += 1
	Label = LabelPrefix + str(OutCount) + ";barcodelabel=" + BarcodeLabel + ";"
	fastq.WriteRec(sys.stdout, Label, Seq[PL:], Qual[PL:])

fastq.ReadRecs(FileName, OnRec)

print >> sys.stderr, "%10u seqs" % SeqCount
print >> sys.stderr, "%10u matched" % OutCount
print >> sys.stderr, "%10u barcode mismatches" % BarcodeMismatchCount
print >> sys.stderr, "%10u primer mismatches" % PrimerMismatchCount
