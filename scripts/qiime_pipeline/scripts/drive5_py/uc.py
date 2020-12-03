import re
import sys
import progress

# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=Label
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NotMatched
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.

MaxError = -1

Type = '?'
ClusterNr = -1
Size = -1
PctId = -1.0
LocalScore = -1.0
Evalue = -1.0
Strand = '.'
QueryStart = -1
SeedStart = -1
Alignment = ""
QueryLabel = ""
TargetLabel = ""
FileName = "?"
Line = ""

def Die(s):
	print >> sys.stderr, "*** ERROR ***", s, sys.argv
	sys.exit(1)

def ProgressFile(File, FileSize):
#	if not sys.stderr.isatty():
#	return
	Pos = File.tell()
	Pct = (100.0*Pos)/FileSize
	Str = "%s %5.1f%%\r" % (FileName, Pct)
	sys.stderr.write(Str)

def Progress(i, N):
#	if not sys.stderr.isatty():
	return
	Pct = (100.0*i)/N
	Str = "%5.1f%%\r" % Pct
	sys.stderr.write(Str)

def PrintLine():
	print Line

def ParseRec(Line):
	global Type
	global ClusterNr
	global Size
	global PctId
	global Strand
	global QueryStart
	global SeedStart
	global Alignment
	global QueryLabel
	global TargetLabel
	global LocalScore
	global Evalue
	
	Fields = Line.split("\t")
	N = len(Fields)
	if N != 9 and N != 10:
		Die("Expected 9 or 10 fields in .uc record, got: " + Line)
	Type = Fields[0]
	
	try:
		ClusterNr = int(Fields[1])
	except:
		ClusterNr = -1
		
	try:	
		Size = int(Fields[2])
	except:
		Size = -1

	Fields2 = Fields[3].split('/')
	LocalScore = -1.0
	Evalue = -1.0
	if len(Fields2) == 3:
		try:
			PctId = float(Fields2[0])
			LocalScore = float(Fields2[1])
			Evalue = float(Fields2[2])
		except:
			PctId = -1.0
	else:
		try:
			PctId = float(Fields[3])
		except:
			PctId = -1.0

	Strand = Fields[4]
	
	try:
		QueryStart = int(Fields[5])
	except:
		QueryStart = -1

	try:
		SeedStart = int(Fields[6])
	except:
		SeedStart = -1

	Alignment = Fields[7]
	QueryLabel = Fields[8]
	
	if len(Fields) > 9:
		TargetLabel = Fields[9]

def GetRec(File, OnRecord):
	global Line
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return 0
		if Line[0] == '#':
			continue
		Line = Line.strip()
		if len(Line) == 0:
			return 1
		ParseRec(Line)
		Ok = OnRecord()
		if Ok != None and Ok == 0:
			return 0
		return 1

def ReadRecs(argFileName, OnRecord, ShowProgress = True):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def ReadRecsOnRec(argFileName, OnRecord, ShowProgress = True):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def GetRecs(argFileName, OnRecord, ShowProgress = True):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def ReadFile(argFileName, OnRecord, ShowProgress = True):
	global FileName
	FileName = argFileName
	File = open(FileName)

	if ShowProgress:
		progress.InitFile(File, FileName)
	while GetRec(File, OnRecord):
		if ShowProgress:
			progress.File()
	if ShowProgress:
		progress.FileDone()
