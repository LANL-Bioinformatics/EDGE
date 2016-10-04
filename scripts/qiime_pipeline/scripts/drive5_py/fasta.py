from die import *
import subprocess
import tempfile
import progress

TRUNC_LABELS=0

def isgap(c):
	return c == '-' or c == '.'

def GetSeqCount(FileName):
	Tmp = tempfile.TemporaryFile()
	try:
		TmpFile = Tmp.file
	except:
		TmpFile = Tmp
	s = subprocess.call([ "grep", "-c", "^>", FileName ], stdout=TmpFile)
	TmpFile.seek(0)
	s = TmpFile.read()
	return int(s)

def GetSeqsDict(FileName):
	return ReadSeqsFast(FileName, False)

def ReadSeqsDict(FileName, Progress = False):
	return ReadSeqsFast(FileName, Progress)

def ReadSeqsOnSeq(FileName, OnSeq, Progress = False):
	ReadSeqs3(FileName, OnSeq, Progress)

def ReadSeqsFastFile(File, Progress = False):
	Seqs = {}
	Id = ""
	N = 0
	while 1:
		if N%10000 == 0 and Progress:
			sys.stderr.write("%u seqs\r" % (N))
		Line = File.readline()
		if len(Line) == 0:
			if Progress:
				sys.stderr.write("%u seqs\n" % (N))
			return Seqs
		if len(Line) == 0:
			continue
		Line = Line.strip()
		if Line[0] == ">":
			N += 1
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			Seqs[Id] = ""
		else:
			if Id == "":
				Die("FASTA file does not start with '>'")
			Seqs[Id] = Seqs[Id] + Line

def ReadSeqsFast(FileName, Progress = True):
	File = open(FileName)
	return ReadSeqsFastFile(File, Progress)

def ReadSeqs(FileName, toupper=False, stripgaps=False, Progress=False):
	if not toupper and not stripgaps:
		return ReadSeqsFast(FileName, False)

	Seqs = {}
	Id = ""
	File = open(FileName)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return Seqs
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			if Id in Seqs.keys():
				Die("Duplicate id '%s' in '%s'" % (Id, FileName))
			Seqs[Id] = ""
		else:
			if Id == "":
				Die("FASTA file '%s' does not start with '>'" % FileName)
			if toupper:
				Line = Line.upper()
			if stripgaps:
				Line = Line.replace("-", "")
				Line = Line.replace(".", "")
			Seqs[Id] = Seqs[Id] + Line

def ReadSeqs2(FileName, ShowProgress = True):
	Seqs = []
	Labels = []
	File = open(FileName)
	if ShowProgress:
		progress.InitFile(File, FileName)
	while 1:
		progress.File()
		Line = File.readline()
		if len(Line) == 0:
			if ShowProgress:
				print >> sys.stderr, "\n"
			return Labels, Seqs
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			Labels.append(Id)
			Seqs.append("")
		else:
			i = len(Seqs)-1
			Seqs[i] = Seqs[i] + Line

def ReadSeqs3(FileName, OnSeq, ShowProgress = True):
	File = open(FileName)
	if ShowProgress:
		progress.InitFile(File, FileName)
	Label = ""
	Seq = ""
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			if Seq != "":
				OnSeq(Label, Seq)
			if ShowProgress:
				print >> sys.stderr, "\n"
			return
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if Seq != "":
				if ShowProgress:
					progress.File()
				if TRUNC_LABELS:
					Label = Label.split()[0]
				OnSeq(Label, Seq)
			Label = Line[1:]
			Seq = ""
		else:
			Seq += Line

def WriteSeq(File, Seq, Label = ""):
	if Label != "":
		print >> File, ">" + Label
	BLOCKLENGTH = 80
	SeqLength = len(Seq)
	BlockCount = int((SeqLength + (BLOCKLENGTH-1))/BLOCKLENGTH)
	for BlockIndex in range(0, BlockCount):
		Block = Seq[BlockIndex*BLOCKLENGTH:]
		Block = Block[:BLOCKLENGTH]
		print >> File, Block

def GetSizeFromLabel(Label, Default = -1):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("size="):
			return int(Field[5:])
	if Default == -1:
		Die("Missing size >" + Label)
	return Default

def StripSizeFromLabel(Label):
	Fields = Label.split(";")
	NewLabel = ""
	for Field in Fields:
		if Field.startswith("size="):
			continue
		if NewLabel != "":
			NewLabel += ";"
		NewLabel += Field
	return NewLabel

def GetQualFromLabel(Label):
	n = Label.find("qual=")
	assert n >= 0
	return Label[n+5:-1]

def StripQualFromLabel(Label):
	n = Label.find("qual=")
	assert n >= 0
	return Label[:n]

def GetField(Label, Name, Default):
	Fields = Label.split(';')
	for Field in Fields:
		if Field.startswith(Name + "="):
			n = len(Name) + 1
			return Field[n:]
	if Default == "":
		Die("Field %s= not found in >%s" % (Name, Label))
	return Default

def GetIntFieldFromLabel(Label, Name, Default):
	return int(GetField(Label, Name, Default))

def GetFieldFromLabel(Label, Name, Default):
	return GetField(Label, Name, Default)

def DeleteFieldFromLabel(Label, Name):
	NewLabel = ""
	Fields = Label.split(';')
	for Field in Fields:
		if len(Field) > 0 and not Field.startswith(Name + "="):
			NewLabel += Field + ';'
	return NewLabel

def ReplaceSize(Label, Size):
	Fields = Label.split(";")
	NewLabel = ""
	Done = False
	for Field in Fields:
		if Field.startswith("size="):
			NewLabel += "size=%u;" % Size
			Done = True
		else:
			if Field != "":
				NewLabel += Field + ";"
	if not Done:
		die.Die("size= not found in >" + Label)
	return NewLabel
