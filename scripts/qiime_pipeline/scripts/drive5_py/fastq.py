import sys
import math

# WARNING
# FASTQ formats differ:
# http://en.wikipedia.org/wiki/FASTQ_format
# The constants below must be hard-coded to match the
# format because the format cannot be reliably detected.

ASCII_Offset = 33
Q_Min = 0
Q_Max = 41

def SetSanger():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 64
	Q_Min = 0
	Q_Max = 41

def SetIllumina():
	global ASCII_Offset, Q_Min, Q_Max
	ASCII_Offset = 33
	Q_Min = 0
	Q_Max = 41
	
def GetLine():
	global File

	Line = File.readline()
	assert len(Line) != 0
	return Line.strip()

def IntQualToChar(iq):
	assert iq >= Q_Min and iq <= Q_Max

	return chr(ASCII_Offset + iq)

def CharToIntQual(c):
	global ASCII_Offset

	ic = ord(c)
	iq = ic - ASCII_Offset
	if iq < Q_Min or iq > Q_Max:
		print >> sys.stderr, "c=%c ic=%d iq=%d Q_Min %d Q_Max %d" % (c, ic, iq, Q_Min, Q_Max)
		assert False

	return iq

def IntQualToProb(q):
	return 10**(-q/10.0)

def ProbToIntQual(P):
	return -10*math.log10(P)

def CharToProb(c):
	ic = ord(c)
	iq = ic - ASCII_Offset
	assert iq >= Q_Min and iq <= Q_Max
	return 10**(-iq/10.0)

def GetRec(File):
	Line = File.readline()
	if len(Line) == 0:
		return "", "", ""

	assert Line[0] == '@'
	Label = Line.strip()[1:]

	Seq = File.readline().strip()
	Plus = File.readline().strip()
	assert Plus[0] == '+'

	Qual = File.readline().strip()
	assert len(Seq) == len(Qual)

	return Label, Seq, Qual

def TruncRec(Seq, Qual, TruncQ, TruncN):
	L = len(Seq)
	for i in range(0, L):
		q = Qual[i]
		iq = CharToIntQual(q)
		if iq <= TruncQ:
			return Seq[:i], Qual[:i]
		if TruncN and Seq[i] == 'N':
			return Seq[:i], Qual[:i]
	return Seq, Qual

def TruncRecRev(Seq, Qual, TruncQ, TruncN):
	L = len(Seq)
	for k in range(0, L):
		i = L - k - 1
		q = Qual[i]
		iq = CharToIntQual(q)
		if iq <= TruncQ:
			return Seq[i:], Qual[i:]
		if TruncN and Seq[i] == 'N':
			return Seq[i:], Qual[i:]
	return Seq, Qual

def GetRecTruncQual(File, MinQ):
	Label, Seq, Qual = GetRec(File)
	L = len(Qual)
	for i in range(0, L):
		q = Qual[i]
		iq = CharToIntQual(q)
		if iq < MinQ:
			return Label, Seq[:i], Qual[:i]
	return Label, Seq, Qual

def ReadRecsOnRec(FileName, OnRec):
	ReadSeqs(FileName, OnRec)

def ReadSeqsOnRec(FileName, OnRec):
	ReadSeqs(FileName, OnRec)

def ReadSeqs(FileName, OnRec):
	global File
	File = open(FileName)
	Label = ""
	Seq = ""
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return

		assert Line[0] == '@'
		Label = Line.strip()[1:]

		Seq = GetLine()
		Plus = GetLine()
		assert Plus[0] == '+'

		Qual = GetLine()
		assert len(Seq) == len(Qual)

		OnRec(Label, Seq, Qual)

def ReadRecs(FileName, OnRec):
	ReadSeqs(FileName, OnRec)

def WriteRec(File, Label, Seq, Qual):
	assert len(Seq) == len(Qual)

	print >> File, "@" + Label
	print >> File, Seq
	print >> File, "+"
	print >> File, Qual

def PrintTable():
	for j in range(1, 11):
		s = ""
		for iq in [ j, 10 + j, 20 + j, 30 + j ]:
			if iq > j:
				s += "  |  "
			c = IntQualToChar(iq)
			s += "%2d  %c  %3d  %7.5f" % (iq, c, ord(c), 10**(-iq/10.0))
		print s

def GetAvgQ(Qual):
	SumQ = 0
	L = len(Qual)
	if L == 0:
		return 0.0

	for q in Qual:
		iq = CharToIntQual(q)
		SumQ += iq

	return float(SumQ)/L

def GetAvgP(Qual):
	SumP = 0
	L = len(Qual)
	if L == 0:
		return 0.0

	for q in Qual:
		P = CharToProb(q)
		SumP += P

	return float(SumP)/L

def GetEE(Qual):
	SumP = 0
	L = len(Qual)
	if L == 0:
		return 0.0

	for q in Qual:
		P = CharToProb(q)
		SumP += P

	return SumP
