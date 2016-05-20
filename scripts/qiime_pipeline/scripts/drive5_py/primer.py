import die

# 	Code	Means				Comp	CompCode
# 	{ 'M', "AC",   'K' },		// GT		K
# 	{ 'R', "AG",   'Y' },		// CT		Y
# 	{ 'W', "AT",   'W' },		// AT		W
# 	{ 'S', "CG",   'S' },		// CG		S
# 	{ 'Y', "CT",   'R' },		// AG		R
# 	{ 'K', "GT",   'M' },		// AC		M
# 	{ 'V', "ACG",  'B' },		// CGT		B
# 	{ 'H', "ACT",  'D' },		// AGT		D
# 	{ 'D', "AGT",  'H' },		// ACT		H
# 	{ 'B', "CGT",  'V' },		// ACG		V
# 	{ 'X', "GATC", 'X' },		// ACGT		X
# 	{ 'N', "GATC", 'N' },		// ACGT		N

LetterToSet = {}
LetterToSet['A'] = "A"
LetterToSet['C'] = "C"
LetterToSet['G'] = "G"
LetterToSet['T'] = "T"
LetterToSet['M'] = "AC"
LetterToSet['R'] = "AG"
LetterToSet['W'] = "AT"
LetterToSet['S'] = "CG"
LetterToSet['Y'] = "CT"
LetterToSet['K'] = "GT"
LetterToSet['V'] = "ACG"
LetterToSet['H'] = "ACT"
LetterToSet['D'] = "AGT"
LetterToSet['B'] = "CGT"
LetterToSet['X'] = "GATC"
LetterToSet['N'] = "GATC"

def MergeChars(a, b):
	global LetterToSet
	if a == b:
		return a
	try:
		sa = LetterToSet[a.upper()]
	except:
		die.Die("Bad letter in primer '%c'" % a)

	try:
		sb = LetterToSet[b.upper()]
	except:
		die.Die("Bad letter in primer '%c'" % b)

	Set = sa
	for x in sb:
		if x not in Set:
			Set += x

	for x in LetterToSet.keys():
		Set2 = LetterToSet[x]
		if len(Set2) != len(Set):
			continue
		Found = True
		for y in Set:
			if y not in Set2:
				Found = False
		if Found:
			return x

	die.Die("MergeChars(%c,%c) failed" % (a, b))

def Expand(Primer):
	Seqs = [ Primer ]
	L = len(Primer)
	for i in range(0, L):
		c = Primer[i].upper()
		Set = LetterToSet[c]
		if len(Set) > 1:
			NewSeqs = [ ]
			for Seq in Seqs:
				for d in Set:
#					print "c=", c, "d=", d, "i=", i, "Seq[:i]", Seq[:i], "Seq[i+1:]", Seq[i+1:]
					Seq = Seq[:i] + d + Seq[i+1:]
					NewSeqs.append(Seq)
			Seqs = NewSeqs
	return Seqs
				
def GetDengen(Primer):
	d = 1
	for c in Primer:
		d *= GetDegenChar(c)
	return d

def MatchLetter(a, b):
	global LetterToSet
	try:
		sa = LetterToSet[a.upper()]
	except:
		return False

	try:
		sb = LetterToSet[b.upper()]
	except:
		return False

	for ca in sa:
		if ca in sb:
			return True

	return False

def MatchPrefix(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	n = PrimerLength
	if L < n:
		n = L
	Diffs = 0
	for i in range(0, n):
		if not MatchLetter(Seq[i], Primer[i]):
			Diffs += 1
	return Diffs

def GetDiffs(Primer1, Primer2):
	assert len(Primer1) == len(Primer2)
	return MatchPrefix(Primer1, Primer2)

def MatchPrefix2(Seq, Primer, MaxDiffs):
	L = len(Seq)
	PrimerLength = len(Primer)
	n = PrimerLength
	if L < n:
		n = L
	Diffs = 0
	for i in range(0, n):
		if not MatchLetter(Seq[i], Primer[i]):
			Diffs += 1
			if Diffs > MaxDiffs:
				return Diffs
	return Diffs

def MatchPos(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	for Pos in range(0, L-PrimerLength):
		if MatchPrefix(Seq[Pos:], Primer) == 0:
			return Pos
	return -1

def Match(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	for Pos in range(0, L-PrimerLength):
		if MatchPrefix(Seq[Pos:], Primer) == 0:
			return Pos
	return -1

def BestMatch(Seq, Primer):
	L = len(Seq)
	PrimerLength = len(Primer)
	BestDiffs = PrimerLength
	BestPos = -1
	for Pos in range(0, L-PrimerLength+1):
		d = MatchPrefix(Seq[Pos:], Primer)
		if d < BestDiffs:
			BestDiffs = d
			BestPos = Pos
	return BestPos, BestDiffs


def BestMatch2(Seq, Primer, MaxDiffs):
	L = len(Seq)
	PrimerLength = len(Primer)
	BestDiffs = PrimerLength
	BestPos = -1
	for Pos in range(0, L-PrimerLength+1):
		d = MatchPrefix2(Seq[Pos:], Primer, MaxDiffs)
		if d < BestDiffs and d <= MaxDiffs:
			BestDiffs = d
			BestPos = Pos
	return BestPos, BestDiffs

def GetDegen(Primer):
	Degen = 1
	print "Primer='" + Primer + '"'
	for Letter in Primer:
		Set = LetterToSet[Letter]
		n = len(Set)
		Degen *= n
	return Degen

# In FASTA or txt
def FromFile(FileName):
	Seqs = []
	File = open(FileName)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return Seqs
		if Line.startswith(">"):
			continue
		Seq = Line[:-1]
		Seqs.append(Seq.upper())
