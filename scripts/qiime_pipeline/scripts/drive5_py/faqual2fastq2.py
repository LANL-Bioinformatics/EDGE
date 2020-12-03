import sys
import fasta
import fastq
import die

VERBOSE=0

FastaFileName = sys.argv[1]
QualFileName = sys.argv[2]

ff = open(FastaFileName)
fq = open(QualFileName)


FastaData = ff.read()
QualData = fq.read()

print >> sys.stderr, "%u bytes %s" % (len(FastaData), FastaFileName)
print >> sys.stderr, "%u bytes %s" % (len(QualData), QualFileName)

assert FastaData[0] == '>'
assert QualData[0] == '>'

def GetNext(Data, k, IsQual):
	Label = ""
	Seq = ""
	Bytes = len(Data)
	InLabel = True
	while k < Bytes:
		c = Data[k]
		k += 1
		if c == '>':
			break

		if c == '\r' or c == '\n':
			InLabel = False
			if c == '\n' and IsQual and Seq != "":
				Seq += " "
			continue

		if InLabel:
			Label += c
		else:
			Seq += c

	return Label.strip(), Seq.strip(), k

kf = 1
kq = 1
while 1:
	Labelf, Seqf, kf = GetNext(FastaData, kf, False)
	Labelq, Seqq, kq = GetNext(QualData, kq, True)
	if VERBOSE:
		print ""
		print "------------------------------------------"
		print "kf=%u kq=%u" % (kf, kq)
		print "Fasta '%s'" % FastaData[kf-8:kf+8]
		print "Qual '%s'" % QualData[kq-8:kq+8]
		print "Labelf='%s'" % Labelf
		print "Labelq='%s'" % Labelq
		print "Seqf='%s'" % Seqf
		print "Seqq='%s'" % Seqq

	if Labelf == "":
		assert Labelq == ""
		break

	if Labelf != Labelq:
		print >> sys.stderr
		print >> sys.stderr, "LABEL MISMATCH"
		print >> sys.stderr, "Labelf:", Labelf
		print >> sys.stderr, "Labelq:", Labelq
		sys.exit(1)

	L = len(Seqf)
	Seqq = Seqq.replace("\r", " ")
	Seqq = Seqq.replace("\n", " ")
	Seqq = Seqq.strip()
	Quals = Seqq.split()
	LQ = len(Quals)
	if LQ != L:
		print >> sys.stderr, "Fasta='%s'" %  FastaData
		print >> sys.stderr, "Qual='%s'" % QualData
		die.Die("L=%u LQ=%u Labelf=%s" % (L, LQ, Labelf))

	q = ""
	for Qual in Quals:
		iq = int(Qual)
		cq = fastq.IntQualToChar(iq)
		q += cq

	assert len(q) == L
	fastq.WriteRec(sys.stdout, Labelf, Seqf, q)
