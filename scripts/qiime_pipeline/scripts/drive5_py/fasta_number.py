#!/usr/bin/python

import sys
import die

Prefix = ""
if len(sys.argv) > 2:
	Prefix = sys.argv[2]

NeedSize = 0
if len(sys.argv) > 3:
	if sys.argv[3] == "-needsize":
		NeedSize = 1
	elif sys.argv[3] == "-nosize":
		NeedSize = 0
	else:
		die.Die("Must specify -needsize or -nosize")

def GetSize(Label):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("size="):
			return int(Field[5:])
	print >> sys.stderr
	print >> sys.stderr, "Size not found in label: " + Label
	sys.exit(1)

File = open(sys.argv[1])
N = 0
while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	Line = Line[:-1]
	if len(Line) == 0:
		continue
	if Line[0] == '>':
		N += 1
		if NeedSize:
			Label = Line[1:].strip()
			Size = GetSize(Label)
			print ">%s%u;size=%u;" % (Prefix, N, Size)
		else:
			print ">%s%u" % (Prefix, N)
	else:
		print Line
