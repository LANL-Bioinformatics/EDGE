import sys
import time

File__ = None
FileSize__ = None
FileName__ = None
Secs__ = None

def InitFile(File, FileName = ""):
	global Secs__, File__, FileSize__, FileName__

	File__ = File
	FileName__ = FileName
	Secs__ = None

	Pos = File.tell()
	File.seek(0, 2)
	FileSize__ = File.tell()
	File.seek(Pos)

def FileDone(Msg = ""):
	global Secs__, File__, FileSize__, FileName__
	Str = "%s 100.0%% %s  \n" % (FileName__, Msg)
	sys.stderr.write(Str)

def File(Msg = ""):
	global Secs__, File__, FileSize__, FileName__

	Secs = time.clock()
	if Secs__ != None and Secs - Secs__ < 1:
		return

	Secs__ = Secs
	Pos = File__.tell()
	Pct = (100.0*Pos)/FileSize__
	Str = "%s %5.1f%% %s  \r" % (FileName__, Pct, Msg)
	sys.stderr.write(Str)

def FileStep(Msg = ""):
	File(Msg)

def Step(Msg, i, N):
	global Secs__, File__, FileSize__, FileName__

	Secs = time.clock()
	if Secs__ != None and Secs - Secs__ < 1:
		return

	Secs__ = Secs
	Pct = (100.0*i)/N
	if i == N-1:
		sys.stderr.write("%5.1f%% %s   \r" % (Pct, Msg))
	else:
		sys.stderr.write("%5.1f%% %s   \n" % (Pct, Msg))
