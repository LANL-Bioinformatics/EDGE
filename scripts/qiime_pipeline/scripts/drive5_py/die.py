import sys
import traceback

def Die(Msg):
	print >> sys.stderr
	print >> sys.stderr

	traceback.print_stack()
	s = ""
	for i in range(0, len(sys.argv)):
		if i > 0:
			s += " "
		s += sys.argv[i]
	print >> sys.stderr, s
	print >> sys.stderr, "**ERROR**", Msg
	print >> sys.stderr
	print >> sys.stderr
	sys.exit(1)
	print "NOTHERE!!"

def Warning(Msg):
	print >> sys.stderr
	print >> sys.stderr, sys.argv
	print >> sys.stderr, "**WARNING**", Msg
