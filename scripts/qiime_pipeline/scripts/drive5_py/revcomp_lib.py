#!/usr/bin/python

import sys
import string
import fasta

Map = {}
p = string.printable
for i in range(0, len(p)):
	c = p[i]
	Map[c] = c

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

Syms = "TGCAKYWSRMBDHVXN"

Map['A'] = 'T'
Map['C'] = 'G'
Map['G'] = 'C'
Map['T'] = 'A'
Map['U'] = 'A'

Map['M'] = 'K'
Map['R'] = 'Y'
Map['W'] = 'W'
Map['S'] = 'S'
Map['Y'] = 'R'
Map['K'] = 'M'
Map['V'] = 'B'
Map['H'] = 'D'
Map['D'] = 'H'
Map['B'] = 'V'
Map['X'] = 'X'
Map['N'] = 'N'

N = len(Syms)
for i in range(0, N):
	C = Syms[i]
	c = C.lower()
	R = Map[C]
	assert Map[R] == C
	r = R.lower()
	Map[r] = c

def RevComp(s):
	global Map
	t = ""
	n = len(s)
	for i in range(0, n):
		c = s[n-i-1]
		t += Map[c]
	return t
