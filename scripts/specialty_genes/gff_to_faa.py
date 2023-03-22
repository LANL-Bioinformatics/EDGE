#!/usr/bin/env python
import sys
from Bio import SeqIO

faa = sys.argv[1];
gff = sys.argv[2];


record_dict = SeqIO.index(faa, "fasta")

inGFFFile = open(gff, 'r')
for line in inGFFFile:
    if line.startswith("##FASTA"):
        break
    elif line.startswith("#"):
        continue

    splitLine = line.split("\t")
    attr = splitLine[8].split(";")
    
    for gffattr in attr:
        if "locus_tag=" in gffattr:
            (tmp,pid) = gffattr.split("=")
            fastaname = pid + "\t"
        if "product=" in gffattr:
            (tmp,product) = gffattr.split("=")
            fastaname += product + ";"
        if "rgi_ARO_categories=" in gffattr:
            (tmp,aro) = gffattr.split("=")
            fastaname += aro + ";"
        if "vf=" in gffattr:
            fastaname += gffattr + ";"
        if "vfclass=" in gffattr:
            fastaname += gffattr + ";"
    print(">%s\n%s" %(fastaname,record_dict[pid].seq))
