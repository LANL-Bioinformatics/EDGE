#!/usr/bin/env python3
from __future__ import print_function
import pysam
import sys
import csv
from copy import copy
from collections import defaultdict
import re

def read_bed_file(fn):
    bedfile = []
    with open(fn) as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab')
        for row in reader:
        # ref start end primername
            bedrow = {}
            bedrow['Reference'] = row[0]
            bedrow['Primer_ID'] = row[3]
            bedrow['direction'] = row[5]
            if bedrow['direction'] == '+':
                bedrow['end'] = int(row[2])
                bedrow['start'] = int(row[1])
            else:
                bedrow['end'] = int(row[1])
                bedrow['start'] = int(row[2])
            bedfile.append(bedrow)
    return bedfile


def check_still_matching_bases(s):
    for flag, length in s.cigartuples:
        if flag == 0:
            return True
    return False

def trim(cigar, s, start_pos, end):
    
    if not end:
        pos = s.pos
    else:
        pos = s.reference_end
    eaten = 0
    while 1:
        ## chomp stuff off until we reach pos
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)

            if args.verbose:
                sys.stderr.write("Chomped a %s, %s" % (flag, length))
        except IndexError:
            sys.stderr.write("Ran out of cigar during soft masking - completely masked read will be ignored")
            break

        if flag == 0:
            ## match
            #to_trim -= length
            eaten += length
            if not end:
                 pos += length
            else:
                 pos -= length
        if flag == 1:
            ## insertion to the ref
            #to_trim -= length
            eaten += length
        if flag == 2:
            ## deletion to the ref
            #eaten += length
            if not end:
                pos += length
            else:
                pos -= length
            pass
        if flag == 4:
            eaten += length
        if not end and pos >= start_pos and flag == 0:
            break
        if end and pos <= start_pos and flag == 0:
            break
        #print >>sys.stderr, "pos:%s %s" % (pos, start_pos)

    extra = abs(pos - start_pos)
    if args.verbose:
        sys.stderr.write("extra %s" % (extra))
    if extra:
        if flag == 0:
            if args.verbose:
                sys.stderr.write("Inserted a %s, %s" % (0, extra))

            if end:
                cigar.append((0, extra))
            else:
                cigar.insert(0, (0, extra))
            eaten -= extra
    if not end:
        s.pos = pos - extra

    if args.verbose:
        sys.stderr.write("New pos: %s" % (s.pos))

    if end:
        cigar.append((4, eaten))
    else:
        cigar.insert(0, (4, eaten))
    oldcigarstring = s.cigarstring
    s.cigartuples = cigar

    #print >>sys.stderr,  s.query_name, oldcigarstring[0:50], s.cigarstring[0:50]

def find_primer(bed, pos, direction, refID):
   # {'Amplicon_size': '1874', 'end': 7651, '#Region': 'region_4', 'start': 7633, 'Coords': '7633', "Sequence_(5-3')": 'GCTGGCCCGAAATATGGT', 'Primer_ID': '16_R'}
   from operator import itemgetter

   ref_alt_ID = re.sub(r'[^\w]','_',refID)
   search_term = re.compile (r'%s|%s' % (refID,ref_alt_ID) )
   
   closest = min([(abs(p['start'] - pos), p['start'] - pos, p) for p in bed if p['direction'] == direction and search_term.search(p['Reference'])], key=itemgetter(0))
   return closest


def go(args):
    if args.report:
        reportfh = open(args.report, "w")

    bed = read_bed_file(args.bedfile)

    counter = defaultdict(int)

    infile = pysam.AlignmentFile("-")
    outfile = pysam.AlignmentFile("-", "wh", template=infile)
    for s in infile:
        cigar = copy(s.cigartuples)
        refname = s.reference_name
        #qname = s.query_name
       
        ## logic - if alignment start site is _before_ but within X bases of
        ## a primer site, trim it off

        if s.is_unmapped:
            sys.stderr.write("%s skipped as unmapped" % (s.query_name))
            continue

        if s.is_supplementary:
            sys.stderr.write("%s skipped as supplementary" % (s.query_name))
            continue

        p1 = find_primer(bed, s.reference_start, '+', refname)
        p2 = find_primer(bed, s.reference_end, '-',refname)

        report = "%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (s.query_name, s.reference_start, s.reference_end, p1[2]['Primer_ID'], p2[2]['Primer_ID'], p1[2]['Primer_ID'], abs(p1[1]), p2[2]['Primer_ID'], abs(p2[1]), s.is_secondary, s.is_supplementary, p1[2]['start'], p2[2]['end'])
#SRR11085797.8084709	14	165	nCoV-2019_1_LEFT_nCoV-2019_1_RIGHT	nCoV-2019_1_LEFT	16	nCoV-2019_1_RIGHT	245	False	False	30	385
        if args.report:
            print(report, file=reportfh)

        if args.verbose:
            sys.stderr.write(report)

        ## if the alignment starts before the end of the primer, trim to that position

        try:
            ## softmask the alignment if left primer start/end inside alignment
           
            if args.start:
                primer_position = p1[2]['start']
            else:
                primer_position = p1[2]['end']
            
            if args.strand:
                if not s.is_reverse and s.reference_start >=  p1[2]['start'] and s.reference_start <=  p1[2]['end']:
                    trim(cigar, s, primer_position, 0)
            else:
                if s.reference_start < primer_position:
                    trim(cigar, s, primer_position, 0)
                else:
                    if args.verbose:
                        sys.stderr.write("ref start %s >= primer_position %s" % (s.reference_start, primer_position))
           
            ## softmask the alignment if right primer start/end inside alignment
            if args.start:
                primer_position = p2[2]['start']
            else:
                primer_position = p2[2]['end']
            
            if args.strand:
                if s.is_reverse and s.reference_start >=  p2[2]['start'] and s.reference_start <=  p2[2]['end']:
                    trim(cigar, s, primer_position, 1)
            else:
                if s.reference_end > primer_position:
                    trim(cigar, s, primer_position, 1)
                else:
                    if args.verbose:
                        sys.stderr.write("ref end %s >= primer_position %s" % (s.reference_end, primer_position))
        except Exception as e:
            sys.stderr.write("problem %s" % (e,))
            pass

        if args.normalise:
            pair = "%s-%s-%d" % (p1[2]['Primer_ID'], p2[2]['Primer_ID'], s.is_reverse)
            counter[pair] += 1

            if counter[pair] > args.normalise:
                continue

        ## if the alignment starts before the end of the primer, trim to that position
#      trim(s, s.reference_start + 40, 0)
#      trim(s, s.reference_end - 40, 1)
#
#      outfile.write(s)
#   except Exception:
#      pass

        if not check_still_matching_bases(s):
             continue

        outfile.write(s)
        
    infile.close()
    outfile.close()   
    sys.stderr.close()
    if args.report:
        reportfh.close()


import argparse

parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
parser.add_argument('bedfile', help='BED file containing the amplicon scheme')
parser.add_argument('--normalise', type=int, help='Subsample to n coverage')
parser.add_argument('--report', type=str, help='Output report to file')
parser.add_argument('--start', action='store_true', help='Trim to start of primers instead of ends')
parser.add_argument('--strand', action='store_true', help='The strand is taken into account while doing the trimming ')
parser.add_argument('--verbose', action='store_true', help='Debug mode')

args = parser.parse_args()
go(args)
