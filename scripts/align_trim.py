#!/usr/bin/env python3
from __future__ import print_function
import pysam
import sys
import csv
from copy import copy
from collections import defaultdict
import re

# consumesReference lookup for if a CIGAR operation consumes the reference sequence
# match, insertion, deletion, ref_skip, soft_clip, hard_clip, pad, equal
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]

def read_bed_file(fn):
    bedfile = []
    with open(fn) as csvfile:
        reader = csv.reader(csvfile, dialect='excel-tab')
        for row in reader:
        # ref start end primername
            bedrow = {}
            bedrow['Reference'] = row[0]
            primerID = row[3].replace('_ext','_R').replace('_lig','_L')
            bedrow['Primer_ID'] = re.sub('_alt\w+','_alt',primerID)
            #print(bedrow['Primer_ID'],file=sys.stderr)
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

def trim(s, start_pos, end):
    """Soft mask an alignment to fit within primer start/end sites.
    Parameters
    ----------
    s : pysam.AlignedSegment
        The aligned segment to mask
    start_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    """
    cigar = copy(s.cigartuples)
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
            sys.stderr.write(f"{s.query_name} ran out of cigar during soft masking - completely masked read will be ignored\n")
            break
 
        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if (consumesReference[flag]):
            if not end:
                 pos += length
            else:
                 pos -= length
                 
        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if (consumesQuery[flag]):
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= start_pos:
            break
        if end and pos <= start_pos:
            break
        #print >>sys.stderr, "pos:%s %s" % (pos, start_pos)
   
    ## primer start_pos in the deletion, it may have amplicon dropout and use other primer in the pool to amplified.
    if flag == 2 and length > 1:
        return 
    
    # calculate how many extra matches are needed in the CIGAR
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
    # softmask the left primer
    if not end:
        # update the position of the leftmost mappinng base
        s.pos = pos - extra
        if flag == 2:
            s.pos = pos
        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                s.pos += length

    if args.verbose:
        sys.stderr.write("New pos: %s" % (s.pos))
        
    if eaten:
        if end:
            cigar.append((4, eaten))
        else:
            cigar.insert(0, (4, eaten))
    # check the new CIGAR and replace the old one
    if args.verbose and (cigar[0][1] <= 0 or cigar[-1][1] <= 0):
        sys.stderr.write("invalid cigar operation created - possibly due to INDEL in primer\n")
    oldcigarstring = s.cigarstring
    s.cigartuples = cigar

    #print >>sys.stderr,  s.query_name, oldcigarstring[0:50], s.cigarstring[0:50]

def find_primer(bed, pos, direction, refID):
   # {'Amplicon_size': '1874', 'end': 7651, '#Region': 'region_4', 'start': 7633, 'Coords': '7633', "Sequence_(5-3')": 'GCTGGCCCGAAATATGGT', 'Primer_ID': '16_R'}
   from operator import itemgetter

   ref_alt_ID = re.sub(r'[^\w]','_',refID)
   search_term = re.compile (r'%s|%s' % (refID,ref_alt_ID) )
   if direction == '+':
       closest = min([(abs(p['start'] - pos), p['start'] - pos, p) for p in bed if p['direction'] == direction and search_term.search(p['Reference'])], key=itemgetter(0))
   else:
       closest = min([(abs(p['end'] - pos), p['end'] - pos, p) for p in bed if p['direction'] == direction and search_term.search(p['Reference'])], key=itemgetter(0))
   return closest


def go(args):
    if args.report:
        reportfh = open(args.report, "w")

    bed = read_bed_file(args.bedfile)

    counter = defaultdict(int)

    infile = pysam.AlignmentFile("-")
    outfile = pysam.AlignmentFile("-", "wh", template=infile)
    for s in infile:
        refname = s.reference_name
        #qname = s.query_name
        qlen = s.infer_read_length()
        amplicon_len = abs(s.template_length)
       
        ## logic - if alignment start site is _before_ but within X bases of
        ## a primer site, trim it off

        if s.is_unmapped:
            if args.verbose:
                sys.stderr.write("%s skipped as unmapped\n" % (s.query_name))
            continue

        if s.is_supplementary:
            if args.verbose:
                sys.stderr.write("%s skipped as supplementary\n" % (s.query_name))
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
        #if 'A00284:333:H755FDSX3:1:1147:28890:22842' in s.query_name:
            #print(f'{s.query_name} {s.cigarstring} {s.is_paired} {s.is_reverse} {amplicon_len} {qlen} {s.reference_start}  {s.reference_end} {p1} {p2}\n' , file=sys.stderr)
            #print(list(p1), file=sys.stderr)
            #print(list(p2), file=sys.stderr)
        try:
            ## softmask the alignment if left primer start/end inside alignment
            
            if args.start:
                primer_position = p1[2]['start'] 
            else:
                primer_position = p1[2]['end'] 
            
            if args.strand:
                if s.is_paired:
                    if amplicon_len > qlen:
                        if not s.is_reverse and s.reference_start >=  (p1[2]['start'] - args.offset) and s.reference_start <  p1[2]['end']:
                            trim(s, primer_position, 0)
                        ### The reads near to an amplified primers set.  ex primer_A_F  primer_A_R, try to trim both primers
                        if s.is_reverse and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1 and s.reference_start <= primer_position and s.reference_end > primer_position:
                            trim(s, primer_position, 0)
                    else:  # short amplicon,  reads length > amplicon size.  check the primer pair's name should be a set for trimming
                        if s.reference_start < primer_position and s.reference_end >= primer_position and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1:
                            trim(s, primer_position, 0)
                            # not proper pair, trim close primer ?
                        if not s.is_reverse and not s.is_proper_pair and s.reference_start >=  (p1[2]['start'] - args.offset) and s.reference_start <  p1[2]['end']:
                            trim(s, primer_position, 0)
                            
                else: ## unpaired reads
                    if not s.is_reverse and s.reference_start >=  (p1[2]['start'] - args.offset)  and s.reference_start <  p1[2]['end']:
                        trim(s, primer_position, 0)
                    ### The reads near to an amplified primers set.  ex primer_A_F  primer_A_R, try to trim both primers
                    if s.is_reverse and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1 and s.reference_start <= primer_position and s.reference_end > primer_position:
                        trim(s, primer_position, 0)
            else:
                # not the correct primer pair ex primer_A_RIGHT  primer_A_LEFT, filter the read
                if levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) >   1:
                    if s.reference_start >=  (p1[2]['start'] - args.offset) and s.reference_start <  p1[2]['end']:
                        trim(s, primer_position, 0)
                elif s.reference_start < primer_position:
                    trim(s, primer_position, 0)
                else:
                    if args.verbose:
                        sys.stderr.write("ref start %s >= primer_position %s" % (s.reference_start, primer_position))
           
            ## softmask the alignment if right primer start/end inside alignment
            if args.start:
                primer_position = p2[2]['start'] 
            else:
                primer_position = p2[2]['end'] 
            
            if args.strand:
                if s.is_paired:
                    if amplicon_len > qlen:
                        if s.is_reverse and s.reference_end >  p2[2]['end'] and s.reference_end <=  (p2[2]['start'] + args.offset):
                            trim(s, primer_position, 1)
                        ### The reads near to an amplified primers set.  ex primer_A_F  primer_A_R, try to trim both primers
                        if not s.is_reverse and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1 and s.reference_end > primer_position and s.reference_start < primer_position :
                            trim(s, primer_position, 1)
                    else:
                        if s.reference_end > primer_position and s.reference_start < primer_position and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1:
                            trim(s, primer_position, 1)
                        if s.is_reverse and not s.is_proper_pair and s.reference_end >  p2[2]['end'] and s.reference_end <=  (p2[2]['start'] + args.offset):
                            trim(s, primer_position, 1)
                            
                else: ## unpaired reads
                    if s.is_reverse and s.reference_end > p2[2]['end'] and s.reference_end <=  (p2[2]['start'] + args.offset):
                        trim(s, primer_position, 1)
                    ### The reads near to an amplified primers set.  ex primer_A_F  primer_A_R, try to trim both primers
                    if not s.is_reverse and levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) <= 1 and s.reference_end >= primer_position and s.reference_start < primer_position :
                        trim(s, primer_position, 1)
            else:
                if levenshtein_distance(p1[2]['Primer_ID'].replace('LEFT','L'), p2[2]['Primer_ID'].replace('RIGHT','R')) > 1 :
                    if s.reference_end > p2[2]['end'] and s.reference_end <=  (p2[2]['start'] + args.offset):
                        trim(s, primer_position, 1)
                elif s.reference_end > primer_position:
                    trim(s, primer_position, 1)
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

def levenshtein_distance(s, t, costs=(1, 1, 1)):
    """ 
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings 
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein 
        distance between the first i characters of s and the 
        first j characters of t
        
        costs: a tuple or a list with three integers (d, i, s)
               where d defines the costs for a deletion
                     i defines the costs for an insertion and
                     s defines the costs for a substitution
    """
    rows = len(s)+1
    cols = len(t)+1
    deletes, inserts, substitutes = costs
    dist = [[0 for x in range(cols)] for x in range(rows)]
    # source prefixes can be transformed into empty strings 
    # by deletions:
    for row in range(1, rows):
        dist[row][0] = row * deletes
    # target prefixes can be created from an empty source string
    # by inserting the characters
    for col in range(1, cols):
        dist[0][col] = col * inserts
        
    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0
            else:
                cost = substitutes
            dist[row][col] = min(dist[row-1][col] + deletes,
                                 dist[row][col-1] + inserts,
                                 dist[row-1][col-1] + cost) # substitution
    return dist[row][col]

import argparse

parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
parser.add_argument('bedfile', help='BED file containing the amplicon scheme')
parser.add_argument('--normalise', type=int, help='Subsample to n coverage')
parser.add_argument('--offset', type=int, default=1, help='extend 5 end primer bases [default: 1]')
parser.add_argument('--report', type=str, help='Output report to file')
parser.add_argument('--start', action='store_true', help='Trim to start of primers instead of ends')
parser.add_argument('--strand', action='store_true', help='The strand is taken into account while doing the trimming ')
parser.add_argument('--verbose', action='store_true', help='Debug mode')

args = parser.parse_args()
go(args)
