#!/usr/bin/env python3
from Bio.Seq import Seq
import math
import re
import json
import os, sys

bin_dir = os.path.abspath(os.path.dirname(__file__))

class NT_AA_POSITION_INTERCHANGE():
    def __init__(self, ref_json=os.path.join(os.path.dirname(bin_dir),"data","SARS-CoV-2.json")):
        """ read the SARS-COV-2 json file"""
        if not os.path.exists(ref_json):
            print(f"Cannot find {ref_json}",file=sys.stderr)
            exit
        self.refdict={}
        with open(ref_json, 'r') as f:
            self.refdict = json.load(f)
    def convert_prot_nt(self,mut):
        """
        Converting protein to nucleotide position
        # mut: S:K467L   return: 22961
        """
        
        m = re.match(r'(\S+):(\D+)(\d+)(\D+)', mut)
        (gene, orig, pos, sub) = m.groups()
        
        return (int(pos)-1)*3 + self.refdict['genes'][gene]['coordinates']['from']

    def convert_nt_prot(self, mut):
        """
        Converting nucleotide changes to protein changes
        # mut: C:21618:G     return: S:T19R
        # mut: del:28362:9   return: N:E31-;N:R32-;N:S33-
        # mut: ins:22205:GAGCCAGAA    return: S:214EPE
        """
        m = re.match(r'(\S+):(\d+):(\S+)', mut)
        (ref, pos, alt) = m.groups()
        codon_num = ''
        return_str = ''
        for g in self.refdict['genes']:
            start =int(self.refdict['genes'][g]['coordinates']['from'])
            end = int(self.refdict['genes'][g]['coordinates']['to'])

            if int(pos) >= start and int(pos) <= end :
                codon_num = int((int(pos) - start) / 3) + 1
                codon_start = int(pos) - ((int(pos) - start)) % 3 - 1
                alt_codon = ''
                if ref == 'del':
                    if ((int(pos) - start + 1)) % 3 == 0:
                        codon = self.refdict['genome'][codon_start : codon_start + 3 + math.ceil((int(alt)-1)/3) * 3 ]
                        alt_codon = self.refdict['genome'][codon_start : codon_start+2] + self.refdict['genome'][codon_start+2+int(alt)]
                    if ((int(pos) - start + 1)) % 3 == 1:
                        codon = self.refdict['genome'][codon_start : codon_start + math.ceil(int(alt)/3) * 3]
                        alt_codon = self.refdict['genome'][codon_start+1+int(alt)] +  self.refdict['genome'][codon_start+2+int(alt)] + self.refdict['genome'][codon_start+3+int(alt)]
                    if ((int(pos) - start + 1)) % 3 == 2:
                        codon = self.refdict['genome'][codon_start : codon_start + 3 + math.ceil(int(alt)/3) * 3  ]
                        alt_codon = self.refdict['genome'][codon_start] + self.refdict['genome'][codon_start+1+int(alt)] + self.refdict['genome'][codon_start+2+int(alt)]  
                    aa = str(Seq(codon).translate())
                    alt_aa = str(Seq(alt_codon).translate())
                    split_aa=list()
                    for i in range(len(aa) - 1):
                        split_aa.append(g + ':' + aa[i] + str(codon_num) + '-')
                        codon_num = codon_num + 1
                    if aa[0] == alt_aa:
                        ## if aa doesnot change, move one aa ahead and make the deletion after 
                        split_aa.append(g + ':' + aa[-1] + str(codon_num) + '-')
                        split_aa = split_aa[1:]
                    else:
                        split_aa.append(g + ':' + aa[-1] + str(codon_num) + alt_aa)
                    return_str = ';'.join(split_aa)
                elif ref == 'ins':
                    if ((int(pos) - start + 1)) % 3 == 0:
                        codon = self.refdict['genome'][codon_start : codon_start + 3]
                        alt_codon = self.refdict['genome'][codon_start : codon_start+2] + alt  + self.refdict['genome'][codon_start+2: codon_start + 2 + (3 - (len(alt)+2)%3)]
                    if ((int(pos) - start + 1)) % 3 == 1:
                        codon = self.refdict['genome'][codon_start : codon_start + 3]
                        alt_codon = alt  + self.refdict['genome'][codon_start : codon_start + (3 - len(alt)%3)] if len(alt)%3 else alt
                    if ((int(pos) - start + 1)) % 3 == 2:
                        codon = self.refdict['genome'][codon_start : codon_start + 3]
                        alt_codon = self.refdict['genome'][codon_start] + alt  + self.refdict['genome'][codon_start+1: codon_start + 1 + (3 - (len(alt)+1)%3)]
                    codon_num = codon_num - 1
                    aa = str(Seq(codon).translate())
                    alt_aa = str(Seq(alt_codon).translate())
                    return_str = g + ":" + str(codon_num) + alt_aa if ((int(pos) - start + 1)) % 3 == 1 else g + ":" + aa + str(codon_num) + alt_aa
                else:
                    codon = self.refdict['genome'][codon_start : codon_start+3]
                    if ((int(pos) - start + 1)) % 3 == 0:
                        alt_codon = self.refdict['genome'][codon_start : codon_start+2] + alt
                    if ((int(pos) - start + 1)) % 3 == 1:
                        alt_codon = alt +  self.refdict['genome'][codon_start+1] +  self.refdict['genome'][codon_start+2] 
                    if ((int(pos) - start + 1)) % 3 == 2:
                        alt_codon = self.refdict['genome'][codon_start] + alt + self.refdict['genome'][codon_start+2] 
                    aa = str(Seq(codon).translate())
                    alt_aa = str(Seq(alt_codon).translate())
                    return_str = g + ":" + aa + str(codon_num) + alt_aa
                break

        return return_str