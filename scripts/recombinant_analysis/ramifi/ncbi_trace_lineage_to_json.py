#!/usr/bin/env python3

import argparse
import json
import re
from collections import defaultdict

import pandas as pd


def setup_argparse():
    parser = argparse.ArgumentParser(description='''Script to convert NCBI TRACK LineageDefinitions to JSON''')
    parser.add_argument('--out', metavar='[FILE]', required=False, type=str, default='lineage_mutation.json',
                        help='output file name [default: lineage_mutation.json]')
    parser.add_argument('--url', metavar='[URL]', required=False, type=str, default='https://ftp.ncbi.nlm.nih.gov/pub/ACTIV-TRACE/20220822-TRACE_LineageDefinitions-v42.1.txt',
                        help='txt file url from one of https://ftp.ncbi.nlm.nih.gov/pub/ACTIV-TRACE/ [default: https://ftp.ncbi.nlm.nih.gov/pub/ACTIV-TRACE/20220822-TRACE_LineageDefinitions-v42.1.txt]')
    argvs = parser.parse_args()
    return argvs

def main():
    argvs = setup_argparse()

    df=pd.read_table(argvs.url)

    lineage_dict = defaultdict(dict)
    df['variation'].replace(to_replace=r'(\d+)', value =r":\1:" , regex=True, inplace=True)
    for i in range(len(df)):
        lineage = df['id'][i]
        
        if len(df['nt_ref'][i]) - len(df['nt_alt'][i]) > 0 :
            nt = "del:" + str(df['nt_pos'][i] + 1) + ":" + str(len(df['nt_ref'][i]) - len(df['nt_alt'][i]))
        elif len(df['nt_ref'][i]) - len(df['nt_alt'][i]) < 0 :
            nt = "ins:" + str(df['nt_pos'][i] + 1) + ":" + df['nt_alt'][i][len(df['nt_ref'][i]):]
        else: 
            nt = df['nt_ref'][i] + ":" + str(df['nt_pos'][i]) + ":" + df['nt_alt'][i]
        aa = df['variation'][i]
        lineage_dict[lineage][nt]=aa

    with open(argvs.out, 'w', encoding ='utf8') as json_file:
        json.dump(lineage_dict, json_file, indent = 4)


if __name__ == '__main__':
	main()
