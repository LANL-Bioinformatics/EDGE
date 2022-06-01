#/usr/bin/env python3

import argparse as ap
import json
import os

bin_dir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def setup_argparse():
    parser = ap.ArgumentParser(prog='recombinant_read_analysis',
                                description='''Script to do variants json to gff coversion ''')
    parser.add_argument(
        '--variantMutation', metavar='[FILE]', required=False, type=str, help=f"variant mutation json file [default: {bin_dir}/lineage_mutation.json]")
    parser.add_argument('--refacc', metavar='[STR]', required=False, type=str, default="NC_045512.2",
                        help='reference accession used in gff first column [default: NC_045512.2]')
    parser.add_argument('--output', metavar='[File]', required=False, type=str, default="variants_mutation.gff",
                        help="output variants mutaions in gff file [default: variants_mutation.gff]")
    argvs = parser.parse_args()
    if not argvs.variantMutation:
        argvs.variantMutation = os.path.join(os.path.dirname(bin_dir),"data", 'variant_mutation.json')
    return argvs

def main():
    argvs = setup_argparse()
    with open(argvs.variantMutation, 'r') as f:
        data = json.load(f)

    refacc=argvs.refacc
    of = open(argvs.output,"w")
    of.write("##gff-version 3\n")
    for k in data:
        for i, nt in enumerate(data[k].keys()):
            varID = f"{k}_{i + 1}"
            ref,pos,alt = nt.split(':')
            start = int(pos)
            end = int(pos)
            var = f"{ref}->{alt}"
            varAA = data[k][nt]
            if ref == 'del':
                var = f"{alt} nt deletion"
                end = int(pos) + int(alt) - 1
            if ref == "ins":
                var = f"{alt} insertion"
            of.write(f"{refacc}\tLANL\tSARS-CoV-2 variants\t{start}\t{end}\t100\t.\t.\tID={varID};Variant={k};mutation={var};AAchange={varAA}\n")

if __name__ == '__main__':
    main()

