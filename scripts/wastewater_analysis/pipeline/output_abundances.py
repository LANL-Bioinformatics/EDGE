#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import json


def main():
    parser = argparse.ArgumentParser(description="Plot abundances from file.")
    parser.add_argument('abundances', type=str, help="abundance file")
    parser.add_argument('--metadata', type=str, required=True, help="metadata file")
    parser.add_argument('-m', dest='min_ab', type=float, default=0.1, help="noise level (% abundance); discard predictions below this threshold")
    parser.add_argument('--voc', dest='voc', type=str, help="comma-separated list of strains of interest, output abundance for these only")
    parser.add_argument('--voc_file', dest='voc_file', type=str, help="json file with VOC info (see example); alternative to using --voc")
    parser.add_argument('-o', dest='out_file', type=str, help="write output to tsv file")
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep='\t', header=0, dtype=str)

    if args.out_file:
        outfile = args.out_file
    else:
        outfilepath = args.abundances.split('/')
        outfilepath[-1] = "predictions.tsv"
        outfile = "/".join(outfilepath)

    # build dictionary specifying output lineages / variants
    if args.voc_file:
        with open(args.voc_file, 'r') as f:
            variants_dict = json.load(f)
    elif args.voc:
        variants_dict = {lineage : [lineage] for lineage in args.voc.split(',')}
    else:
        print("No VOC list (--voc) or VOC file (--voc_file) specified")
        print("Outputting abundance for ALL reference sequences")
        variants_dict = {}

    # print(variants_dict)

    lineage_to_variant = {}
    if args.voc or args.voc_file:
        abundance_dict = {variant : [0, 0] for variant in variants_dict.keys()}
        for variant, lineages in variants_dict.items():
            for lineage in lineages:
                assert lineage not in lineage_to_variant.keys()
                lineage_to_variant[lineage] = variant
    else:
        lineages = df["Pango lineage"].unique()
        lineage_to_variant = {lineage : lineage for lineage in lineages}
        abundance_dict = {lineage : [0, 0] for lineage in lineages}

    abundance_format = ""
    total_ab = 0
    with open(args.abundances, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[-1] == "tpm":
                # kallisto header
                abundance_format = "kallisto"
                continue
            elif line[-1] == "NumReads":
                # salmon header
                abundance_format = "salmon"
                continue
            if abundance_format == "":
                print("ERROR: abundance file format not recognized as kallisto or salmon")
                sys.exit(1)
            seqname = line[0].split('|')[0]
            lineage = df.loc[df["Virus name"] == seqname]["Pango lineage"]
            lineage = lineage.iloc[0]
            if abundance_format == "kallisto":
                tpm = float(line[-1])
            else: # salmon format
                tpm = float(line[-2])
            abundance = tpm / 10**6
            if abundance >= args.min_ab / 100:
                total_ab += abundance
                try:
                    variant = lineage_to_variant[lineage]
                    abundance_dict[variant][0] += tpm
                    abundance_dict[variant][1] += abundance
                except KeyError as e:
                    # lineage not required for output
                    pass

    # compute corrected abundances
    with open(outfile, 'w') as f:
        f.write("## evaluating {}\n".format(args.abundances))
        f.write("## {}\n".format(' '.join(sys.argv)))
        f.write("# variant\ttpm\tfreq(%)\tadj_freq(%)\n")
        for variant, values in abundance_dict.items():
            tpm, ab = values
            if total_ab > 0:
                corrected_ab = ab / total_ab
            else:
                corrected_ab = 0
            f.write("{}\t{:.0f}\t{:.2f}\t{:.2f}\n".format(
                    variant, tpm, ab * 100, corrected_ab * 100))

    return



if __name__ == "__main__":
    sys.exit(main())
