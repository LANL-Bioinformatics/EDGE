#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import glob


def main():
    parser = argparse.ArgumentParser(description="Preprocess reference collection: randomly select samples and write into individual files in lineage-specific directories.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, help="metadata tsv file for full sequence database")
    parser.add_argument('-n, --nonN_counts', dest='nonN_counts', type=str, help="txt file with the number of nonambiguous characters per sequence")
    parser.add_argument('-f, --fasta', dest='fasta_in', type=str, help="fasta file representing full sequence database")
    parser.add_argument('-k', dest='select_k', type=int, default=1000, help="randomly select 1000 sequences per lineage")
    parser.add_argument('--min_len', dest='min_len', type=int, help="minimal number of non-ambiguous nucleotides in sequence to be selected")
    parser.add_argument('--country', dest='country', type=str, help="only consider sequences found in specified country")
    parser.add_argument('--state', dest='state', type=str, help="only consider sequences found in specified state")
    parser.add_argument('--seed', dest='seed', default=0, type=int, help="random seed for sequence selection")
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default="seqs_per_lineage", help="output directory")
    args = parser.parse_args()

    # create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # read metadata
    metadata_df = read_metadata(args.metadata, args.nonN_counts)
    lineages = metadata_df["pangolin_lineage"].unique()

    # select sequences
    selection_dict = {}
    lineages_with_sequence = []
    for lin_id in lineages:
        # create lineage directory
        try:
            os.mkdir("{}/{}".format(args.outdir, lin_id))
        except FileExistsError:
            # empty existing directory
            old_files = glob.glob("{}/{}/*".format(args.outdir, lin_id))
            for f_trash in old_files:
                os.remove(f_trash)
        # filter for lineage, country and length
        samples = metadata_df.loc[metadata_df["pangolin_lineage"] == lin_id]
        if args.country:
            samples = samples.loc[samples["country"] == args.country]
        if args.state:
            samples = samples.loc[samples["division"] == args.state]
        if args.min_len:
            samples = samples.loc[samples["nonN_count"] > args.min_len]
        # randomly select sequences
        select_n = min(len(samples), args.select_k)
        selection = samples.sample(n=select_n, random_state=args.seed)
        if select_n == 0:
            print("WARNING: no sequences satisfying country and length restrictions for lineage {}".format(lin_id))
            continue
        elif select_n == 1:
            gisaid_id = selection["gisaid_epi_isl"].item()
            seq_name =  selection["strain"].item()
            selection_dict[seq_name] = (lin_id, gisaid_id)
        else:
            gisaid_ids = list(selection["gisaid_epi_isl"])
            seq_names = list(selection["strain"])
            for i, seq_name in enumerate(seq_names):
                gisaid_id =  gisaid_ids[i]
                selection_dict[seq_name] = (lin_id, gisaid_id)
        lineages_with_sequence.append(lin_id)

    print("{} sequences selected".format(len(selection_dict.keys())))

    # write sequences to separate files
    with open(args.fasta_in, 'r') as f_in:
        keep_line = False
        line_idx = 0
        selection_idx = 0
        for line in f_in:
            if line[0] == '>':
                # sequence identifier
                line_idx += 1
                if line_idx % 100000 == 0:
                    print("{} sequences from input fasta processed".format(line_idx))
                    print("{} sequences from selection found".format(selection_idx))
                seq_id = line.rstrip('\n').lstrip('>')
                try:
                    lin_id, gisaid_id = selection_dict[seq_id]
                    keep_line = True
                    selection_idx += 1
                    # print("keeping sequence {}".format(seq_id))
                except KeyError as e:
                    # item not found as sequence was not selected
                    # print("not keeping sequence {}".format(seq_id))
                    keep_line = False
            elif keep_line:
                # write nucleotide sequence
                outfile = "{}/{}/{}.fa".format(args.outdir, lin_id, gisaid_id)
                with open(outfile, 'w') as f_out:
                    f_out.write(">{}\n".format(seq_id))
                    f_out.write(line)
        print("{} sequences from input fasta processed".format(line_idx))
        print("{} sequences from selection found".format(selection_idx))

    # write lineages
    with open("{}/lineages.txt".format(args.outdir), 'w') as f:
        for lin_id in sorted(lineages_with_sequence):
            f.write("{}\n".format(lin_id))

    return


def read_metadata(metadata_file, nonN_count_file=None):
    """Read metadata from tsv into dataframe"""
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)
    # add field with number of N's in sequence
    if nonN_count_file:
        df = add_nonN_count(df, nonN_count_file)
    else:
        df["nonN_count"] = "."
    # adjust date representation in dataframe
    df["date"] = df["date"].str.replace('-XX','-01')
    df["date"] = pd.to_datetime(df.date, yearfirst=True)
    # remove samples wich have no pangolin lineage assigned (NaN or None)
    df = df[df["pangolin_lineage"].notna()]
    df = df[df["pangolin_lineage"] != "None"]
    return df


def add_nonN_count(df, nonN_count_file):
    """Count number of nonambiguous nucleotides per sequence and
    add counts to dataframe"""
    count_dict = {}
    with open(nonN_count_file, 'r') as f:
        for line in f:
            id, count = line.rstrip('\n').split('\t')
            count_dict[id] = int(count)
    assert len(df.index) == len(count_dict)
    count_list = [count_dict[id] for id in df["strain"]]
    df["nonN_count"] = count_list
    return df


if __name__ == "__main__":
    sys.exit(main())
