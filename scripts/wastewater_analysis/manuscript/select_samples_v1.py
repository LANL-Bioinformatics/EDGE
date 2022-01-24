#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import vcf


def main():
    parser = argparse.ArgumentParser(description="Build reference set consisting of a selection of samples per pangolin lineage.")
    parser.add_argument('--vcf', required=True, type=str, nargs='+', help="vcf files per lineage")
    parser.add_argument('--freq', required=True, type=str, nargs='+', help="allele frequency files per lineage")
    parser.add_argument('--min_aaf', default=0.5, type=float, help="minimal alternative allele frequency (AAF) to consider variation")
    parser.add_argument('--max_per_lineage', default=100, type=int, help="select at most k sequences per lineage")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, help="metadata tsv file for full sequence database")
    parser.add_argument('-f, --fasta', dest='fasta_in', type=str, help="fasta file representing full sequence database")
    parser.add_argument('-n, --nonN_counts', dest='nonN_counts', type=str, help="txt file with the number of nonambiguous characters per sequence")
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default='.', help="output directory")
    args = parser.parse_args()

    # Create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # Select references per pango lineage
    full_df = read_metadata(args.metadata, args.nonN_counts)
    selection_df = select_ref_genomes(full_df, args.max_per_lineage, args.vcf,
                                      args.freq, args.min_aaf)
    # Write metadata of selected samples to new tsv
    metadata_out = args.outdir + "/metadata.tsv"
    selection_df.to_csv(metadata_out, sep='\t', index=False)
    print("Metadata for selected sequences is in {}".format(metadata_out))
    # Filter fasta according to selection and write new fasta
    fasta_out = args.outdir + "/sequences.fasta"
    filter_fasta(args.fasta_in, fasta_out, selection_df)
    print("Selected sequences written to {}".format(fasta_out))
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


def select_ref_genomes(metadata_df, max_per_lineage, vcf_list, freq_list, min_aaf):
    """For every pangolin lineage, select exactly one sample."""
    # check which lineages are present
    lineages = metadata_df["pangolin_lineage"].unique()
    lineage_counts = metadata_df["pangolin_lineage"].value_counts()
    print("# lineages = {}".format(len(lineages)))
    # assign vcfs to lineages, assuming vcfs are in current directory and named after the corresponding lineage
    vcf_dict = {vcf.split('/')[-1].split('_')[0] : vcf for vcf in vcf_list}
    freq_dict = {fname.split('/')[-1].split('_')[0] : fname for fname in freq_list}
    # select samples for every lineage
    selection_ids = []
    for lin_id in lineages:
        samples = metadata_df.loc[metadata_df["pangolin_lineage"] == lin_id]
        samples = samples.sort_values(by=["nonN_count", "date"], ascending=False)
        # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
        try:
            allele_freq_file = freq_dict[lin_id]
        except KeyError as e:
            print("WARNING: skipping lineage {}, allele frequency info missing".format(lin_id))
            continue
        variant_positions = []
        with open(allele_freq_file, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[0] == "CHROM":
                    continue
                ref_info = line[4]
                ref_allele, freq = ref_info.split(':')
                ref_allele_freq = float(freq)
                alt_allele_freq = 1 - ref_allele_freq
                if alt_allele_freq > min_aaf:
                    variant_positions.append(int(line[1]))
        print("{} total # sites with alt allele frequency > {} = {}".format(
                lin_id, min_aaf, len(variant_positions)))
        # read vcf and process samples
        try:
            vcf_file = vcf_dict[lin_id]
        except KeyError as e:
            print("WARNING: skipping lineage {}, VCF info missing".format(lin_id))
            continue
        vcf_reader = vcf.Reader(open(vcf_file, 'rb'))
        samples = vcf_reader.samples
        sample_patterns = {sample : [] for sample in samples}
        for record in vcf_reader:
            if record.POS in variant_positions:
                alleles = [record.REF] + [str(x) for x in record.ALT]
                for sample in samples:
                    genotype = record.genotype(sample)['GT']
                    allele_idx = int(genotype[0])
                    allele = alleles[allele_idx]
                    sample_patterns[sample].append(allele)
        variation_seen = {pos : [] for pos in variant_positions}
        selection_count = 0
        if len(variant_positions) == 0:
            selection_ids.append(samples[0])
            selection_count += 1
        else:
            for sample in samples:
                select = False
                variation = sample_patterns[sample]
                for i, pos in enumerate(variant_positions):
                    allele = variation[i]
                    if allele not in variation_seen[pos]:
                        select = True
                        variation_seen[pos].append(allele)
                if select:
                    selection_ids.append(sample)
                    selection_count += 1
                    if selection_count == max_per_lineage:
                        break
        print("{} sequences selected for lineage {}".format(selection_count,
                                                            lin_id))
        if selection_count == 0:
            print("ERROR: no sequences selected for lineage {}".format(lin_id))
            sys.exit(1)

    print("{} sequences selected in total".format(len(selection_ids)))
    selection_df = metadata_df.loc[
                        metadata_df["gisaid_epi_isl"].isin(selection_ids)]
    return selection_df


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


def filter_fasta(fasta_in, fasta_out, selection_df):
    """Filter fasta according to selected metadata"""
    keep_line = False
    selection_identifiers = selection_df["strain"].unique()
    with open(fasta_in, 'r') as f_in:
        with open(fasta_out, 'w') as f_out:
            for line in f_in:
                if line[0] == '>':
                    # sequence identifier
                    seq_id = line.rstrip('\n').lstrip('>')
                    if seq_id in selection_identifiers:
                        f_out.write(line)
                        keep_line = True
                    else:
                        keep_line = False
                elif keep_line:
                    # nucleotide sequence
                    f_out.write(line)
    return


if __name__ == "__main__":
    sys.exit(main())
