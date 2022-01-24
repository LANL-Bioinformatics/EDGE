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
    parser.add_argument('-o, --outdir', dest='outdir', type=str, default='.', help="output directory")
    args = parser.parse_args()

    # Create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    # Select references per pango lineage
    full_df = read_metadata(args.metadata)
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


def read_metadata(metadata_file):
    """Read metadata from tsv into dataframe and filter for completeness"""
    df = pd.read_csv(metadata_file, sep='\t', header=0, dtype=str)
    # adjust date representation in dataframe
    df["date"] = df["Collection date"].str.replace('-XX','-01')
    df["date"] = pd.to_datetime(df.date, yearfirst=True)
    # remove samples wich have no pangolin lineage assigned (NaN or None)
    df = df.loc[df["Pango lineage"].notna()]
    df = df.loc[df["Pango lineage"] != "None"]
    # remove samples which are marked as incomplete or N-Content > 1%
    df = df.astype({"Is complete?" : 'bool',
                    "N-Content" : 'float'})
    df = df.loc[(df["Is complete?"] == True) & (df["N-Content"] <= 0.01)]
    return df


def select_ref_genomes(metadata_df, max_per_lineage, vcf_list, freq_list, min_aaf):
    """For every pangolin lineage, select exactly one sample."""
    # check which lineages are present
    lineages = metadata_df["Pango lineage"].unique()
    lineage_counts = metadata_df["Pango lineage"].value_counts()
    print("# lineages = {}".format(len(lineages)))
    # assign vcfs to lineages, assuming vcfs are in current directory and named after the corresponding lineage
    vcf_dict = {vcf.split('/')[-1].split('_')[0] : vcf for vcf in vcf_list}
    freq_dict = {fname.split('/')[-1].split('_')[0] : fname for fname in freq_list}
    # select samples for every lineage
    selection_ids = []
    for lin_id in lineages:
        samples = metadata_df.loc[metadata_df["Pango lineage"] == lin_id]
        samples = samples.sort_values(by=["N-Content", "Collection date"],
                                      ascending=[True, False])
        # read allele frequencies and extract sites with AAF >= minimal alt allele frequency
        try:
            allele_freq_file = freq_dict[lin_id]
        except KeyError as e:
            print("Skipping lineage {}, allele frequency info missing".format(lin_id))
            continue
        variant_positions = []
        variant_alleles = {}
        with open(allele_freq_file, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[0] == "CHROM":
                    continue
                ref_pos = int(line[1])
                ref_info = line[4]
                ref_allele, freq = ref_info.split(':')
                ref_allele_freq = float(freq)
                alt_alleles = line[5:]
                keep_alleles = []
                for alt_info in alt_alleles:
                    alt_allele, freq = alt_info.split(':')
                    alt_allele_freq = float(freq)
                    if alt_allele_freq > min_aaf:
                        keep_alleles.append(alt_allele)
                if len(keep_alleles) > 0:
                    variant_positions.append(ref_pos)
                    variant_alleles[ref_pos] = keep_alleles
        print("Filtering mutations for alt allele frequencies > {}".format(
                                                                    min_aaf))
        print("{} total # mutation sites kept = {}".format(
                lin_id, len(variant_positions)))
        # read vcf and process samples
        try:
            vcf_file = vcf_dict[lin_id]
        except KeyError as e:
            print("Skipping lineage {}, VCF info missing".format(lin_id))
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
                    if allele in variant_alleles[pos]:
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
                        metadata_df["Accession ID"].isin(selection_ids)]
    return selection_df


def filter_fasta(fasta_in, fasta_out, selection_df):
    """Filter fasta according to selected metadata"""
    keep_line = False
    selection_identifiers = list(selection_df["Virus name"].unique())
    with open(fasta_in, 'r') as f_in:
        with open(fasta_out, 'w') as f_out:
            for line in f_in:
                if line[0] == '>':
                    # sequence identifier
                    seq_id = line.rstrip('\n').lstrip('>').split('|')[0]
                    if len(selection_identifiers) == 0:
                        break
                    elif seq_id in selection_identifiers:
                        f_out.write(line)
                        keep_line = True
                        selection_identifiers.remove(seq_id)
                    else:
                        keep_line = False
                elif keep_line:
                    # nucleotide sequence
                    f_out.write(line)
    return


if __name__ == "__main__":
    sys.exit(main())
