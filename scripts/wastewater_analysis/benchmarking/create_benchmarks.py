#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import pandas as pd
from random import randint

from select_samples import filter_fasta, read_metadata, add_nonN_count


def main():
    parser = argparse.ArgumentParser(description="Create wastewater benchmarks.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, help="metadata tsv file for full sequence database")
    parser.add_argument('-s, --state', dest='state', type=str, default="Connecticut", help="sample location")
    parser.add_argument('-d, --date', dest='date', type=str, default="2021-02-11", help="sample date")
    parser.add_argument('-fr, --fasta_ref', dest='fasta_ref', required=True, type=str, help="fasta file representing full sequence database")
    parser.add_argument('-fv, --fasta_voc', dest='fasta_VOC', required=True, type=str, help="comma-separated list of fasta files for Variants Of Concern (VOC)")
    parser.add_argument('-o, --outdir', dest='outdir', required=True, type=str, help="output directory")
    parser.add_argument('--voc_perc', dest='voc_perc', required=True, type=str, help="comma-separated list of VOC frequencies (%) to be simulated")
    parser.add_argument('--total_cov', dest='total_cov', default=10000, type=int, help="total sequencing depth to be simulated")
    parser.add_argument('--data_exploration_only', action='store_true', help="exit after sequence selection")
    parser.add_argument('--spike_only', action='store_true', help="simulate reads for spike region only")
    args = parser.parse_args()

    # create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    VOC_frequencies = args.voc_perc.split(',')
    total_cov = args.total_cov
    VOC_files = args.fasta_VOC.split(',')
    VOC_names = [filepath.split('/')[-1] for filepath in VOC_files]
    exclude_list = [name.split('_')[0] for name in VOC_names]

    full_df = read_metadata(args.metadata)
    selection_df = select_benchmark_genomes(full_df, args.state, args.date,
                                            exclude_list)
    # filter fasta according to selection and write new fasta
    fasta_selection = args.outdir + "/sequences.fasta"
    filter_fasta(args.fasta_ref, fasta_selection, selection_df)
    print("Selected sequences written to {}".format(fasta_selection))
    # write corresponding metadata to tsv
    metadata_out = args.outdir + "/metadata.tsv"
    selection_df.to_csv(metadata_out, sep='\t', index=False)
    print("Metadata for selected sequences is in {}".format(metadata_out))

    if args.data_exploration_only:
        sys.exit()

    if args.spike_only:
        # trim sequences to select spike region
        print("\nTrimming genomes around spike region (21063--25884)")
        trimmed_selection = args.outdir + "/sequences.trimmed.fasta"
        subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(fasta_selection, trimmed_selection), shell=True)
        # also trim VOC sequences
        for filename in VOC_files:
            VOC_name = filename.rstrip('.fasta').split('/')[-1]
            trimmed_file = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
            subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(filename, trimmed_file), shell=True)
        fasta_selection = trimmed_selection
        print("\nSpike sequences ready\n")

    # simulate reads
    for VOC_freq in VOC_frequencies:
        # simulate reads for background sequences
        VOC_cov =  total_cov * float(VOC_freq)/100
        VOC_cov = round(VOC_cov, 2)
        background_cov = round((total_cov - VOC_cov) / len(selection_df.index), 2)
        print("Simulating reads from {} at {}x coverage".format(fasta_selection,
                                                                background_cov))
        subprocess.check_call("art_illumina -ss HS25 -i {0} -l 150 -f {1} -p -o {2}/background_{1}x -m 250 -s 10".format(fasta_selection, background_cov, args.outdir), shell=True)
        # simulate reads for VOC, then merge and shuffle
        for filename in VOC_files:
            VOC_name = filename.rstrip('.fasta').split('/')[-1]
            if args.spike_only:
                voc_fasta = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
            else:
                voc_fasta = filename
            print("Simulating reads from {} at {}x coverage".format(VOC_name,
                                                                    VOC_cov))
            subprocess.check_call("art_illumina -ss HS25 -i {0} -l 150 -f {1} -p -o {2}/{3}_{1}x -m 250 -s 10".format(voc_fasta, VOC_cov, args.outdir, VOC_name), shell=True)
            print("\nMerging fastqs...")
            subprocess.check_call("cat {0}/background_{1}x1.fq {0}/{2}_{3}x1.fq > {0}/tmp1.fq".format(args.outdir, background_cov, VOC_name, VOC_cov), shell=True)
            subprocess.check_call("cat {0}/background_{1}x2.fq {0}/{2}_{3}x2.fq > {0}/tmp2.fq".format(args.outdir, background_cov, VOC_name, VOC_cov), shell=True)
            print("Shuffling reads...")
            subprocess.check_call("shuffle.sh in={0}/tmp1.fq in2={0}/tmp2.fq out={0}/wwsim_{1}_ab{2}_1.fastq out2={0}/wwsim_{1}_ab{2}_2.fastq overwrite=t fastawrap=0".format(args.outdir, VOC_name, VOC_freq), shell=True)
        print("\nBenchmarks with a VOC frequency of {}% are ready!\n\n".format(VOC_freq))
    # clean up temporary files
    os.remove("{}/tmp1.fq".format(args.outdir))
    os.remove("{}/tmp2.fq".format(args.outdir))
    return


def select_benchmark_genomes(df, state, date, exclude_list):
    """Select genomes by location and date"""
    state_df = df.loc[df["division"] == state]
    selection_df = state_df.loc[state_df["date"] == date]
    print("\nLineage counts for {} on {}:".format(state, date))
    print(selection_df["pangolin_lineage"].value_counts())
    print("\nExcluding VOC lineages {} from selection\n".format(exclude_list))
    selection_df = selection_df.loc[
                        ~selection_df["pangolin_lineage"].isin(exclude_list)]
    # # show number of samples per date
    # samples_per_date = state_df["date"].value_counts().sort_index()
    # print("Samples per date:")
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #     print(samples_per_date)

    # # show lineages per date
    # grouped_mass_df = state_df.groupby(["date"])
    # print(grouped_mass_df.get_group(date)["strain", "pangolin_lineage"])
    # for key, item in grouped_mass_df:
    #     print(key, grouped_mass_df.get_group(key), "\n\n")
    return selection_df


if __name__ == "__main__":
    sys.exit(main())
