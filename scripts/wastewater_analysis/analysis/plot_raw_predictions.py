#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot abundances from file.")
    parser.add_argument('abundances', type=str, help="abundance file")
    parser.add_argument('--metadata', type=str, help="metadata file")
    parser.add_argument('-o', dest='outdir', type=str, default='.', help="write output here")
    parser.add_argument('--outprefix', type=str, default="")
    parser.add_argument('--outformat', type=str, default="png")
    args = parser.parse_args()

    abundance_df = pd.read_csv(args.abundances, sep='\t', header=0,
                               dtype = {'target_id': 'str',
                                        'length': 'int',
                                        'eff_length': 'int',
                                        'est_counts': 'float',
                                        'tpm': 'float'}
                              )
    abundance_df['freq'] = abundance_df['tpm'] / 10000
    # print(abundance_df)

    # plot distribution
    plt.figure()
    plt.gcf().set_size_inches(10, 3)
    plt.hist(abundance_df['freq'], bins=20)
    plt.xlim(0, 100)
    plt.xlabel("Abundance (%)")
    plt.ylabel("# predictions")
    plt.yscale('log')
    plt.tight_layout()
    for fmt in args.outformat.split(','):
        outfile = "{}/{}raw_freq_hist_0_100.{}".format(args.outdir,
                                                       args.outprefix,
                                                       fmt)
        plt.savefig(outfile)

    plt.xlim(1, 100)
    for fmt in args.outformat.split(','):
        outfile = "{}/{}raw_freq_hist_1_100.{}".format(args.outdir,
                                                       args.outprefix,
                                                       fmt)
        plt.savefig(outfile)

    plt.figure()
    plt.gcf().set_size_inches(10, 3)
    plt.hist(abundance_df['freq'], bins=100, range=(0, 1))
    plt.xlim(0, 1)
    plt.xlabel("Abundance (%)")
    plt.ylabel("# predictions")
    plt.yscale('log')
    plt.tight_layout()
    for fmt in args.outformat.split(','):
        outfile = "{}/{}raw_freq_hist_0_1.{}".format(args.outdir,
                                                      args.outprefix,
                                                      fmt)
        plt.savefig(outfile)
    return



if __name__ == "__main__":
    sys.exit(main())
