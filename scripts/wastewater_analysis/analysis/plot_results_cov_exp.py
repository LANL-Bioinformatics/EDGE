#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Plot predictions with error bars from coverage experiment.")
    parser.add_argument('--predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--voc', type=str, required=True)
    parser.add_argument('--repeats', type=int, required=True)
    parser.add_argument('-o,--outfile', dest='outfile', required=True)
    args = parser.parse_args()

    # read predictions
    data_dict = {}
    for file in args.predictions:
        with open(file, 'r') as f:
            p_dir = file.split('/')[-3]
            sample_id = file.split('/')[-4].lstrip("cov_exp_")
            prediction = -1
            for line in f:
                if line[0] == "#":
                    # header
                    continue
                [variant, tmp, freq, adj_freq] = line.strip().split('\t')
                if variant == args.voc:
                    prediction = float(freq)
                    break
            prediction = max(prediction, 0)
            p_dir = p_dir.lstrip("p")
            if p_dir not in data_dict:
                data_dict[p_dir] = [prediction]
            else:
                data_dict[p_dir].append(prediction)

    # store in dataframe
    df = pd.DataFrame()
    for p_dir, predictions in data_dict.items():
        if not len(predictions) == args.repeats:
            print("WARNING: number of predictions found for {} < {}".format(
                    p_dir, args.repeats))
        df[p_dir] = predictions
    # print(df)


    # increase font size for all figures
    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 10,
                         'legend.title_fontsize': 10})

    ax = sns.boxplot(data=df) #, color='grey'
    # plt.ylim(0, 100)
    # plt.title(sample_id)
    plt.xlabel("Genome coverage (%)")
    plt.ylabel("Predicted {} abundance (%)".format(args.voc))
    plt.gcf().set_size_inches(4, 5)
    plt.tight_layout()
    plt.savefig(args.outfile)
    return


if __name__ == "__main__":
    sys.exit(main())
