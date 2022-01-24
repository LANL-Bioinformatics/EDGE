#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import statistics

def main():
    parser = argparse.ArgumentParser(description="Create sequencing depth plots.")
    parser.add_argument('depth_files', type=argparse.FileType('r'), nargs='+', help="samtools depth file(s) to be plotted")
    parser.add_argument('-o,--outfile', dest='outfile', type=str, default="plot.png", help="output figure name")
    parser.add_argument('-p1', type=int, default=0, help="start plotting at this position")
    parser.add_argument('-p2', type=int, default=-1, help="end plotting at this position")
    parser.add_argument('--samples', type=str, help="comma-separated list of sample names corresponding to input BAMs")
    parser.add_argument('--normalize', action="store_true", help="plot depth normalized by median")
    args = parser.parse_args()

    plot_file = args.outfile

    # read samtools depth file into dataframe
    df_dict = {}
    depth_files = []
    for file in args.depth_files:
        df = pd.read_csv(file,
                         sep='\t',
                         # engine='python',
                         names=["ref", "pos", "depth"]
                         )
        try:
            median_depth = statistics.median(df["depth"])
        except statistics.StatisticsError:
            continue
        # print("median depth {}: {}".format(file, median_depth))
        if median_depth == 0 and args.normalize:
            print("Median depth 0 for sample {} --- disabling normalization for this sample".format(file))
            median_depth = 1
        df["pos"] = df["pos"].astype('int')
        df["norm_depth"] = df["depth"] / median_depth
        df_dict[file] = df
        depth_files.append(file)

    pos1 = args.p1
    if args.p2 >= 0:
        assert args.p2 > args.p1
        pos2 = args.p2
    else:
        pos2 = max([max(df["pos"]) for df in df_dict.values()])

    # plot depth
    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 12,
                         'legend.title_fontsize': 12,
                         'figure.titlesize': 16})
    plt.figure(figsize=(15,3))
    for i, file in enumerate(depth_files):
        df = df_dict[file]
        if args.samples:
            samples = args.samples.split(',')
            name = samples[i]
        else:
            name = str(file).rstrip(".depth")
        selected_data = df.loc[(df["pos"] >= pos1) & (df["pos"] <= pos2)]
        pos_data = selected_data["pos"]
        if args.normalize:
            plt.plot(pos_data, selected_data["norm_depth"],
                     alpha=0.7,
                     label=name)
        else:
            plt.plot(pos_data, selected_data["depth"],
                     alpha=0.7,
                     label=name)

    if args.normalize:
        min_depth = min([min(df["norm_depth"]) for df in df_dict.values()])
        max_depth = max([max(df["norm_depth"]) for df in df_dict.values()])
        plt.ylabel("normalized depth")
    else:
        min_depth = min([min(df["depth"]) for df in df_dict.values()])
        max_depth = max([max(df["depth"]) for df in df_dict.values()])
        plt.ylabel("depth")
    # plt.ylim(min_depth, max_depth)
    # plt.ylim(-5000, 105000)
    if len(depth_files) <= 5:
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.xlabel("position")
    # plt.axvline(x=3925744, linestyle=':', color="red", label="oriC")
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.show()


    plt.ylim(0, 10000)
    outfile = ".".join(plot_file.split(".")[:-1] + ["ylim10000.png"])
    plt.savefig(outfile)

    if len(depth_files) > 1:
        # plot number of samples with coverage above threshold per position
        df = df_dict[depth_files[0]]
        df = df.set_index("pos")
        # print(df)
        combi_df = pd.DataFrame(index=df.index)
        combi_df[depth_files[0]] = df["depth"]
        for file in depth_files[1:]:
            df = df_dict[file]
            df = df.set_index("pos")["depth"]
            # print(df)
            combi_df = pd.concat([combi_df, df], axis=1)
        plt.figure()
        thresholds = [20, 100, 1000, 10000, 50000]
        # print(combi_df.index)
        # print(depth_files)
        for threshold in thresholds:
            counts = combi_df[combi_df > threshold].count(axis=1)
            plt.plot(combi_df.index,
                     counts,
                     label="threshold = {}x".format(threshold))
        plt.gcf().set_size_inches(15, 3)
        plt.xlabel("position")
        plt.ylabel("# samples with depth > threshold")
        plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
        plt.tight_layout()
        outfile = ".".join(plot_file.split(".")[:-1] + ["counts.png"])
        plt.savefig(outfile)

if __name__ == "__main__":
    sys.exit(main())
