#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
import scipy.stats
import h5py
import json
import seaborn as sns

colors = {
    'B.1.1.7': (0.4980392156862745, 0.788235294117647, 0.4980392156862745, 1.0),
    'B.1.351': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314, 1.0),
    'B.1.427': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.429': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253, 1.0),
    'P.1': (0.4, 0.4, 0.4, 1.0),
    'B.1.427/B.1.429': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.526': 'gold'}

def main():
    parser = argparse.ArgumentParser(description="Plot biobot results for a given VOC.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--run_info', type=str, nargs='+', default=[], help="run info files")
    parser.add_argument('--confidence', type=float, default=95, help="confidence interval width")
    parser.add_argument('--metadata', type=str, required=True, help="metadata for samples")
    parser.add_argument('--metadata_ref', type=str, required=True, help="metadata for references")
    parser.add_argument('--gisaid_freqs', type=str, help="VOC frequencies per sampling date observed in GISAID")
    parser.add_argument('--voc', type=str, required=True, help="VOC to plot")
    parser.add_argument('--min_ab', default=0, type=float, help="set predictions <= min_ab to 0")
    parser.add_argument('--min_aln', type=int, default=0, help="remove samples with < min_aln aligned reads")
    parser.add_argument('--max_ct', type=float, default=34, help="remove samples with CT > max_ct")
    parser.add_argument('--min_pop', type=int, default=0, help="remove samples with population < min_pop")
    parser.add_argument('--quick', action='store_true', help="skip confidence interval computation for any predictions <= min_ab")
    parser.add_argument('--fig_width', type=float, help="figure width")
    parser.add_argument('--fig_height', default=4, type=float, help="figure height")
    parser.add_argument('--sep_states', action="store_true", help="separate different states by vertical lines in plot")
    parser.add_argument('--add_population', action="store_true", help="add population size on second y-axis")
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--outprefix', default="plot_bootstrap", help="add prefix to output figure names")
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    parser.add_argument('-v,--verbose', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.min_aln > 0 and not args.run_info:
        print("ERROR: can't filter for minimal number of aligned reads without run_info files")
        sys.exit(1)

    # read sample metadata
    df = pd.read_csv(args.metadata, sep=',', header=0, dtype=str, parse_dates=["Date"])

    voc_list = args.voc.split(',')
    for variant in voc_list:
        df[variant] = "."

    # read run info per sample
    df["n_reads"] = 0
    df["n_aligned"] = 0
    df["n_unique"] = 0
    for filename in args.run_info:
        # assumes that predictions are in directory named after sample
        sample_id = filename.split('/')[-2]
        with open(filename, 'r') as f:
            run_info = json.load(f)
            df.loc[df["ID"] == sample_id, "n_processed"] = run_info["n_processed"]
            df.loc[df["ID"] == sample_id, "n_aligned"] = run_info["n_pseudoaligned"]
            df.loc[df["ID"] == sample_id, "n_unique"] = run_info["n_unique"]

    # read prediction files per sample
    for filename in args.predictions:
        # assumes that predictions are in directory named after sample
        sample_id = filename.split('/')[-2]
        with open(filename, 'r') as f:
            freq_dict = {voc : 0 for voc in voc_list}
            for line in f:
                if line[0] == "#":
                    # header
                    continue
                [variant, tmp, freq, adj_freq] = line.strip().split('\t')
                if variant in voc_list:
                    freq_dict[variant] = float(freq)
            for voc, freq in freq_dict.items():
                row = df.loc[df["ID"] == sample_id]
                if float(row["Cq"]) <= args.max_ct:
                    if row["n_aligned"].item() >= args.min_aln:
                        if float(row["population"]) >= args.min_pop:
                            df.loc[df["ID"] == sample_id, voc] = freq



    if len(voc_list) > 1:
        joint_voc_list = voc_list
        joint_voc_name = '/'.join(joint_voc_list)
        df[joint_voc_name] = df[joint_voc_list].sum(axis=1)
        voc_list = [x for x in voc_list if x not in joint_voc_list]
        voc_list.append(joint_voc_name)
        colors[joint_voc_name] = colors[joint_voc_list[0]]

    voc_list = sorted(voc_list)
    colorlist = [colors[voc] for voc in voc_list]

    # plot composition per sample, sorted by state
    plot_df = df.loc[df[voc_list[0]] != "."]
    plot_df = plot_df.sort_values(by=["State", "Sampling_location_(uploaded_from_processed_data)", "Date"])
    states = plot_df["State"].unique()
    if args.sep_states:
        xticklabels = ["{0}, {1}".format(row["State"],
                      str(row["Date"]).split(" ")[0])
                      for index, row in plot_df.iterrows()]
        xtickinfo = ["{0}, {1}, {2}".format(row["State"],
                    str(row["Date"]).split(" ")[0],
                    row["Sampling_location_(uploaded_from_processed_data)"])
                    for index, row in plot_df.iterrows()]
    else:
        xticklabels = ["{1}".format(row["State"], str(row["Date"]).split(" ")[0])
                            for index, row in plot_df.iterrows()]

    # fig, ax = plt.subplots(figsize=(5, 20))
    voc = voc_list[0]
    metadata_ref = pd.read_csv(args.metadata_ref, sep='\t', header=0, dtype=str)
    seqnames = list(metadata_ref.loc[metadata_ref["pangolin_lineage"] == voc,
                               "strain"])

    if args.gisaid_freqs:
        gisaid_df = pd.read_csv(args.gisaid_freqs, sep='\t', header=0, parse_dates=["Date"])
        voc_freqs = []
        for index, row in plot_df.iterrows():
            freqs = gisaid_df.loc[(gisaid_df["State"] == row["State"]) &
                                 (gisaid_df["Date"] == row["Date"])][voc]
            freq = freqs.max()
            assert freq == freqs.min()
            voc_freqs.append(float(freq))
        gisaid_col = "gisaid-{}".format(voc)
        plot_df[gisaid_col] = voc_freqs

    # plot predictions with error bars representing confidence intervals
    plot_df = plot_df.reset_index(drop=True)
    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 12,
                         'legend.title_fontsize': 12}) # increase font size
    plt.figure()
    gisaid_col = "gisaid-{}".format(voc)
    ax = plot_df.plot(x="State",
                      y=[voc, gisaid_col],
                      kind="bar",
                      # color=colors[voc],
                      legend=True,
                      capsize=2)
    ax.legend(["Wastewater", "GISAID"])
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels, fontsize=10, rotation=45, ha="right")
    # plt.xlabel("Sample location (state) and date")
    plt.xlabel("")
    plt.ylabel("Abundance (%)")
    plt.ylim(0, 32)

    # add vertical lines to separate locations
    if args.sep_states:
        xticklocs, labels = plt.xticks()
        prev_location = ""
        prev_state = ""
        for i, label in enumerate(xtickinfo):
            [state, date, location] = label.split(",")
            if prev_state != "" and state != prev_state:
                line_loc = (xticklocs[i] + xticklocs[i-1]) / 2
                plt.axvline(x=line_loc, color='k', alpha=0.5, lw=0.5)
            elif prev_location != "" and location != prev_location:
                line_loc = (xticklocs[i] + xticklocs[i-1]) / 2
                plt.axvline(x=line_loc, color='k', alpha=0.5, lw=0.5, ls=(0, (5, 10)))
            prev_location = location
            prev_state = state

    # add population size per sample
    if args.add_population:
        ax2 = ax.twinx()
        ax2.spines['right'].set_position(('axes', 1.0))
        plot_df["population"] = plot_df["population"].astype("float")
        # print(plot_df)
        plot_df.plot(y="population", ax=ax2, x_compat=True, color="k", style="+", legend=False)
        ax2.set_ylabel("Population size per sample (+)")
        ax2.set_yscale("log")

    if args.fig_width:
        plt.gcf().set_size_inches(args.fig_width, args.fig_height)
    plt.title(voc, fontsize=16)
    ax.yaxis.grid(alpha=0.2)
    # ax.legend(ncol=len(voc_list), loc='upper right')
    plt.tight_layout()
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    ax.set_ylim(0, 105)
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}_ylim100.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    ax.set_ylim(0, 1.05)
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}_ylim1.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    # write predictions to file
    outfile = "{}/predictions_{}_m{}_a{}.csv".format(args.outdir, voc,
                                                     args.min_ab, args.min_aln)
    plot_df.to_csv(outfile)

    # plot population versus accuracy
    plot_df["population"] = plot_df["population"].astype("float")
    plot_df["diff"] = abs(plot_df[voc] - plot_df[gisaid_col])
    selection = plot_df.loc[plot_df[gisaid_col] > 1]
    plt.figure()
    plt.rcParams.update({'font.size': 12}) # increase font size
    selection.plot.scatter(x="population", y="diff")
    plt.tight_layout()
    plt.savefig(args.outdir + "/population_vs_accuracy_{}.png".format(voc))

if __name__ == "__main__":
    sys.exit(main())
