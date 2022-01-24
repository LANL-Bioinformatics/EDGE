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
# from matplotlib import cm
# from matplotlib.colors import ListedColormap, LinearSegmentedColormap

colors = {
    'B.1.1.7': (0.4980392156862745, 0.788235294117647, 0.4980392156862745, 1.0),
    'B.1.351': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314, 1.0),
    'B.1.427': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.429': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253, 1.0),
    'P.1': (0.4, 0.4, 0.4, 1.0),
    'B.1.427/B.1.429': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.526': 'gold'}


def main():
    parser = argparse.ArgumentParser(description="Plot predictions along with clinical frequencies for a given VOC.")
    parser.add_argument('--voc', type=str, required=True, help="VOC to plot")
    parser.add_argument('--predictions', type=str, required=True, help="prediction csv")
    parser.add_argument('--clinical_data', type=str, required=True, help="add clinical data to plot")
    parser.add_argument('--aln_stats', type=str, nargs='+', help="basic alignment stats by samtools (including coverage info)")
    parser.add_argument('--sample_rna_levels', type=str, help="csv file containing sample RNA levels (copies/mL)")
    parser.add_argument('--fig_width', type=float, help="figure width")
    parser.add_argument('--fig_height', default=4, type=float, help="figure height")
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--outprefix', default="plot_clinical", help="add prefix to output figure names")
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    args = parser.parse_args()

    voc = args.voc
    column_types = {voc : float,
                    "lower_err" : float,
                    "upper_err" : float
                    }

    # read predictions
    df = pd.read_csv(args.predictions, sep=',', header=0, parse_dates=["Date"],
                     dtype=column_types)
    df["rolling_av"] = df[voc].rolling(3, center=True, min_periods=1).mean()
    conf_intervals = [list(df["lower_err"]), list(df["upper_err"])]

    # read RNA level per sample
    df["copies/mL"] = df["copies/mL"].astype("float")
    df["CT undiluted"] = df["CT undiluted"].astype("float")


    # read clinical data
    clinical_df = pd.read_csv(args.clinical_data,
                              sep=',',
                              header=0,
                              dtype={"Number of sequenced cases" : int},
                              parse_dates=["First day of week"])
    # add rows with zero cases where no clinical abundance is reported
    extra_rows = []
    for key, item in clinical_df.groupby("First day of week"):
        if voc not in item["Lineage"].values:
            row_dict = {"First day of week" : key,
                        "Lineage" : voc,
                        "Number of sequenced cases" : 0}
            extra_rows.append(row_dict)
    tmp_df = pd.DataFrame(extra_rows)
    clinical_df = clinical_df.append(tmp_df)

    for key, item in clinical_df.groupby("First day of week"):
        if voc not in item["Lineage"].values:
            print("WARNING: no abundance reported for {} at {}".format(voc, key))
    # compute relative VOC abundance per week
    total_cases = clinical_df.groupby("First day of week")["Number of sequenced cases"].sum()
    clinical_df = clinical_df.loc[clinical_df["Lineage"] == voc]
    clinical_df = clinical_df.sort_values(by="First day of week")
    clinical_df["total_cases"] = 0
    clinical_df["abundance(%)"] = 0.0
    for index, row in clinical_df.iterrows():
        n_cases = row["Number of sequenced cases"]
        abundance = n_cases / total_cases[row["First day of week"]] * 100
        clinical_df.loc[index, "abundance(%)"] = abundance
        clinical_df.loc[index, "total_cases"] = total_cases[row["First day of week"]]
    # print(clinical_df)
    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 12,
                         'legend.title_fontsize': 12}) # increase font size
    plt.figure()
    plt.plot_date(clinical_df["First day of week"],
                  clinical_df["abundance(%)"],
                  'g-',
                  marker='^',
                  ms=4,
                  label="Clinical samples")
    # print(voc)
    # print(clinical_df.head(100))
    # plot predictions
    predictions = df[voc]
    prediction_dates = pd.to_datetime(df["Date"])
    plt.bar(prediction_dates, predictions, label="Wastewater samples")
    ax = plt.gca()
    # change xticks: ticks and date at data points
    # plt.xticks(prediction_dates,
    #            [date.date() for date in prediction_dates],
    #            rotation=45,
    #            ha="right")
    # add tick markers at data points
    plt.plot_date(prediction_dates, [-2 for x in prediction_dates], marker=3, color='navy')
    plt.plot_date(prediction_dates, df["rolling_av"], '-', color='navy',
                  label="Wastewater rolling average (window=3)")
    # plt.errorbar(prediction_dates, predictions, yerr=conf_intervals)

    # add RNA levels per sample
    if args.sample_rna_levels:
        ax2 = ax.twinx()
        ax2.spines['right'].set_position(('axes', 1.0))
        df["copies/mL"] = df["copies/mL"].astype("float")
        # print(plot_df)
        ax2.plot_date(prediction_dates, df["copies/mL"], '+k', label="copies/mL")
        # df.plot(y="copies/mL", ax=ax2, x_compat=True, color="k", style="+", legend=False)
        ax2.set_ylabel("SARS-CoV-2 copies/mL sludge(+)")
    elif args.aln_stats:
        df["cov"] = 0.0
        for filename in args.aln_stats:
            # assumes that predictions are in directory named after sample
            sample_id = filename.split('/')[-1].split('.')[0]
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith("percentage of target"):
                        line = line.split('\t')
                        cov = float(line[1].rstrip('\n'))
                df.loc[df["ID"] == sample_id, "cov"] = cov
        ax2 = ax.twinx()
        ax2.spines['right'].set_position(('axes', 1.0))
        df["cov"] = df["cov"].astype("float")
        # print(plot_df)
        ax2.plot_date(prediction_dates, df["cov"], '+k', label="coverage")
        # df.plot(y="copies/mL", ax=ax2, x_compat=True, color="k", style="+", legend=False)
        ax2.set_ylabel("Percent genome with >20x coverage")


    if args.fig_width:
        plt.gcf().set_size_inches(args.fig_width, args.fig_height)
    plt.title(voc, fontsize=16)
    plt.grid(axis='y', which='major', alpha=0.2)
    plt.legend(loc='upper left')
    plt.ylabel("Abundance (%)")
    # sns.despine()
    plt.tight_layout()
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    ax.set_ylim(-5, 105)
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}_ylim100.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    ax.set_ylim(0, 1.05)
    for format in args.output_format.split(','):
        outfile = "{}/{}_{}_ylim1.{}".format(args.outdir, args.outprefix, voc, format)
        plt.savefig(outfile)

    # write data to tsv
    outfile = "{}/{}_{}.tsv".format(args.outdir, args.outprefix, voc)
    df.to_csv(outfile, sep='\t', index=False,
              columns=["Date", "ID", voc, "rolling_av"])

if __name__ == "__main__":
    sys.exit(main())
