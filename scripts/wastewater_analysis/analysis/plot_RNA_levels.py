#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import json
import datetime
import seaborn as sns
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Plot RNA levels.")
    parser.add_argument('--case_rate_info', type=str, help="csv file containing case rates and RNA levels (copies/mL)")
    parser.add_argument('--run_info', type=str, nargs='+', default=[], help="run info files")
    parser.add_argument('--aln_stats', type=str, nargs='+', required=True, help="basic alignment stats by samtools (including coverage info)")
    parser.add_argument('--aln_stats_spike', type=str, nargs='+', help="alignment stats by samtools for spike region")
    parser.add_argument('--metadata', type=str, required=True, help="metadata file")
    parser.add_argument('--spike', action='store_true', help="choose aln bins suitable for spike-only analysis")
    parser.add_argument('--biobot', action='store_true', help="choose rna level bins suitable for biobot data")
    parser.add_argument('-o,--outdir', dest='outdir', required=True)
    parser.add_argument('--outsuffix', type=str, default="")
    parser.add_argument('--outformat', type=str, default='png', help="comma-separated list of desired output formats")
    args = parser.parse_args()

    # read sample metadata
    df = pd.read_csv(args.metadata,
                     sep=',',
                     header=0,
                     dtype=str,
                     parse_dates=["Date"])
    df = df.dropna()

    # increase font size for all figures
    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 10,
                         'legend.title_fontsize': 10})

    # read run info per sample
    df["n_processed"] = -1
    df["n_aligned"] = 0
    df["n_unique"] = 0
    df["p_unique"] = 0.0
    for filename in args.run_info:
        # assumes that predictions are in directory named after sample
        sample_id = filename.split('/')[-2]
        with open(filename, 'r') as f:
            run_info = json.load(f)
            df.loc[df["ID"] == sample_id, "n_processed"] = int(run_info["n_processed"])
            df.loc[df["ID"] == sample_id, "n_aligned"] = int(run_info["n_pseudoaligned"])
            df.loc[df["ID"] == sample_id, "n_unique"] = int(run_info["n_unique"])
            df.loc[df["ID"] == sample_id, "p_unique"] = float(run_info["p_unique"])

    df = df.loc[df["n_processed"] >= 0]
    # print(df)

    # divide read counts into bins
    if args.spike:
        n_aligned_bins = [0, 50000, 100000, 200000, 1000000]
        bin_labels = ["0-50K", "50-100K", "100-200K", ">200K"]
    else:
        n_aligned_bins = [0, 500000, 1000000, 2000000, 10000000]
        bin_labels = ["0-0.5M", "0.5-1M", "1-2M", ">2M"]
    df["n_aligned_bin"] = pd.cut(df["n_aligned"],
                                 bins=n_aligned_bins,
                                 include_lowest=True,
                                 labels=bin_labels)

    # read RNA level per sample and divide into bins
    if "copies/mL" in df.columns:
        df["copies/mL"] = df["copies/mL"].astype("float")
        if args.biobot:
            rna_level_bins = [0, 500, 1000, 2000, 100000]
            bin_labels = ["0-500", "500-1000", "1000-2000", ">2000"]
        else:
            rna_level_bins = [0, 20000, 40000, 80000, 10000000]
            bin_labels = ["0-20K", "20-40K", "40-80K", ">80K"]
        df["copies/mL_bin"] = pd.cut(df["copies/mL"],
                                     bins=rna_level_bins,
                                     include_lowest=True,
                                     labels=bin_labels)

    # read CT values
    if "CT undiluted" in df.columns:
        df["CT"] = df["CT undiluted"].astype("float")
    elif "Cq" in df.columns:
        df["CT"] = df["Cq"].astype("float")
    mean_ct = df["CT"].mean()
    median_ct = df["CT"].median()
    print("Mean Ct value:", mean_ct)
    print("Median Ct value:", median_ct)

    # read alignment stats (coverage)
    df["cov"] = 0.0
    cov_threshold = -1
    for filename in args.aln_stats:
        # assumes that predictions are in directory named after sample
        sample_id = filename.split('/')[-1].split('.')[0]
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith("percentage of target"):
                    line = line.split('\t')
                    if cov_threshold == -1:
                        cov_threshold = line[0].split(" ")[-2]
                    else:
                        assert cov_threshold == line[0].split(" ")[-2]
                    cov = float(line[1].rstrip('\n'))
            df.loc[df["ID"] == sample_id, "cov"] = cov

    if args.aln_stats_spike:
        df["spike_cov"] = 0.0
        for filename in args.aln_stats_spike:
            # assumes that predictions are in directory named after sample
            sample_id = filename.split('/')[-1].split('.')[0]
            with open(filename, 'r') as f:
                for line in f:
                    if line.startswith("percentage of target"):
                        line = line.split('\t')
                        assert cov_threshold == line[0].split(" ")[-2]
                        cov = float(line[1].rstrip('\n'))
                df.loc[df["ID"] == sample_id, "spike_cov"] = cov

        plt.figure()
        sns.regplot(data=df, x="cov", y="spike_cov")
        plt.xlabel("Percent genome with >{}x coverage".format(cov_threshold))
        plt.ylabel("Percent spike with >{}x coverage".format(cov_threshold))
        plt.tight_layout()
        for fmt in args.outformat.split(','):
            plt.savefig("{}/cov_genome_vs_spike{}.{}".format(args.outdir,
                                                             args.outsuffix,
                                                             fmt))

    df = df.sort_values(by='Date')
    dates = df["Date"]

    # plot CT values vs coverage
    if "CT" in df.columns:
        ax = sns.relplot(data=df, x="CT", y="cov",
                   hue="n_aligned_bin", palette="rocket_r",
                   # fit_reg=False,
                   legend="brief",
                   facet_kws=dict(despine=False))
        # df.plot.scatter(x="copies/mL", y="cov")
        plt.xlabel("Ct")
        plt.ylabel("Percent genome with >{}x coverage".format(cov_threshold))
        plt.ylim(-5, 105)
        # plt.legend(bbox_to_anchor=(1, 0.3), borderaxespad=0., title="# reads aligned")
        ax._legend.set_bbox_to_anchor([0.97, 0.83])
        ax._legend.set_title("# reads aligned")
        ax._legend.set_frame_on(True)
        ax._legend.get_frame().set_linewidth(1)
        plt.gcf().set_size_inches(6, 5)
        plt.tight_layout()
        for fmt in args.outformat.split(','):
            plt.savefig("{}/CT_scatter_cov{}.{}".format(args.outdir,
                                                        args.outsuffix,
                                                        fmt))

    # plot readcount vs coverage
    plt.figure()
    # sns.regplot(data=df, x="n_aligned", y="cov")
    ax = sns.relplot(data=df, x="n_aligned", y="cov",
               hue="copies/mL_bin", palette="rocket_r",
               # fit_reg=False,
               legend="brief",
               facet_kws=dict(despine=False))
    plt.xlabel("# reads aligned")
    plt.ylabel("Percent genome with >{}x coverage".format(cov_threshold))
    plt.ylim(-5, 105)
    ax._legend.set_bbox_to_anchor([0.95, 0.25])
    ax._legend.set_title("SARS-CoV-2 copies/mL")
    ax._legend.set_frame_on(True)
    ax._legend.get_frame().set_linewidth(1)
    plt.gcf().set_size_inches(6, 5)
    plt.tight_layout()
    for fmt in args.outformat.split(','):
        plt.savefig("{}/aln_vs_cov{}.{}".format(args.outdir,
                                                args.outsuffix,
                                                fmt))

    outfile = args.outdir + "/raw_data_Ct_vs_cov.tsv"
    df.to_csv(outfile, sep='\t', index=False)


    if args.case_rate_info:
        # plot case rate versus wastewater RNA levels
        case_rate_df = pd.read_csv(args.case_rate_info, sep=',', header=0,
                                   dtype={"case rate per 100K" : float,
                                          "SARS-CoV-2 copies/ml sludge" : float
                                         },
                                   parse_dates=["Date"])
        case_rate_df["av_case_rate"] = case_rate_df["case rate per 100K"].rolling(7, center=True, min_periods=1).mean()
        case_rate_df["av_copies_ml"] = case_rate_df["SARS-CoV-2 copies/ml sludge"].rolling(7, center=True, min_periods=1).mean()

        plot_df = case_rate_df

        dates = plot_df["Date"]

        plt.rcParams.update({'font.size': 11})
        fig = plt.figure()
        ax1 = plt.gca()
        months = matplotlib.dates.MonthLocator()
        ax1.xaxis.set_major_locator(months)
        year_month = matplotlib.dates.DateFormatter('%Y-%m-%d')
        ax1.xaxis.set_major_formatter(year_month)
        ax2 = ax1.twinx()
        ax2.spines['right'].set_position(('axes', 1.0))
        ax2.xaxis.set_major_locator(months)
        ax2.xaxis.set_major_formatter(year_month)


        ax1.bar(dates,
                plot_df["SARS-CoV-2 copies/ml sludge"],
                label="SARS-CoV-2 copies/ml sludge",
                alpha=0.5)

        ax1.plot_date(dates, plot_df["av_copies_ml"], '-',
                      label="7-day rolling average",
                      color="navy")

        ax2.bar(dates,
                plot_df["case rate per 100K"],
                label="Case rate per 100K",
                color="lightgreen",
                alpha=0.5)

        ax2.plot_date(dates, plot_df["av_case_rate"], '-',
                    label="7-day rolling average",
                    color="green")

        ax1.set_zorder(1)
        ax1.set_frame_on(False)

        ax1.set_ylim(0,250000)
        ax2.set_ylim(0,175)

        plt.gcf().set_size_inches(10, 3)
        handles1, labels1 = ax1.get_legend_handles_labels()
        ax1.legend(handles1[::-1], labels1[::-1], loc="upper left")
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(handles2[::-1], labels2[::-1], loc="upper right")

        plt.tight_layout()
        for fmt in args.outformat.split(','):
            plt.savefig("{}/case_rates_rna_levels{}.{}".format(args.outdir,
                                                               args.outsuffix,
                                                               fmt))


    return




if __name__ == "__main__":
    sys.exit(main())
