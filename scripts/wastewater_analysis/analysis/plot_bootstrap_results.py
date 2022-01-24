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
# import seaborn as sns
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
    parser = argparse.ArgumentParser(description="Plot bootstrap results for a given VOC.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--run_info', type=str, nargs='+', default=[], help="run info files")
    parser.add_argument('--bootstraps', type=str, nargs='+', default=[], help="bootstrap hdf5 files")
    parser.add_argument('--confidence', type=float, default=95, help="confidence interval width")
    parser.add_argument('--metadata', type=str, required=True, help="metadata for samples")
    parser.add_argument('--metadata_ref', type=str, required=True, help="metadata for references")
    parser.add_argument('--voc', type=str, required=True, help="VOC to plot")
    parser.add_argument('--biobot', action='store_true', help="adapt plot settings suitable for biobot data")
    parser.add_argument('--max_bootstraps', type=int, help="maximum number of bootstrap experiments to be taken into account")
    parser.add_argument('--min_ab', default=0, type=float, help="set predictions <= min_ab to 0")
    parser.add_argument('--min_aln', type=int, default=0, help="remove samples with < min_aln aligned reads")
    parser.add_argument('--max_ct', type=float, default=34, help="remove samples with Ct > max_ct")
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
    if args.biobot:
        df["CT undiluted"] = df["Cq"].astype("float")
    else:
        df["CT undiluted"] = df["CT undiluted"].astype("float")
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
                assert len(row["CT undiluted"].values) == 1
                if row["CT undiluted"].values[0] <= args.max_ct:
                    if row["n_aligned"].item() >= args.min_aln:
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
    if args.biobot:
        plot_df = plot_df.sort_values(by=["State", "Sampling_location_(uploaded_from_processed_data)", "Date"])
    else:
        plot_df = plot_df.sort_values(by=["Date"])
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
    # print(states)
    # print(xticklabels)
    # print(plot_df)

    # fig, ax = plt.subplots(figsize=(5, 20))
    voc = voc_list[0]
    metadata_ref = pd.read_csv(args.metadata_ref, sep='\t', header=0, dtype=str)
    seqnames = list(metadata_ref.loc[metadata_ref["pangolin_lineage"] == voc,
                               "strain"])
    # plot_df["conf_interval"] = pd.Series([[0, 0]], dtype='object')
    plot_df["lower_err"] = 0
    plot_df["upper_err"] = 0
    # plot_df = plot_df.astype({"conf_interval": 'object'})
    # read bootstrap results
    for filename in args.bootstraps:
        print(filename)
        sample_id = filename.split('/')[-2]
        prediction = df.loc[df["ID"] == sample_id, voc].item()
        if prediction == ".":
            continue
        elif float(prediction) < args.min_ab:
            prediction = 0
            plot_df.loc[plot_df["ID"] == sample_id, voc] = 0
        if args.quick and prediction == 0:
            # plot_df.loc[plot_df["ID"] == sample_id, "conf_interval"] = [[0, 0]]
            continue
        f = h5py.File(filename, 'r')
        # ids = [id.decode("utf-8") for id in f["aux"]["ids"]]
        # eff_lengths = [int(l) for l in f["aux"]["eff_lengths"]]
        # voc_idx = ids.index(seqname)
        # voc_len = eff_lengths[voc_idx]
        bs_gen = f["bootstrap"]
        bootstrap_df = pd.DataFrame()
        bootstrap_df["id"] = [id.decode("utf-8") for id in f["aux"]["ids"]]
        bootstrap_df["eff_len"] = np.array(f["aux"]["eff_lengths"])
        bootstrap_results = []
        for boots_num, bs in enumerate(bs_gen):
            if args.max_bootstraps and boots_num >= args.max_bootstraps:
                break
            abundance = 0
            for seqname in seqnames:
                voc_idx = bootstrap_df.index[bootstrap_df["id"] == seqname]
                voc_len = bootstrap_df.loc[voc_idx, "eff_len"].item()
                # convert est_counts to % abundance (tpm / 10e4)
                bootstrap_df[bs] = np.array(f["bootstrap"][bs])
                normalization = np.array(bootstrap_df[bs] /
                                    bootstrap_df["eff_len"]).sum()
                voc_count = bootstrap_df.loc[voc_idx, bs].item()
                if voc_count > 0:
                    seq_ab = 100 * (voc_count / voc_len) / normalization
                    if seq_ab >= args.min_ab:
                        abundance += seq_ab
            assert abundance <= 100
            bootstrap_results.append(abundance)
        f.close()
        # compute confidence intervals
        ordered_bootstraps = sorted(bootstrap_results)
        lower = prediction - np.percentile(ordered_bootstraps,
                                           (100-args.confidence)/2)
        lower = min(max(lower, 0), prediction)
        upper = np.percentile(ordered_bootstraps,
                              args.confidence + ((100-args.confidence)/2)) \
                              - prediction
        upper = min(max(upper, 0), 100-prediction)
        assert lower >= 0
        assert upper >= 0
        assert prediction - lower <= 100
        assert prediction + upper <= 100
        # plot_df.loc[plot_df["ID"] == sample_id, "conf_interval"] = [[lower, upper]]
        plot_df.loc[plot_df["ID"] == sample_id, "lower_err"] = lower
        plot_df.loc[plot_df["ID"] == sample_id, "upper_err"] = upper

    # conf_intervals = list(plot_df["conf_interval"])
    # conf_intervals = list(zip(*conf_intervals))
    lower_err_list = list(plot_df["lower_err"])
    upper_err_list = list(plot_df["upper_err"])
    conf_intervals = [[lower_err_list, upper_err_list]]
    if args.verbose:
        print(plot_df)
    # plot predictions with error bars representing confidence intervals
    plot_df = plot_df.reset_index(drop=True)
    plt.figure()
    plt.rcParams.update({'font.size': 12}) # increase font size
    if args.bootstraps:
        ax = plot_df.plot(x="State",
                          y=voc,
                          kind="bar",
                          yerr=conf_intervals,
                          color=colors[voc],
                          legend=False,
                          capsize=2)
    else:
        ax = plot_df.plot(x="State",
                          y=voc,
                          kind="bar",
                          color=colors[voc],
                          legend=False,
                          capsize=2)
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels, fontsize=10, rotation=45, ha="right")
    # plt.xlabel("Sample location (state) and date")
    plt.xlabel("")
    plt.ylabel("Abundance (%)")

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
    plt.title(voc, fontsize=14)
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

if __name__ == "__main__":
    sys.exit(main())
