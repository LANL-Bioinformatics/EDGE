#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt

colors = {
    'B.1.1.7': (0.4980392156862745, 0.788235294117647, 0.4980392156862745, 1.0),
    'B.1.351': (0.9921568627450981, 0.7529411764705882, 0.5254901960784314, 1.0),
    'B.1.427': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.429': (0.7490196078431373, 0.3568627450980392, 0.09019607843137253, 1.0),
    'P.1': (0.4, 0.4, 0.4, 1.0),
    'B.1.427/B.1.429': (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0),
    'B.1.526': 'gold'}

def main():
    parser = argparse.ArgumentParser(description="Plot nucleotide diversity per VOC.")
    parser.add_argument('--site_pi_files', type=str, nargs='+', help="nucleotide diversity files (site-pi)")
    parser.add_argument('--allele_freq_files', type=str, nargs='+')
    parser.add_argument('--voc_names', type=str, nargs='+')
    parser.add_argument('--ref_size', required=True, type=int)
    parser.add_argument('--min_af', default=0, type=float)
    parser.add_argument('--outdir', default='sequence_diversity')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    # create output directory
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass

    diversity_dict = {}
    for i, file in enumerate(args.site_pi_files):
        voc = args.voc_names[i]
        pos_list = []
        pi_list = []
        with open(file, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[0] == "CHROM":
                    continue
                pos_list.append(int(line[1]))
                pi_list.append(float(line[2]))
        diversity_dict[voc] = [pos_list, pi_list]

    plt.rcParams.update({'font.size': 14,
                         'legend.fontsize': 12,
                         'legend.title_fontsize': 12,
                         'figure.titlesize': 16})

    plot_nuc_diversity(diversity_dict, args.ref_size, args.outdir)
    plot_nuc_diversity_subplots(diversity_dict, args.ref_size, args.outdir)

    allele_freq_dict = {}
    for i, file in enumerate(args.allele_freq_files):
        voc = args.voc_names[i]
        filtered = 0
        pos_list = []
        alt_allele_freq_list = []
        with open(file, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[0] == "CHROM":
                    continue
                ref_info = line[4]
                allele, freq = ref_info.split(':')
                ref_allele_freq = float(freq)
                alt_allele_freq = 1 - ref_allele_freq
                if alt_allele_freq > args.min_af:
                    if args.verbose:
                        print(line)
                    filtered += 1
                pos_list.append(int(line[1]))
                alt_allele_freq_list.append(alt_allele_freq)
        print("{} total # sites with alt allele frequency > {} = {}".format(
                voc, args.min_af, filtered))
        allele_freq_dict[voc] = [pos_list, alt_allele_freq_list]

    plot_allele_freq_subplots(allele_freq_dict, args.ref_size, args.outdir)

    return


def plot_nuc_diversity(diversity_dict, ref_size, outdir):
    """Plot nucleotide diversity per VOC"""
    outfile = outdir + "/nucleotide_diversity.png"
    plt.figure()
    for voc, diversity in diversity_dict.items():
        pos_list, pi_list = diversity
        pi_list_all = []
        pos0 = 0
        for i, pos1 in enumerate(pos_list):
            if pos1 == pos0:
                pi_list_all[-1] += pi_list[i]
            else:
                pi_list_all.extend([0]*(pos1-pos0-1))
                pi_list_all.append(pi_list[i])
            pos0 = pos1
        pi_list_all.extend([0]*(ref_size-pos1))
        # print(pi_list_all)
        assert len(pi_list_all) == ref_size
        plt.plot(range(1, ref_size+1), pi_list_all, label=voc,
                 color=colors[voc], alpha=0.7)
        print("{} done".format(voc))

    plt.xlabel("Reference position")
    plt.ylabel("Nucleotide diversity")
    plt.gcf().set_size_inches(15, 3)
    plt.legend(bbox_to_anchor=(1.1, 1))
    plt.tight_layout()
    plt.savefig(outfile)
    return

def plot_nuc_diversity_subplots(diversity_dict, ref_size, outdir):
    """Plot nucleotide diversity per VOC using subplots"""
    outfile = outdir + "/nucleotide_diversity_subplots.png"
    n_plots = len(diversity_dict.keys())
    fig, axs = plt.subplots(n_plots, 1, sharex=True, sharey=False)
    ax_idx = 0
    for voc, diversity in diversity_dict.items():
        pos_list, pi_list = diversity
        pi_list_all = []
        pos0 = 0
        for i, pos1 in enumerate(pos_list):
            if pos1 == pos0:
                pi_list_all[-1] += pi_list[i]
            else:
                pi_list_all.extend([0]*(pos1-pos0-1))
                pi_list_all.append(pi_list[i])
            pos0 = pos1
        pi_list_all.extend([0]*(ref_size-pos1))
        # print(pi_list_all)
        assert len(pi_list_all) == ref_size
        ax = axs[ax_idx]
        ax.plot(range(1, ref_size+1), pi_list_all, label=voc,
                      color=colors[voc], alpha=0.7)
        ax.set_ylim(0, 0.5)
        ax.set_ylabel("diversity")
        ax.set_title(voc)
        ax_idx += 1

    ax.set_xlabel("Reference position")
    plt.gcf().set_size_inches(15, 8)
    # plt.legend(bbox_to_anchor=(1.1, 1))
    plt.tight_layout()
    plt.savefig(outfile)
    return

def plot_allele_freq_subplots(allele_freq_dict, ref_size, outdir):
    """Plot nucleotide diversity per VOC using subplots"""
    outfile = outdir + "/allele_freq_subplots.png"
    n_plots = len(allele_freq_dict.keys())
    fig, axs = plt.subplots(n_plots, 1, sharex=True, sharey=False)
    ax_idx = 0
    for voc, freq_info in allele_freq_dict.items():
        pos_list, freq_list = freq_info
        freq_list_all = []
        pos0 = 0
        for i, pos1 in enumerate(pos_list):
            if pos1 == pos0:
                freq_list_all[-1] += freq_list[i]
            else:
                freq_list_all.extend([0]*(pos1-pos0-1))
                freq_list_all.append(freq_list[i])
            pos0 = pos1
        freq_list_all.extend([0]*(ref_size-pos1))
        # print(pi_list_all)
        assert len(freq_list_all) == ref_size
        ax = axs[ax_idx]
        ax.plot(range(1, ref_size+1), freq_list_all, label=voc,
                      color=colors[voc], alpha=0.7)
        ax.set_ylim(0, 1)
        ax.set_ylabel("AAF")
        ax.set_title(voc)
        ax_idx += 1

    ax.set_xlabel("Reference position")
    plt.gcf().set_size_inches(15, 8)
    # plt.legend(bbox_to_anchor=(1.1, 1))
    plt.tight_layout()
    plt.savefig(outfile)
    svg_outfile = outdir + "/allele_freq_subplots.svg"
    plt.savefig(svg_outfile)
    return


if __name__ == "__main__":
    sys.exit(main())
