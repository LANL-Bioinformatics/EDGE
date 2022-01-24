#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def main():
    parser = argparse.ArgumentParser(description="Compare predicted frequencies from kallisto and salmon.")
    parser.add_argument('--kallisto', type=str, nargs='+', help="kallisto prediction files")
    parser.add_argument('--salmon', type=str, nargs='+', help="salmon prediction files")
    parser.add_argument('--voc', dest='voc', type=str, required=True, help="comma-separated list of strains of interest")
    parser.add_argument('-o,--outdir', dest='outdir', required=True)
    parser.add_argument('-s,-suffix', dest='suffix', default="", help="add suffix to output figure names")
    parser.add_argument('-m', dest='min_ab', default=0, type=float, help="minimal abundance (any samples with true abundance below this threshold are skipped; any predictions below this threshold are considered absent)")
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    args = parser.parse_args()

    plt.rcParams.update({'font.size': 12,
                         'legend.fontsize': 12,
                         'legend.title_fontsize': 12,
                         'figure.titlesize': 14})

    voc_list = args.voc.split(',')
    output_formats = args.output_format.split(',')

    # read predictions
    kallisto_tups = read_predictions(args.kallisto, voc_list, args.min_ab)
    salmon_tups = read_predictions(args.salmon, voc_list, args.min_ab)

    # sort error tuples by voc frequency
    kallisto_tups.sort(key = lambda x : x[1])
    salmon_tups.sort(key = lambda x : x[1])
    voc_list = sorted(voc_list)

    # fix color per voc
    colormap = cm.get_cmap('Accent', len(voc_list))
    colors = {voc : colormap((i)/len(voc_list))
                for i, voc in enumerate(voc_list)}

    # plot kallisto and salmon prediction vs truth
    plt.figure()
    # freq_values = list(set([x[1] for x in kallisto_tups]))
    # kallisto_values = []
    # for freq in freq_values:
    #     estimates = [x[3] for x in kallisto_tups if x[1] == freq and x[3] > 0]
    #     av_est = 0 if len(estimates) == 0 else sum(estimates)/len(estimates)
    #     kallisto_values.append(av_est)
    freq_values = [x[1] for x in kallisto_tups]
    kallisto_values = [x[3] for x in kallisto_tups]
    plt.scatter(freq_values, kallisto_values, label="kallisto", alpha=0.7, s=20)
    # freq_values = list(set([x[1] for x in salmon_tups]))
    # salmon_values = []
    # for freq in freq_values:
    #     estimates = [x[3] for x in salmon_tups if x[1] == freq and x[3] > 0]
    #     av_est = 0 if len(estimates) == 0 else sum(estimates)/len(estimates)
    #     salmon_values.append(av_est)
    freq_values = [x[1] for x in salmon_tups]
    salmon_values = [x[3] for x in salmon_tups]
    plt.scatter(freq_values, salmon_values, label="salmon", alpha=0.7, s=20)
    plt.xlim(-5, 105)
    plt.ylim(-5, 105)
    plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("True VOC frequency (%)")
    plt.ylabel("Estimated VOC frequency (%)")
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_salmon_vs_truth_scatter{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    plt.xlim(-0.02, 5.02)
    plt.ylim(-0.02, 5.02)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_salmon_vs_truth_scatter_zoom{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.07, 150)
    plt.ylim(0.07, 150)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_salmon_vs_truth_scatter_loglog{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    # plot results per VOC
    for voc in voc_list:
        plt.figure()
        freq_values = [x[1] for x in kallisto_tups if x[0] == voc]
        kallisto_values = [x[3] for x in kallisto_tups if x[0] == voc]
        plt.scatter(freq_values, kallisto_values, label="kallisto", alpha=0.7, s=20)
        freq_values = [x[1] for x in salmon_tups if x[0] == voc]
        salmon_values = [x[3] for x in salmon_tups if x[0] == voc]
        plt.scatter(freq_values, salmon_values, label="salmon", alpha=0.7, s=20)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(0.07, 150)
        plt.ylim(0.07, 150)
        plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
        plt.legend()
        plt.title(voc)
        plt.grid(which="both", alpha=0.2)
        plt.xlabel("True VOC frequency (%)")
        plt.ylabel("Estimated VOC frequency (%)")
        plt.tight_layout()
        for format in output_formats:
            plt.savefig("{}/kallisto_salmon_vs_truth_scatter_loglog_{}{}.{}".format(args.outdir,
                                                             voc,
                                                             args.suffix,
                                                             format))

    # plot salmon vs kallisto predictions on a scatterplot
    plt.figure()
    for voc in voc_list:
        kallisto_values = [x[3] for x in kallisto_tups if x[0] == voc]
        salmon_values = [x[3] for x in salmon_tups if x[0] == voc]
        plt.scatter(kallisto_values, salmon_values, label=voc, alpha=0.7,
                    color=colors[voc], s=20)
    plt.xlim(-5, 105)
    plt.ylim(-5, 105)
    plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("Abundance estimate by kallisto (%)")
    plt.ylabel("Abundance estimate by salmon (%)")
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_vs_salmon_scatter{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    plt.xlim(-0.02, 5.02)
    plt.ylim(-0.02, 5.02)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_vs_salmon_scatter_zoom{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.07, 150)
    plt.ylim(0.07, 150)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/kallisto_vs_salmon_scatter_loglog{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    return


def read_predictions(filenames, voc_list, min_ab):
    err_tups = []
    for filename in filenames:
        dir_name = filename.split('/')[-2]
        voc_name = dir_name.split('_')[0]
        voc_freq = float(dir_name.split('_')[-1].lstrip('ab'))
        if voc_name not in voc_list:
            continue
        elif voc_freq < min_ab:
            continue
        with open(filename, 'r') as f:
            variant_found = False
            for line in f:
                if line[0] == "#":
                    continue
                [variant, tpm, ab, corrected_ab] = line.rstrip('\n').split('\t')
                if variant not in voc_list:
                    continue
                ab = float(ab)
                abs_err = abs(ab - voc_freq)
                if ab < min_ab:
                    continue
                if variant == voc_name:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, abs_err, ab))
            if not variant_found:
                # add zero estimate to error list
                err_tups.append((voc_name, voc_freq, voc_freq, 0))
    return err_tups


if __name__ == "__main__":
    sys.exit(main())
