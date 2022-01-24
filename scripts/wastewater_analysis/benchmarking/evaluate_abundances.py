#!/usr/bin/env python3

import sys
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def main():
    parser = argparse.ArgumentParser(description="Evaluate predicted frequencies.")
    parser.add_argument('predictions', type=str, nargs='+', help="prediction files")
    parser.add_argument('--voc', dest='voc', type=str, required=True, help="comma-separated list of strains of interest")
    parser.add_argument('-o,--outdir', dest='outdir', required=True)
    parser.add_argument('-s,-suffix', dest='suffix', default="", help="add suffix to output figure names")
    parser.add_argument('-v,--verbose', dest='verbose', action='store_true')
    parser.add_argument('-m', dest='min_ab', default=0, type=float, help="minimal abundance (any samples with true abundance below this threshold are skipped; any predictions below this threshold are considered absent)")
    parser.add_argument('--no_plots', action='store_true')
    parser.add_argument('--joint_eval', dest='joint_eval', type=str, default="", help="comma-separated list of VOCs to be evaluated jointly (compare sum of estimates to sum of true frequencies)")
    parser.add_argument('--joint_average', action='store_true')
    parser.add_argument('--output_format', dest='output_format', default='png', help="comma-separated list of desired output formats")
    args = parser.parse_args()

    false_pos_count = 0
    false_neg_count = 0
    true_pos_count = 0
    true_neg_count = 0
    err_list = []
    variant_set = set()
    voc_list = args.voc.split(',')
    joint_voc_list = args.joint_eval.split(',')
    output_formats = args.output_format.split(',')

    # read predictions
    for filename in args.predictions:
        dir_name = filename.split('/')[-2]
        voc_name = dir_name.split('_')[0]
        voc_freq = float(dir_name.split('_')[-1].lstrip('ab'))
        if voc_name not in voc_list:
            continue
        elif voc_freq < args.min_ab:
            continue
        variant_set.add(voc_name)
        with open(filename, 'r') as f:
            variant_found = False
            err_tups = []
            positives = []
            for line in f:
                if line[0] == "#":
                    continue
                [variant, tpm, ab, corrected_ab] = line.rstrip('\n').split('\t')
                if variant not in voc_list:
                    continue
                ab = float(ab)
                abs_err = abs(ab - voc_freq)
                if ab < args.min_ab:
                    continue
                positives.append(variant)
                if variant == voc_name:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, abs_err, ab))
                elif variant in joint_voc_list and voc_name in joint_voc_list:
                    variant_found = True
                    err_tups.append((voc_name, voc_freq, abs_err, ab))
                else:
                    if args.joint_average and (variant in joint_voc_list
                                               or voc_name in joint_voc_list):
                        false_pos_count += 1/len(joint_voc_list)
                    else:
                        false_pos_count += 1
                    if args.verbose:
                        print("False positive: {} predicted at {}% in {}".format(
                                variant, ab, filename))
            if variant_found:
                if args.joint_average and (variant in joint_voc_list
                                           or voc_name in joint_voc_list):
                    true_pos_count += 1/len(joint_voc_list)
                else:
                    true_pos_count += 1
                if len(err_tups) == 1:
                    err_list.append(err_tups[0])
                else:
                    voc_name = err_tups[0][0]
                    voc_freq = err_tups[0][1]
                    ab = sum([x[3] for x in err_tups])
                    abs_err = abs(ab - voc_freq)
                    err_list.append((voc_name, voc_freq, abs_err, ab))
            else:
                if args.joint_average and voc_name in joint_voc_list:
                    false_neg_count += 1/len(joint_voc_list)
                else:
                    false_neg_count += 1
                if args.verbose:
                    print("VOC not found in {}".format(filename))
                # add zero estimate to error list?
                # err_list.append((voc_name, voc_freq, voc_freq, 0))
            for variant in voc_list:
                if variant not in positives and variant != voc_name:
                    # true negative
                    if args.joint_average and (variant in joint_voc_list
                                               or voc_name in joint_voc_list):
                        true_neg_count += 1/len(joint_voc_list)
                    else:
                        true_neg_count += 1
            true_neg_count += len([x for x in voc_list if
                                    x not in positives and x != voc_name ])


    # compute stats
    average_rel_err = sum([x[2]/x[1]*100 for x in err_list]) / len(err_list)
    average_rel_err_tp = (sum([x[2]/x[1]*100 for x in err_list if x[3] > 0])
                            / len(err_list))
    # print("average relative error: {}%".format(average_rel_err))
    print("average relative error of true positives: {}%".format(
                                                            average_rel_err_tp))
    print("total # true positives: {}".format(true_pos_count))
    print("total # true negatives: {}".format(true_neg_count))
    print("total # false positives: {}".format(false_pos_count))
    print("total # false negatives: {}".format(false_neg_count))

    fpr = false_pos_count / (false_pos_count + true_neg_count)
    fnr = false_neg_count / (false_neg_count + true_pos_count)
    recall = true_pos_count / (true_pos_count + false_neg_count)
    precision = true_pos_count / (true_pos_count + false_pos_count)
    print("FPR = {}".format(fpr))
    print("FNR = {}".format(fnr))
    print("Precision = {}".format(precision))
    print("Recall = {}\n".format(recall)) # sensitivity

    if args.no_plots:
        sys.exit()

    # sort error tuples by voc frequency
    err_list.sort(key = lambda x : x[1])
    variant_list = sorted(list(variant_set))

    # fix color per voc
    colormap = cm.get_cmap('Accent', len(variant_list))
    colors = {voc : colormap((i)/len(variant_list))
                for i, voc in enumerate(variant_list)}

    if args.joint_average:
        # compute average error for jointly evaluated VOCs
        if joint_voc_list != [""]:
            err_tups = [x for x in err_list if x[0] in joint_voc_list]
            new_err_list = [x for x in err_list if x[0] not in joint_voc_list]
            joint_voc_name = '/'.join(joint_voc_list)
            voc_freq = 0
            voc_freq_errs = []
            voc_freq_abs = []
            for tup in err_tups:
                if voc_freq > 0 and voc_freq != tup[1]:
                    av_err = sum(voc_freq_errs) / len(voc_freq_errs)
                    av_ab = sum(voc_freq_abs) / len(voc_freq_abs)
                    new_err_list.append((joint_voc_name, voc_freq, av_err, av_ab))
                    voc_freq_errs = []
                    voc_freq_abs = []
                voc_freq = tup[1]
                voc_freq_errs.append(tup[2])
                voc_freq_abs.append(tup[3])
            # finally add last tuple
            if voc_freq_errs:
                av_err = sum(voc_freq_errs) / len(voc_freq_errs)
                av_ab = sum(voc_freq_abs) / len(voc_freq_abs)
                new_err_list.append((joint_voc_name, voc_freq, av_err, av_ab))
            new_err_list.sort(key = lambda x : x[1])
            err_list = new_err_list
            variant_list = [voc for voc in variant_list if voc not in joint_voc_list]
            variant_list.append(joint_voc_name)
            colors[joint_voc_name] = colors[joint_voc_list[0]]

    plt.rcParams.update({'font.size': 14}) # increase font size
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc]
        err_values = [x[2]/x[1]*100 for x in err_list if x[0] == voc]
        plt.plot(freq_values, err_values, label=voc, color=colors[voc])
    plt.legend()
    plt.grid(which="both", alpha=0.2)
    plt.ylim(-5, 105)
    plt.xlabel("True VOC frequency (%)")
    plt.ylabel("Relative prediction error (%)")
    # plt.gcf().set_size_inches(4, 3)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot{}.{}".format(args.outdir,
                                                     args.suffix,
                                                     format))

    # also plot on log scale
    plt.xscale('log')
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_error_plot_logscale{}.{}".format(args.outdir,
                                                              args.suffix,
                                                              format))

    # plot true vs estimated frequencies on a scatterplot
    plt.figure()
    for voc in variant_list:
        freq_values = [x[1] for x in err_list if x[0] == voc]
        est_values = [x[3] for x in err_list if x[0] == voc]
        plt.scatter(freq_values, est_values, label=voc, alpha=0.7,
                    color=colors[voc], s=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.07, 150)
    plt.ylim(0.07, 150)
    plt.plot([0, 100], [0, 100], 'k-', lw=0.75)
    plt.legend(prop={'size': 12}) #ncol=len(variants_list),
    plt.grid(which="both", alpha=0.2)
    plt.xlabel("True VOC frequency (%)")
    plt.ylabel("Estimated VOC frequency (%)")
    # # Hide the right and top spines
    # ax = plt.gca()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    plt.tight_layout()
    for format in output_formats:
        plt.savefig("{}/freq_scatter_loglog{}.{}".format(args.outdir,
                                                         args.suffix,
                                                         format))

    return



if __name__ == "__main__":
    sys.exit(main())
