#!/usr/bin/env python3
import argparse as ap
import json
import os
import sys
import re
import shutil
import shlex
import subprocess
import time
import csv
from collections import defaultdict

import importlib_resources
import pandas as pd
import pysam

import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots
# standardize the logging output
import logging

toolname = os.path.basename(__file__)
bin_dir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()
sys.path.append(bin_dir)

try:
    from  __init__ import __version__
except:
    from  .__init__ import __version__

import translate


class SmartFormatter(ap.HelpFormatter):
    def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return ap.HelpFormatter._split_lines(self, text, width)


def setup_argparse():
    parser = ap.ArgumentParser(prog=toolname,
                                description='''Script to do recombinant read analysis''',
                                formatter_class=SmartFormatter )

    #parser.add_argument('--recombinant_variants', metavar='[list]', required=False, nargs=2, default=[
    #                    'Omicron', 'Delta'], help="list of two parents variants for recombinant checking [default: Omicron Delta]")
    parser.add_argument('--refacc', metavar='[STR]', required=False, type=str, default="NC_045512.2",
                        help='reference accession used in bam [default: NC_045512.2]')
    parser.add_argument('--minMixAF',metavar='[FLOAT]',required=False, type=float, default=0.2, help="minimum alleic frequency for checking mixed mutations on vcf [default:0.2]")
    parser.add_argument('--maxMixAF',metavar='[FLOAT]',required=False, type=float, default=0.8, help="maximum alleic frequency for checking mixed mutations on vcf [default:0.8]")
    parser.add_argument('--minMixed_n',metavar='[INT]',required=False, type=int, default=6, help="threshold of mixed mutations count for vcf.")
    parser.add_argument('--minReadCount',metavar='[INT]',required=False, type=int, default=10, help="threshold of read with variant count when no vcf provided.")
    parser.add_argument('--minRecombReadCount',metavar='[INT]',required=False, type=int, default=5, help="minimum recombinant reads count in a region [default:5]")
    parser.add_argument('--minRecombRate',metavar='[FLOAT]',required=False, type=float, default=0.1, help="minimum recombinant rate in a region [default:0.1]")
    parser.add_argument(
        '--lineageMutation', metavar='[FILE]', required=False, type=str, help=f"lineage mutation json file [default: variant_mutation.json]")
    parser.add_argument(
        '--variantMutation', metavar='[FILE]', required=False, type=str, help=f"variant mutation json file [default: lineage_mutation.json]")
    parser.add_argument(
        '--mutations_af_plot', action='store_true',  help="generate mutations_af_plot (when --vcf provided)")

    parser.add_argument('--verbose', action='store_true', 
                        help='Show more infomration in log')
    parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))

    inGrp = parser.add_argument_group('Input')
    inGrp.add_argument(
        '--bam', metavar='[FILE]', required=True, type=str, help="<Required> bam file")
    inGrp.add_argument(
        '--vcf',metavar='[File]', required=False, type=str, help="<Optional> vcf file which will infer the two parents of recombinant_variants")
    outGrp = parser.add_argument_group('Output')

    outGrp.add_argument('--tsv', metavar='[FILE]', required=False, type=str,
                        help='output file name [default: recombinant_reads.tsv]')
    outGrp.add_argument('--outbam', metavar='[File]', required=False, type=str,
                        help="output recombinant reads in bam file [default: recombinant_reads.bam]")

    EC19Grp = parser.add_argument_group('EDGE COVID-19 Options','options specific used for EDGE COVID-19')
    EC19Grp.add_argument('-eo', '--ec19_projdir',
                        metavar='[PATH]', required=False, type=str,  help="ec-19 project directory")
    EC19Grp.add_argument(
        '--igv', metavar='[PATH]', required=False, type=str,  help="igv.html relative path")
    EC19Grp.add_argument(
        '--igv_variants', action='store_true',  help="add variants igv track")

    argvs = parser.parse_args()

    if not argvs.lineageMutation:
        try:
            argvs.lineageMutation = importlib_resources.files(toolname).joinpath('data/lineage_mutation.json')
        except:
            argvs.lineageMutation = os.path.join(bin_dir,"data", 'lineage_mutation.json')
    if not argvs.variantMutation:
        try:
            argvs.variantMutation = importlib_resources.files(toolname).joinpath('data/variant_mutation.json')
        except:
            argvs.variantMutation = os.path.join(bin_dir,"data", 'variant_mutation.json')
    if not argvs.tsv:
        tsv_filename = "recombinant_reads.tsv"
        if argvs.ec19_projdir:
            argvs.tsv = os.path.join(argvs.ec19_projdir, 'ReadsBasedAnalysis',
                                     'readsMappingToRef', 'MixedInfectionAnalysis', tsv_filename)
        else:
            argvs.tsv = tsv_filename
    if argvs.ec19_projdir and not argvs.igv:
        argvs.igv = os.path.join('..', '..', '..','IGV', 'ref_tracks', 'igv.html')
    if not argvs.outbam:
        argvs.outbam = os.path.splitext(argvs.tsv)[0] + ".bam"
    return argvs


def load_var_mutation(file):
    with open(file, 'r') as f:
        v_data = json.load(f)
    delta_uniq_nt = list(
        set(v_data['Delta'].keys()) - set(v_data['Omicron'].keys())) if 'Delta' in v_data and 'Omicron' in v_data else list()
    omicron_uniq_nt = list(
        set(v_data['Omicron'].keys()) - set(v_data['Delta'].keys())) if 'Delta' in v_data and 'Omicron' in v_data else list()

    nt_to_aa = dict()
    nt_to_variant = dict()
    for k in v_data:
        nt_to_aa.update(v_data[k])
        for nt in v_data[k]:
            if nt not in nt_to_variant:
                nt_to_variant[nt] = [k]
            else:
                nt_to_variant[nt].append(k)

    return(delta_uniq_nt, omicron_uniq_nt, nt_to_variant, nt_to_aa)

def process_cmd(cmd, msg=None, stdout_log=None,  shell_bool=False, timeout = 60):
    if msg:
        logging.info(msg)

    time_out = time.time() + 60 * timeout   # 60 minutes in default
    logging.debug("CMD: %s" %(" ".join(cmd)) )
    proc = subprocess.Popen(shlex.split(" ".join(cmd)), shell=shell_bool, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    
    if stdout_log:
        f=open (stdout_log, "a") 
    while proc.poll() is None:  ## not finished the proc., if it finished will be not None.
        output = proc.stdout.readline()
        if output:
            stdout_msg=output.decode()
            sys.stdout.write(stdout_msg)
            if stdout_log:
                f.write(stdout_msg)
        # Wait a little before looping
        time.sleep(0.01)
        if time.time() > time_out:
            try:
                output = proc.communicate(timeout=5)[0]  # try 5 sec if still not end..., kill the proc and capture the final stream
            except subprocess.TimeoutExpired:
                proc.kill()
                output = proc.communicate()[0]
            if output:
                stdout_msg=output.decode()
                sys.stdout.write(stdout_msg)
                if stdout_log:
                    f.write(stdout_msg)
            break
    if stdout_log:
        f.close()
    proc.stdout.close()
    rc = proc.returncode
    if rc != 0: 
        logging.error("Failed %d %s" % (rc, " ".join(cmd)))
    
    return 0

def load_lineage_mutation(file):
    with open(file, 'r') as f:
        l_data = json.load(f)
    nt_to_lineage = dict()
    nt_to_aa = dict()
    for k in l_data:
        nt_to_aa.update(l_data[k])
        for nt in l_data[k]:
            if nt not in nt_to_lineage:
                nt_to_lineage[nt] = [k]
            else:
                nt_to_lineage[nt].append(k)
    return nt_to_lineage

def parse_vcf(argvs,nt_to_variants):
    filtered_nt_to_variants=dict()
    filtered_nt_to_variants_af=dict()
    filtered_nt_to_variants_dp=dict()
    parents_v=defaultdict(dict)

    mix_count=0
    mix_count_in_known_variants=dict()
    mutations_count=0
    with open(argvs.vcf,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            content = line.strip().split('\t')
            AFreq = 0
            depth = 0
            if "ALT_FREQ" in line:
                content2 = content[-1].split(':')
                AFreq = float(content2[-1])
            if 'DP=' in line:
                m = re.search(r'DP=(\d+)', line)
                depth=m[1]
            if 'DP4=' in line:
                m = re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', line)
                (ref_foward,ref_reverse,alt_forward,alt_reverse) = m.groups()
                AFreq = float((int(alt_forward) + int(alt_reverse))/int(depth))
            ref_bases = content[3]
            alt_bases = content[4]
            mutations_count += 1
            mut_list=[]
            if ',' in str(alt_bases) or (AFreq > argvs.minMixAF and AFreq < argvs.maxMixAF):
                mix_count += 1

            if ',' in alt_bases:
                alt_list = alt_bases.split(',')
                if len(ref_bases) == 1:
                    for alt in alt_list:
                        if len(alt) == 1:
                            mut = ref_bases + ":" + str(int(content[1])) + ":" + str(alt)
                        if len(alt) > 1:
                            mut = "ins" + ":" + str(int(content[1]) + 1) + ":" + str(alt[1:])
                        mut_list.append(mut)
                else: # len(ref_bases) > 1:
                    for alt in alt_list:
                        if len(alt) < len(ref_bases):
                            mut = 'del' + ":" + str(int(content[1]) + 1) + ":" + str(len(ref_bases)-len(alt))
                        if len(alt) > len(ref_bases):
                            mut = "ins" + ":" + str(int(content[1]) + 1) + ":" + str(alt.replace(ref_bases,''))
                        if len(alt) == len(ref_bases):
                            if alt[0] == ref_bases[0]:
                                alt=alt[1:]
                                ref_bases=ref_bases[1:]
                            mut = ref_bases + ":" + str(int(content[1]) + 1) + ":" + str(alt)
                        mut_list.append(mut)
            else:
                if len(ref_bases) - len(alt_bases) != 0:
                    ref_bases = 'del' if len(content[3]) - len(content[4]) > 0 else 'ins'
                    alt_bases = abs(len(content[3]) - len(content[4])) if ref_bases == 'del' else content[4].replace(content[3],'')
                    mut = ref_bases + ":" + str(int(content[1]) + 1) + ":" + str(alt_bases)
                else:
                    mut = ref_bases + ":" + str(int(content[1])) + ":" + str(alt_bases)
                mut_list.append(mut)

            for nt_v in mut_list:
                if nt_v in nt_to_variants:
                    if ',' in str(alt_bases) or (AFreq > argvs.minMixAF and AFreq < argvs.maxMixAF):
                        mix_count_in_known_variants[nt_v] = 1
                    ## scan first time find unique mutations variants
                    if len(nt_to_variants[nt_v]) == 1:
                        if ''.join(nt_to_variants[nt_v]) in parents_v:
                            parents_v[''.join(nt_to_variants[nt_v])]['uniq'] += 1 
                        else:
                            parents_v[''.join(nt_to_variants[nt_v])]['uniq'] = 1
                            parents_v[''.join(nt_to_variants[nt_v])]['all'] = 0
                            parents_v[''.join(nt_to_variants[nt_v])]['filtered'] = 0
                    if nt_v == mut_list[0]:
                        filtered_nt_to_variants[nt_v]=nt_to_variants[nt_v]
                        filtered_nt_to_variants_af[nt_v]=AFreq
                        filtered_nt_to_variants_dp[nt_v]=depth
                else:
                    filtered_nt_to_variants[nt_v]=[]
                    filtered_nt_to_variants_af[nt_v]=AFreq
                    filtered_nt_to_variants_dp[nt_v]=depth
    # use the first scan, uniq list to count all mutataions with the variants and variants between AF range
    for v in parents_v:
        AFavg = 0 
        for nt_v in filtered_nt_to_variants:
            if v in filtered_nt_to_variants[nt_v]:
                AFavg += filtered_nt_to_variants_af[nt_v]
                parents_v[v]['all'] += 1
                if filtered_nt_to_variants_af[nt_v] > argvs.minMixAF and filtered_nt_to_variants_af[nt_v] < argvs.maxMixAF:
                    parents_v[v]['filtered'] += 1
        AFavg = AFavg / parents_v[v]['all']
        parents_v[v]['avgAF'] = AFavg

    logging.debug(f"{mix_count}/{mutations_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.")
    logging.debug(f"{len(mix_count_in_known_variants)}/{mix_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF} and known variants.")
    ## sort by observerd variants counts
    parents_v = dict(sorted(parents_v.items(), key=lambda item: item[1]['all'], reverse=True))
    logging.debug(f"All probable parents, mutation count: {parents_v}")
    logging.debug(f"All probable parents, mutation count with AF {argvs.minMixAF}-{argvs.maxMixAF}: {parents_v}")
    if len(parents_v.keys())<2:
        logging.error(f'no two parents variants detected. {parents_v}')
        #sys.exit(1)
    if mix_count < argvs.minMixed_n:
        logging.error(f'count of mixed mutations with AF between {argvs.minMixAF} and {argvs.maxMixAF} is less than {argvs.minMixed_n}.')
        sys.exit(1)
    if len(parents_v.keys()) > 0:
        parent1 = list(parents_v.keys())[0]
        logging.info(f"Parent 1 ({parent1}) mean AF: {parents_v[parent1]['avgAF']}")
        if len(parents_v.keys()) > 1:
            parent2 = list(parents_v.keys())[1]
            logging.info(f"Parent 2 ({parent2}) mean AF: {parents_v[parent2]['avgAF']}")

    return parents_v, filtered_nt_to_variants, filtered_nt_to_variants_af, filtered_nt_to_variants_dp

def find_read_with_variants(nt_to_variant, argvs):
    if not os.path.exists(argvs.bam + ".bai"):
        pysam.index(argvs.bam)
    samfile = pysam.AlignmentFile(argvs.bam, "rb")
    ## contig, mapped, unmapped, total
    idx_stats = samfile.get_index_statistics()
    stats = {'total': idx_stats[0].total,
             'mapped': idx_stats[0].mapped, 'unmapped': idx_stats[0].unmapped}
    logging.info(f"Total Reads: {stats['total']}")
    logging.info(f"Mapped Reads: {stats['mapped']}")
    mutation_reads = defaultdict(lambda: defaultdict(list))
    reads_coords = defaultdict(dict)

    logging.info("Finding reads with variant mutations")
    for i in nt_to_variant:
        ref, pos, alt = i.split(':')
        pos = int(pos)
        check_pos = pos - 1
        if (ref == 'del' or ref == 'ins'):
            check_pos = pos - 2
        # pos is 1-index and substract 1 to get zero index for pileup
        for pileupcolumn in samfile.pileup(contig=argvs.refacc, start=pos-1, stop=pos):              
            if pileupcolumn.pos == check_pos:
                #print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                for pileupread in pileupcolumn.pileups:
                    if ref == 'del' or ref == 'ins':
                        if pileupread.indel < 0 and ref == 'del':
                            #print('\tpos %s %s indel base in read %s = %s %s' % (pos,pileupcolumn.pos, pileupread.alignment.query_name, pileupread.indel , pileupread.query_position))
                            if (int(alt) + int(pileupread.indel)) == 0:
                                for mut in nt_to_variant[i]:
                                    if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                        mutation_reads[pileupread.alignment.query_name][pos].extend(
                                            nt_to_variant[i])
                            else:
                                if 'other' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                    mutation_reads[pileupread.alignment.query_name][pos].extend(['other'])
                        if pileupread.indel > 0 and ref == 'ins':
                            if (len(alt) - int(pileupread.indel)) == 0:
                                for mut in nt_to_variant[i]:
                                    if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                        mutation_reads[pileupread.alignment.query_name][pos].extend(
                                            nt_to_variant[i])
                            else:
                                if 'other' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                    mutation_reads[pileupread.alignment.query_name][pos].extend(['other'])
                        if pileupread.indel == 0:
                            if 'rev' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['rev'] + nt_to_variant[i])
                    else:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if alt == pileupread.alignment.query_sequence[pileupread.query_position]:
                                for mut in nt_to_variant[i]:
                                    if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                        mutation_reads[pileupread.alignment.query_name][pos].extend([mut])
                            elif ref == pileupread.alignment.query_sequence[pileupread.query_position]:
                                if 'rev' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                    mutation_reads[pileupread.alignment.query_name][pos].extend(
                                        ['rev'] + nt_to_variant[i])
                            else:
                                other_mut = ref + ":" + \
                                    str(pos) + ":" + \
                                    pileupread.alignment.query_sequence[pileupread.query_position]

                                if other_mut in nt_to_variant:
                                    for mut in nt_to_variant[other_mut]:
                                        if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                            #print(nt_to_variant[other_mut])
                                            mutation_reads[pileupread.alignment.query_name][pos].extend([mut])
                                elif 'other' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                    mutation_reads[pileupread.alignment.query_name][pos].extend(['other'])
                    if pileupread.alignment.query_name in mutation_reads:
                        if 'start' in reads_coords[pileupread.alignment.query_name]:
                            if pileupread.alignment.reference_start < reads_coords[pileupread.alignment.query_name]['start']:
                                reads_coords[pileupread.alignment.query_name]['start'] = pileupread.alignment.reference_start
                        else:
                            reads_coords[pileupread.alignment.query_name]['start'] = pileupread.alignment.reference_start
                        if 'end' in reads_coords[pileupread.alignment.query_name]:
                            if pileupread.alignment.reference_end > reads_coords[pileupread.alignment.query_name]['end']:
                                reads_coords[pileupread.alignment.query_name]['end'] = pileupread.alignment.reference_end
                        else:
                            reads_coords[pileupread.alignment.query_name]['end'] = pileupread.alignment.reference_end
    samfile.close()
    return mutation_reads, reads_coords, stats

def find_parents_by_mutation_reads(mutation_reads, argvs):
    parents_v=defaultdict(dict)
    for name in mutation_reads:
        for pos in mutation_reads[name]:
            if 'rev' not in mutation_reads[name][pos] and 'unknown' not in mutation_reads[name][pos] and len(mutation_reads[name][pos]) == 1:
                if ''.join(mutation_reads[name][pos]) in parents_v:
                    parents_v[''.join(mutation_reads[name][pos])]['uniq'] += 1
                
                else:
                    parents_v[''.join(mutation_reads[name][pos])]['uniq'] = 0
                    parents_v[''.join(mutation_reads[name][pos])]['all'] = 0
                
    for v in parents_v:
        for name in mutation_reads:
            for pos in mutation_reads[name]:
                if 'rev' not in mutation_reads[name][pos] and v in mutation_reads[name][pos]:
                    parents_v[v]['all'] += 1

    parents_v = dict(sorted(parents_v.items(), key=lambda item: item[1]['all'], reverse=True))
    if (parents_v[list(parents_v.keys())[0]]['all'] < argvs.minReadCount or parents_v[list(parents_v.keys())[1]]['all'] < argvs.minReadCount):
        logging.error(f'Not have enough reads ({argvs.minReadCount}) to support both parents variants: {parents_v}')
        sys.exit(1)

    if len(parents_v.keys())<2:
        logging.error(f'no two parents variants detected. {parents_v}')
        #sys.exit(1)

    
    return parents_v

def find_recomb(mutation_reads, reads_coords, two_parents_list, reads_stats, argvs):
    logging.info("Finding recombinant reads")
    logging.info(f"Parents: {','.join(two_parents_list)} ")
    ## Can expand the list to check other recombinant???
    recomb1_reads = []
    recomb2_reads = []
    recombx_reads = []  # recombinant happend 2x more in a read
    parent1_reads = []
    parent2_reads = []
    reads_has_two_and_more_mutations_count = 0
    cross_region=defaultdict(lambda: defaultdict(list))
    parent_pos_count = dict()
    total_dist = 0
    for read in mutation_reads:
        list_for_check_recomb=[]
        list_for_check_recomb_pos=[]
        var1_count = 0
        var2_count = 0
        if len(mutation_reads[read]) > 1:
            reads_has_two_and_more_mutations_count += 1
        for k, v in sorted(mutation_reads[read].items()):
            var1 = two_parents_list[0]
            var2 = two_parents_list[1]
            if 'rev' in v:
                if var1 in v and var2 not in v:
                    var2_count += 1
                    list_for_check_recomb.append("P2")
                    list_for_check_recomb_pos.append(k)
                if var2 in v and var1 not in v:
                    var1_count += 1
                    list_for_check_recomb.append("P1")
                    list_for_check_recomb_pos.append(k)
            else:
                if var1 in v and var2 not in v:
                    var1_count += 1
                    list_for_check_recomb.append("P1")
                    list_for_check_recomb_pos.append(k)
                if var2 in v and var1 not in v:
                    var2_count += 1
                    list_for_check_recomb.append("P2")
                    list_for_check_recomb_pos.append(k)
        if len(list_for_check_recomb) > 1:
            if 'P1' not in list_for_check_recomb:
                for i in range(reads_coords[read]['start']+1,reads_coords[read]['end']+1):
                    parent_pos_count[i] = parent_pos_count[i] + 1 if i in parent_pos_count else 1 
                parent2_reads.append(read)
            if 'P2' not in list_for_check_recomb:
                for i in range(reads_coords[read]['start']+1,reads_coords[read]['end']+1):
                    parent_pos_count[i] = parent_pos_count[i] + 1 if i in parent_pos_count else 1
                parent1_reads.append(read)
            if 'P1' in list_for_check_recomb and 'P2' in list_for_check_recomb:
                if 'P1' == list_for_check_recomb[0] and (list_for_check_recomb == sorted(list_for_check_recomb) ):
                    recomb1_reads.append(read)
                    recomb1_start=list_for_check_recomb_pos[0]
                    recomb1_end=list_for_check_recomb_pos[1]
                    for r in range(1,len(list_for_check_recomb)):
                        if list_for_check_recomb[r] == 'P1':
                            recomb1_start = list_for_check_recomb_pos[r]
                        else:
                            recomb1_end = list_for_check_recomb_pos[r]
                            break
                    cross_region[f"{recomb1_start}-{recomb1_end}"]['recomb1'].append(read)
                    dist = int(recomb1_end) - int(recomb1_start) + 1
                    total_dist = dist + total_dist
                elif 'P1' == list_for_check_recomb[0] and (list_for_check_recomb != sorted(list_for_check_recomb) ):
                    recombx_reads.append(read)
                if 'P2' == list_for_check_recomb[0] and (list_for_check_recomb == sorted(list_for_check_recomb, reverse=True) ):
                    recomb2_reads.append(read)
                    recomb2_start=list_for_check_recomb_pos[0]
                    recomb2_end=list_for_check_recomb_pos[1]
                    for r in range(1,len(list_for_check_recomb)):
                        if list_for_check_recomb[r] == 'P2':
                            recomb2_start = list_for_check_recomb_pos[r]
                        else:
                            recomb2_end = list_for_check_recomb_pos[r]
                            break
                    cross_region[f"{recomb2_start}-{recomb2_end}"]['recomb2'].append(read)
                    dist = int(recomb2_end) - int(recomb2_start) + 1
                    total_dist = dist + total_dist
                elif 'P2' == list_for_check_recomb[0] and (list_for_check_recomb != sorted(list_for_check_recomb, reverse=True) ):
                    recombx_reads.append(read)

    start_list = [reads_coords[i]['start'] for i in (recomb1_reads + recomb2_reads + recombx_reads + parent1_reads + parent2_reads)]
    end_list = [reads_coords[i]['end'] for i in (recomb1_reads + recomb2_reads + recombx_reads + parent1_reads + parent2_reads)]
    mutaions_list = []
    for read in (recomb1_reads + recomb2_reads + recombx_reads + parent1_reads + parent2_reads):
        tmp = dict()
        for k, v in sorted(mutation_reads[read].items()):
            if 'rev' in v:
                v = [" ".join(v).replace("rev", "rev of ")]

            tmp[k] = v
        mutaions_list.append(tmp)
    
    if argvs.ec19_projdir:
        rel_link = argvs.igv
        igv_list = [f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + '">recombinant 1</a>' for i in recomb1_reads]
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + '">recombinant 2</a>' for i in recomb2_reads])
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + '">recombinant x</a>' for i in recombx_reads])
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + f'">parent {two_parents_list[0]}</a>' for i in parent1_reads])
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + f'">parent {two_parents_list[1]}</a>' for i in parent2_reads])
        df = pd.DataFrame(list(zip((recomb1_reads + recomb2_reads + recombx_reads + parent1_reads + parent2_reads), start_list, end_list, mutaions_list, igv_list)), columns=[
            'read_name', 'start', 'end', 'mutaions_json', 'note'])
    else:
        note_list = [ 'recombinant 1' for i in recomb1_reads ]
        note_list.extend( [ 'recombinant 2' for i in recomb2_reads ] )
        note_list.extend( [ 'recombinant x' for i in recombx_reads ] )
        note_list.extend( [  'parent ' + two_parents_list[0]  for i in parent1_reads] )
        note_list.extend( [  'parent ' + two_parents_list[1]  for i in parent2_reads] )
        df = pd.DataFrame(list(zip((recomb1_reads + recomb2_reads + recombx_reads + parent1_reads + parent2_reads), start_list, end_list, mutaions_list, note_list)), columns=[
            'read_name', 'start', 'end', 'mutaions_json','note'])
    df.sort_values(by=['start'], inplace=True)
    df.to_csv(argvs.tsv, sep="\t",index=False)


    cr_df = pd.DataFrame(((region, json.dumps(dict(sorted(reads.items())))) for region, reads in cross_region.items()), columns=['Cross_region', 'Reads'])
    cr_df['Reads'] = cr_df['Reads']
    cr_df.sort_values(by=['Cross_region'],inplace=True, key = lambda col: [ int(x.split("-")[0]) for x in col] )
    cr_df.to_csv(os.path.splitext(argvs.tsv)[0] + "_by_cross_region.tsv", sep="\t",index=False, quoting=csv.QUOTE_NONE)
    if argvs.ec19_projdir:
        old_width = pd.get_option('display.max_colwidth')
        pd.set_option('display.max_colwidth', None)
        html_file = os.path.splitext(argvs.tsv)[0] + ".html"
        df.to_html(html_file, index=False, escape=False)
        cr_df['Cross_region'] = cr_df['Cross_region'].apply(lambda x: '<a target="_blank" href="{0}?locus=NC_045512_2:{1}">{1}</a>'.format(rel_link,x))
        html_file = os.path.splitext(argvs.tsv)[0] + "_by_cross_region.html"
        cr_df.to_html(html_file, index=False, escape=False)
        pd.set_option('display.max_colwidth', old_width)

    reads_stats['mutation_reads'] = reads_has_two_and_more_mutations_count
    reads_stats['parents'] = ','.join(two_parents_list)
    reads_stats['recomb1_reads'] = len(recomb1_reads)
    reads_stats['recomb2_reads'] = len(recomb2_reads)
    reads_stats['recombx_reads'] = len(recombx_reads)
    reads_stats['parent1_reads'] = len(parent1_reads)
    reads_stats['parent2_reads'] = len(parent2_reads)
    total_recomb_paraents_reads = len(parent1_reads) + len(parent2_reads) + len(recomb1_reads) + len(recomb2_reads) + len(recombx_reads)
    reads_stats['recomb1_perc'] = len(recomb1_reads)/total_recomb_paraents_reads * 100 if total_recomb_paraents_reads > 0 else 0
    reads_stats['recomb2_perc'] = len(recomb2_reads)/total_recomb_paraents_reads * 100 if total_recomb_paraents_reads > 0 else 0
    reads_stats['recombx_perc'] = len(recombx_reads)/total_recomb_paraents_reads * 100 if total_recomb_paraents_reads > 0 else 0
    #reads_stats['mean_dist'] = total_dist/(len(recomb1_reads) + len(recomb2_reads)) if (len(recomb1_reads) + len(recomb2_reads)) > 0 else 0
    logging.info(f"Reads with Variants Mutations (2+) : {reads_stats['mutation_reads']}")
    logging.info(f"Recombinants and Parents Reads count : {total_recomb_paraents_reads}")
    logging.info(f"Recombinant 1 Reads: {reads_stats['recomb1_reads']} ({reads_stats['recomb1_perc']:.2f}%)")
    logging.info(f"Recombinant 2 Reads: {reads_stats['recomb2_reads']} ({reads_stats['recomb2_perc']:.2f}%)")
    logging.info(f"Recombinant x Reads: {reads_stats['recombx_reads']} ({reads_stats['recombx_perc']:.2f}%)")
    logging.info(f"Parent 1 ({two_parents_list[0]}) Reads: {reads_stats['parent1_reads']}")
    logging.info(f"Parent 2 ({two_parents_list[1]}) Reads: {reads_stats['parent2_reads']}")

    return df, reads_stats, parent_pos_count

def write_recombinant_bam(recomb_reads, two_parents_list, argvs):
    samfile = pysam.AlignmentFile(argvs.bam, "rb")
    logging.info("Index input Sam/BAM-file by query name. The index is kept in memory and can be substantial.")
    name_indexed = pysam.IndexedReads(samfile)
    name_indexed.build()
    tmp_bam = os.path.splitext(argvs.outbam)[0] +'.tmp.bam'
    if os.path.exists(tmp_bam):
        os.remove(tmp_bam)
    
    reads_dict=dict()
    reads_dict['recomb1'] = recomb_reads[recomb_reads['note'].str.contains('recombinant 1')]
    reads_dict['recomb2'] = recomb_reads[recomb_reads['note'].str.contains('recombinant 2')]
    reads_dict['recombx'] = recomb_reads[recomb_reads['note'].str.contains('recombinant x')]
    reads_dict['parent1'] = recomb_reads[recomb_reads['note'].str.contains(f'parent {two_parents_list[0]}' )]
    reads_dict['parent2'] = recomb_reads[recomb_reads['note'].str.contains(f'parent {two_parents_list[1]}' )]

    for key in reads_dict:
        outbam_name = os.path.splitext(argvs.outbam)[0] + '.' + key + '.bam'
        logging.info(f"Writting {key} reads to {outbam_name}")
        outbam = pysam.AlignmentFile(tmp_bam, "wb", template=samfile)
        df = reads_dict[key]
        for name in df['read_name']:
            try:
                name_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = name_indexed.find(name)
                for x in iterator:
                    outbam.write(x)
        outbam.close()
        pysam.sort("-o", outbam_name, tmp_bam)
        pysam.index(outbam_name)
        os.remove(tmp_bam)
    samfile.close()

def write_stats(stats, argvs):
    outstats = os.path.splitext(argvs.tsv)[0] + ".stats"
    logging.info(f"Writting recombinant stats to {outstats}")
    with open(outstats,"w") as of:
        of.write("\t".join(stats.keys()) + "\n")
        of.write(str("\t".join(str(x) for x in stats.values())) + "\n")

def mutations_af_plot(parents,nt_to_variants,nt_to_variants_af,nt_to_variants_dp, nt_to_aa_class, reads_stats, recomb_read_count_by_index, parent_read_count_by_index,parent_pos_count,ec19_lineage, ec19_config, nt_to_lineage, argvs):
    output = os.path.splitext(argvs.tsv)[0] + ".mutations_af_plot.html"
    logging.info(f"Generating mutations AF plots and save to {output}")
    all_mut_nt_pos = [ i.split(":")[1] for i in list(nt_to_variants.keys())]
    all_mut_nt = [ i.split(":")[0] + ":<b>" + i.split(":")[1] + "</b>:" + i.split(":")[2] for i in list(nt_to_variants.keys())]
    for i in nt_to_variants:
        if i not in nt_to_lineage and len(nt_to_variants[i])>1:
            nt_to_lineage[i] = nt_to_variants[i]
    igv_url = argvs.igv if argvs.igv else "igv.html" 
    color1='blue'
    color2='red'
    if parents[0] == 'Omicron' or parents[1] == 'Delta':
        color1='blue'
        color2='red'
    if parents[1] == 'Omicron' or parents[0] == 'Delta':
        color1='red'
        color2='blue'
    if recomb_read_count_by_index:
        total_recomb_paraents_reads= reads_stats['parent1_reads'] + reads_stats['parent2_reads'] + reads_stats['recomb1_reads'] + reads_stats['recomb2_reads'] + reads_stats['recombx_reads']
        recomb_rate_by_index = { i+0.5 : recomb_read_count_by_index[i+0.5]/(recomb_read_count_by_index[i+0.5] + parent_read_count_by_index[i+0.5]) if (i + 0.5) in recomb_read_count_by_index else 0 for i in range(len(all_mut_nt_pos)-1)} 
        recomb_count_by_index = [ recomb_read_count_by_index[i+0.5] if (i + 0.5) in recomb_read_count_by_index else 0 for i in range(len(all_mut_nt_pos)-1)]
        recombinant_rate = list(recomb_rate_by_index.values())
        fig = make_subplots(rows=2,cols=2, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.01 ,specs=[[{"type": "xy"}, {"type": "pie"}],[{"type": "xy"}, {"type": "pie"}]], row_heights=[200,600],column_widths=[0.9,0.2])
        #fig = make_subplots(rows=2,cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[200,600])
        main_row=2
        ## recombinate rate track
        fig.update_yaxes(range=[-0.01, max(recombinant_rate)+0.01],title='Recombinate Rate', row=1, col=1 )
        fig.add_trace(go.Bar(
            x=list(recomb_rate_by_index.keys()), 
            y=[ i if i > 0 else None for i in recombinant_rate], 
            #mode='lines',
            #line=dict(color='black'),
            marker_color= 'black',marker_line_color='black',
            name='Recombinant Rate',
            width=[0.3 if i > 0.01 else 1 for i in recombinant_rate],
            hovertemplate = 'Rate: %{y:.6f}<br>'+'%{hovertext}<extra></extra>',
            hovertext=[ 'Count: ' + str(recomb_count_by_index[x]) + '<br>CR: '+ all_mut_nt[x] + "-" + all_mut_nt[x+1] for x in range(len(all_mut_nt)-1)],
            #customdata=[ url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if parents[0] in nt_to_variants[x] and parents[1] not in nt_to_variants[x] else None for x in nt_to_variants_af ],
            showlegend=False,
            ),row=1,col=1)
        for i in range(len(all_mut_nt_pos)-1):
            if ((i + 0.5) in recomb_rate_by_index and recomb_rate_by_index[i+0.5] > 0) or (int(all_mut_nt_pos[i]) in parent_pos_count and int(all_mut_nt_pos[i+1]) in parent_pos_count) :
                fig.add_vrect(
                    x0=i, x1=i+1,
                    layer="below", line_width=0, fillcolor="grey", opacity=0.2,
                    row=1,col=1
                )
    else:
        fig = make_subplots(rows=1,cols=2,specs=[[{"type": "xy"}, {"type": "pie"}]],column_widths=[0.9,0.2])
        main_row=1
    ## mutation scatter plot
    fig.add_trace(go.Scatter(
                x=list(range(len(all_mut_nt))), 
                y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants], 
                mode='markers',
                marker=dict(color=color1,size=10),
                name=parents[0],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants],
                customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants ],
                showlegend=True,
                ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
                x=list(range(len(all_mut_nt))), 
                y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants], 
                mode='markers',
                marker=dict(color=color2,size=10),
                name=parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants],
                customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants ],
                showlegend=True,
                ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
                x=list(range(len(all_mut_nt))), 
                y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants], 
                mode='markers',
                marker=dict(color='purple',size=10),
                name=parents[0] + ', ' + parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants],
                customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants ],
                showlegend=True,
                ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
                x=list(range(len(all_mut_nt))), 
                y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants], 
                mode='markers',
                marker=dict(color='green',size=10),
                name='Not ' + parents[0] + ', Not ' + parents[1],
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ str(nt_to_variants_dp[x]) + '<br>' + 'AA: ' + nt_to_aa_class.convert_nt_prot(x) + '<br>Var: ' + ','.join(nt_to_lineage[x]) if x in nt_to_lineage and nt_to_lineage[x] and  parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants],
                customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants ],
                showlegend=True,
                ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
                x=list(range(len(all_mut_nt))), 
                y=[ float(nt_to_variants_af[x])  if x not in nt_to_lineage else None for x in nt_to_variants], 
                mode='markers',
                marker=dict(color='grey',size=10),
                name='Undefined Mutations',
                hovertemplate = 'Mut: %{x}<br>' + 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
                text=[ 'DP: '+ str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x not in nt_to_lineage else None for x in nt_to_variants],
                customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x not in nt_to_lineage else None for x in nt_to_variants ],
                showlegend=True,
                ),row=main_row,col=1)
    # pie chart
    pie_content={parents[0]:0, parents[1]:0, parents[0] + ', ' + parents[1]:0, 'Not ' + parents[0] + ', Not ' + parents[1]:0,'Undefined_Mutations':0}
    for x in nt_to_variants:
        if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x]:
            pie_content[parents[0]] += 1
        if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x]:
            pie_content[parents[1]] += 1
        if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x]:
            pie_content[parents[0] + ', ' + parents[1]] += 1
        if x in nt_to_lineage and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x]:
            pie_content['Not ' + parents[0] + ', Not ' + parents[1]] += 1
        if x not in nt_to_lineage:
            pie_content['Undefined_Mutations'] += 1
    fig.add_trace(go.Pie(labels=list(pie_content.keys()),
                        values=list(pie_content.values()),
                        textinfo='value+percent',
                        marker=dict(colors=[color1,color2,'purple','green','grey'] ),
                        showlegend=False,
                        hovertemplate = 'count: %{value}<br>' + 'percent: %{percent:.1%}<br>'+'<extra></extra>',
                        ), row=main_row, col=2
        )
    fig.update_yaxes(range=[0, 1.1],title='Alternative Frequency', row=main_row, col=1)
    fig.update_xaxes(tickvals = list(range(len(all_mut_nt))), ticktext=all_mut_nt, tickfont=dict(size=10),tickangle=-60)
    if 'projname' in ec19_config:
        fig.update_layout(title=ec19_config['projname'] + ' (' + ec19_lineage + ')')
    fig.update_layout(legend=dict(
        yanchor="top",
        y=1.05,
        xanchor="right",
        x=0.95,
        bgcolor="rgba(255,255,255,0.1)",
    ))
    # Get HTML representation of plotly.js and this figure
    plot_div = plot(fig, output_type='div', include_plotlyjs=True)

    # Get id of html div element that looks like
    # <div id="301d22ab-bfba-4621-8f5d-dc4fd855bb33" ... >
    res = re.search('<div id="([^"]*)"', plot_div)
    div_id = res.groups()[0]

    # Build JavaScript callback for handling clicks
    # and opening the URL in the trace's customdata 
    js_callback = """
    <script>
    var plot_element = document.getElementById("{div_id}");
    plot_element.on('plotly_click', function(data){{
        var point = data.points[0];
        if (point) {{
            window.open(point.customdata);
        }}
    }})
    </script>
    """.format(div_id=div_id)

    # Build HTML string
    html_str = """
    <html>
    <body>
    {plot_div}
    {js_callback}
    </body>
    </html>
    """.format(plot_div=plot_div, js_callback=js_callback)
    with open(output, 'w') as f:
        f.write(html_str)
    fig.write_image(output+'.png', width=1980, height=1080)
    return output

def mutations_af_plot_genome(parents,nt_to_variants,nt_to_variants_af,nt_to_variants_dp, nt_to_aa_class, reads_stats, recomb_read_count_by_pos, parent_pos_count, ec19_lineage, ec19_config, nt_to_lineage, argvs):
    output = os.path.splitext(argvs.tsv)[0] + ".mutations_af_plot_genomeview.html"
    logging.info(f"Generating mutations AF plots Genome View and save to {output}")
    all_mut_nt_pos = [ i.split(":")[1] for i in list(nt_to_variants.keys())]
    all_mut_nt = [ i.split(":")[0] + ":<b>" + i.split(":")[1] + "</b>:" + i.split(":")[2] for i in list(nt_to_variants.keys())]
    for i in nt_to_variants:
        if i not in nt_to_lineage and len(nt_to_variants[i])>1:
            nt_to_lineage[i] = nt_to_variants[i]
    igv_url = argvs.igv if argvs.igv else "igv.html" 
    color1='blue'
    color2='red'
    if parents[0] == 'Omicron' or parents[1] == 'Delta':
        color1='blue'
        color2='red'
    if parents[1] == 'Omicron' or parents[0] == 'Delta':
        color1='red'
        color2='blue'
    if recomb_read_count_by_pos: 
        total_recomb_paraents_reads= reads_stats['parent1_reads'] + reads_stats['parent2_reads'] + reads_stats['recomb1_reads'] + reads_stats['recomb2_reads'] + reads_stats['recombx_reads'] 
        for i in recomb_read_count_by_pos:
            if i not in parent_pos_count:
                parent_pos_count[i]=0
        recomb_rate_by_pos = { i:recomb_read_count_by_pos[i]/ (recomb_read_count_by_pos[i] + parent_pos_count[i]) for i in recomb_read_count_by_pos}
        recombinant_rate = list(recomb_rate_by_pos.values())
        #fig = make_subplots(rows=2,cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[200,600])
        fig = make_subplots(rows=2,cols=2, shared_xaxes=True, vertical_spacing=0.02, horizontal_spacing=0.01 ,specs=[[{"type": "xy"}, {"type": "pie"}],[{"type": "xy"}, {"type": "pie"}]], row_heights=[200,600],column_widths=[0.9,0.2])

        main_row=2
        ## recombinate rate track
        fig.update_yaxes(range=[0,max(recombinant_rate)+0.01],title='Recombinate Rate', row=1,col=1)
        fig.add_trace(go.Bar(
            x=[ int(x) for x in recomb_rate_by_pos], 
            y=recombinant_rate, 
            marker_color= 'black',marker_line_color='black',
            name='Recombinant Rate',
            width=1,
            hovertemplate = 'Pos: %{x}<br>'+'Rate: %{y:.6f}<br>'+'Count: %{hovertext}<extra></extra>',
            hovertext=list(recomb_read_count_by_pos.values()),
            showlegend=False), 
            row=1,col=1)
        
        recomb_read_count_pos_list = list(set(list(recomb_rate_by_pos.keys()) + list(parent_pos_count.keys())))
        shade_start = recomb_read_count_pos_list[0]
        for i in range(1,len(recomb_read_count_pos_list)-1):
        
            if recomb_read_count_pos_list[i] + 1 != recomb_read_count_pos_list[i+1]:
                fig.add_vrect(
                    x0=shade_start, x1=recomb_read_count_pos_list[i],
                    layer="below", line_width=0, fillcolor="grey", opacity=0.2,
                    row=1,col=1
                )
                shade_start = recomb_read_count_pos_list[i+1]

        fig.add_vrect(
                x0=shade_start, x1=recomb_read_count_pos_list[-1],
                layer="below", line_width=0, fillcolor="grey", opacity=0.2,
                row=1,col=1
            )
    else:
        fig = make_subplots(rows=1,cols=2, horizontal_spacing=0.01 ,specs=[[{"type": "xy"}, {"type": "pie"}]],column_widths=[0.9,0.2])
        main_row=1

    ## mutation scatter plot
    fig.add_trace(go.Scatter(
            x=[ int(x) for x in all_mut_nt_pos], 
            y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants], 
            mode='markers',
            marker=dict(color=color1,size=10),
            #marker=dict(color=color1,size=10,line_color='pink',line_width=0.8, symbol=[ 'circle-x' if parents[0] in nt_to_variants[x] and parents[1] not in nt_to_variants[x] and x.split(':')[1] in cr_coord else 'circle' for x in nt_to_variants]), 
            name=parents[0],
            hovertemplate = 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
            text=[ 'Mut: '  + x + '<br>DP: ' + str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants],
            customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants ],
            showlegend=True,
            ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
            x=[ int(x) for x in all_mut_nt_pos], 
            y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants], 
            mode='markers',
            marker=dict(color=color2,size=10),
            #marker=dict(color=color2,size=10,line_color='Skyblue',line_width=0.8, symbol=[ 'circle-x' if  parents[1] in nt_to_variants[x] and parents[0] not in nt_to_variants[x] and x.split(':')[1] in cr_coord else 'circle' for x in nt_to_variants]), 
            name=parents[1],
            hovertemplate = 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
            text=[ 'Mut: '  + x + '<br>DP: ' + str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants],
            customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] else None for x in nt_to_variants ],
            showlegend=True,
            ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
            x=[ int(x) for x in all_mut_nt_pos], 
            y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants], 
            mode='markers',
            marker=dict(color='purple',size=10),
            #marker=dict(color='purple',size=10,line_color='Yellow',line_width=0.8, symbol=[ 'circle-x' if parents[0] in nt_to_variants[x] and parents[1] in nt_to_variants[x] and x.split(':')[1] in cr_coord else 'circle' for x in nt_to_variants]), 
            name=parents[0] + ', ' + parents[1],
            hovertemplate = 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
            text=[ 'Mut: '  + x + '<br>DP: ' + str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants],
            customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x] else None for x in nt_to_variants ],
            showlegend=True,
            ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
            x=[ int(x) for x in all_mut_nt_pos], 
            y=[ float(nt_to_variants_af[x])  if x in nt_to_lineage and nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants], 
            mode='markers',
            marker=dict(color='green',size=10),
            name='Not ' + parents[0] + ', Not ' + parents[1],
            hovertemplate = 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
            text=[ 'Mut: '  + x + '<br>DP: ' + str(nt_to_variants_dp[x]) + '<br>' + 'AA: ' + nt_to_aa_class.convert_nt_prot(x) + '<br>Var: ' + ','.join(nt_to_lineage[x]) if x in nt_to_lineage and nt_to_lineage[x] and  parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants],
            customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x in nt_to_lineage and nt_to_lineage[x] and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x] else None for x in nt_to_variants ],
            showlegend=True,
            ),row=main_row,col=1)
    fig.add_trace(go.Scatter(
            x=[ int(x) for x in all_mut_nt_pos], 
            y=[ float(nt_to_variants_af[x])  if x not in nt_to_lineage else None for x in nt_to_variants],  
            mode='markers',
            marker=dict(color='grey',size=10),
            name='Undefined Mutations',
            hovertemplate = 'AF: %{y:.2f}<br>'+'%{text}<extra></extra>',
            text=[ 'Mut: '  + x + '<br>DP: ' + str(nt_to_variants_dp[x]) + '<br>AA: ' + nt_to_aa_class.convert_nt_prot(x) if x not in nt_to_lineage else None for x in nt_to_variants],
            customdata=[ igv_url + '?locus=NC_045512_2:' + str(int(x.split(':')[1]) - 100) + '-' +  str(int(x.split(':')[1]) + 100) if x not in nt_to_lineage else None for x in nt_to_variants ],
            showlegend=True,
            ),row=main_row,col=1)
    ### Genome annotation
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=266,
        y=1.02,
        ax=13468,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=13468,
        y=1.02,
        ax=21555,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=21563,
        y=1.02,
        ax=25384,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=25393,
        y=1.02,
        ax=26220,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=26245,
        y=1.02,
        ax=26472,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=26523,
        y=1.02,
        ax=27191,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=27202,
        y=1.02,
        ax=27387,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=27394,
        y=1.02,
        ax=27759,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=27756,
        y=1.02,
        ax=27887,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=27894,
        y=1.02,
        ax=28259,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=28274,
        y=1.02,
        ax=29533,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        axref='x', 
        x=29558,
        y=1.02,
        ax=29674,
        ay=1.02,
        showarrow=True,
        arrowhead=2,
        arrowside='start',
        arrowsize=1,
        arrowwidth=1.5
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=6000,
        y=1.06,
        text="ORF1a",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=17500,
        y=1.06,
        text="ORF1b",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=23500,
        y=1.06,
        text="Spike",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=25800,
        y=1.06,
        text="3a",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=26350,
        y=1.06,
        text="E",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=26800,
        y=1.06,
        text="M",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=27300,
        y=1.06,
        text="6",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=27700,
        y=1.06,
        text="7ab",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=28100,
        y=1.06,
        text="8",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=29000,
        y=1.06,
        text="N",
        showarrow=False,
        font_size=12,
    )
    fig.add_annotation(
        xref='x', yref='paper',
        x=29600,
        y=1.06,
        text="10",
        showarrow=False,
        font_size=12,
    )
    fig.add_vline(x=13468, line_width=1, line_color="white",col=1)
    fig.add_vline(x=21556, line_width=1, line_color="white",col=1)
    fig.add_vline(x=25385, line_width=1, line_color="white",col=1)
    fig.add_vline(x=26221, line_width=1, line_color="white",col=1)
    fig.add_vline(x=26473, line_width=1, line_color="white",col=1)
    fig.add_vline(x=27192, line_width=1, line_color="white",col=1)
    fig.add_vline(x=27388, line_width=1, line_color="white",col=1)
    fig.add_vline(x=27888, line_width=1, line_color="white",col=1)
    fig.add_vline(x=28260, line_width=1, line_color="white",col=1)
    fig.add_vline(x=29554, line_width=1, line_color="white",col=1)
    fig.add_vline(x=29675, line_width=1, line_color="white",col=1)
    # pie chart
    pie_content={parents[0]:0, parents[1]:0, parents[0] + ', ' + parents[1]:0, 'Not ' + parents[0] + ', Not ' + parents[1]:0,'Undefined_Mutations':0}
    for x in nt_to_variants:
        if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x]:
            pie_content[parents[0]] += 1
        if x in nt_to_lineage and parents[1] in nt_to_lineage[x] and parents[0] not in nt_to_lineage[x]:
            pie_content[parents[1]] += 1
        if x in nt_to_lineage and parents[0] in nt_to_lineage[x] and parents[1] in nt_to_lineage[x]:
            pie_content[parents[0] + ', ' + parents[1]] += 1
        if x in nt_to_lineage and parents[0] not in nt_to_lineage[x] and parents[1] not in nt_to_lineage[x]:
            pie_content['Not ' + parents[0] + ', Not ' + parents[1]] += 1
        if x not in nt_to_lineage:
            pie_content['Undefined_Mutations'] += 1
    fig.add_trace(go.Pie(labels=list(pie_content.keys()),
                        values=list(pie_content.values()),
                        textinfo='value+percent',
                        marker=dict(colors=[color1,color2,'purple','green','grey'] ),
                        showlegend=False,
                        hovertemplate = 'count: %{value}<br>' + 'percent: %{percent:.1%}<br>'+'<extra></extra>',
                        ), row=main_row, col=2
        )
    # Set axes properties
    fig.update_xaxes(range=[1, 29903], showgrid=False)
    fig.update_yaxes(range=[0, 1.1],title='Alternative Frequency', row=main_row, col=1)
    if 'projname' in ec19_config:
        fig.update_layout(title=ec19_config['projname'] + ' (' + ec19_lineage + ')')
    fig.update_layout(legend=dict(
        yanchor="top",
        y=1.05,
        xanchor="right",
        x=0.95,
        bgcolor="rgba(255,255,255,0.1)",
    ))
    # Get HTML representation of plotly.js and this figure
    plot_div = plot(fig, output_type='div', include_plotlyjs=True)

    # Get id of html div element that looks like
    # <div id="301d22ab-bfba-4621-8f5d-dc4fd855bb33" ... >
    res = re.search('<div id="([^"]*)"', plot_div)
    div_id = res.groups()[0]

    # Build JavaScript callback for handling clicks
    # and opening the URL in the trace's customdata 
    js_callback = """
    <script>
    var plot_element = document.getElementById("{div_id}");
    plot_element.on('plotly_click', function(data){{
        var point = data.points[0];
        if (point) {{
            window.open(point.customdata);
        }}
    }})
    </script>
    """.format(div_id=div_id)

    # Build HTML string
    html_str = """
    <html>
    <body>
    {plot_div}
    {js_callback}
    </body>
    </html>
    """.format(plot_div=plot_div, js_callback=js_callback)
    with open(output, 'w') as f:
        f.write(html_str)
    fig.write_image(output+'.png', width=1980, height=1080)
    return output

def update_igv_html(two_parents_list,reads_stats, argvs):
    igv_html_file = os.path.join(os.path.dirname(argvs.tsv), argvs.igv)
    if not os.path.exists(igv_html_file):
        return
    update_igv_html_file = os.path.splitext(igv_html_file)[0] + '.recombreads.html'
    variants_gff = os.path.join(os.path.dirname(igv_html_file), 'variant_mutation.gff')
    if argvs.igv_variants:
        cmd=[ os.path.join(bin_dir,'varJSonToGFF.py'), '--refacc', argvs.refacc, '--variantMutation', argvs.variantMutation, '--out', variants_gff ]
        process_cmd(cmd)
    logging.info(f"Adding recombinant reads as track in IGV view to {igv_html_file}")
    variants_track = {'name':"Variants Mutations", 'format':"gff3", 'displayMode':"expanded", 'height': 150, 'url': "./variant_mutation.gff", 'indexed': False, 'visibilityWindow':32000, 'colorBy':"Variant"}
    recomb1_track = {
        'name': 'Recombinant 1', 
        'type':'alignment', 
        'format': 'bam', 
        'colorBy': 'strand', 
        'url': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recomb1.bam', 
        'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recomb1.bam.bai',
        'squishedRowHeight': 10,
        'height': 250,
        'displayMode': 'SQUISHED' }
    recomb2_track = {
        'name': 'Recombinant 2', 
        'type':'alignment', 
        'format': 'bam', 
        'colorBy': 'strand', 
        'url': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recomb2.bam', 
        'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recomb2.bam.bai',
        'squishedRowHeight': 10,
        'height': 250,
        'displayMode': 'SQUISHED' }
    recombx_track = {
        'name': 'Recombinant x', 
        'type':'alignment', 
        'format': 'bam', 
        'colorBy': 'strand', 
        'url': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recombx.bam', 
        'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.recombx.bam.bai',
        'squishedRowHeight': 10,
        'height': 250,
        'displayMode': 'SQUISHED' }
    parent1_track = {
        'name': f'Parent {two_parents_list[0]}', 
        'type':'alignment', 
        'format': 'bam', 
        'colorBy': 'strand', 
        'url': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.parent1.bam', 
        'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.parent1.bam.bai',
        'squishedRowHeight': 10,
        'height': 250,
        'displayMode': 'SQUISHED' }
    parent2_track = {
        'name': f'Parent {two_parents_list[1]}', 
        'type':'alignment', 
        'format': 'bam', 
        'colorBy': 'strand', 
        'url': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.parent2.bam', 
        'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/MixedInfectionAnalysis/recombinant_reads.parent2.bam.bai',
        'squishedRowHeight': 10,
        'height': 250,
        'displayMode': 'SQUISHED' }
    of = open(update_igv_html_file,'w')
    with open(igv_html_file, 'r') as f:
        for line in f:
            if 'options =' in line:
                if 'options = {"tr' in line or 'options = {"re' in line or 'options = {"locus' in line:
                    options = [i.strip().replace(";","") for i in line.split("=")]
                    options_dict = json.loads(options[1])
                    createBrowser=""
                elif 'options = {' in line:
                    nextline = f.readline()
                    options = "{"
                    while 'igv.createBrowser' not in nextline:
                        nextline = re.sub(r'(\w+):', r'"\1":', nextline)
                        options = options + nextline
                        nextline = f.readline()
                        if 'body' in nextline or 'html' in nextline:
                            break
                    createBrowser=nextline
                    options=options.replace(";","").replace('""NC_045512_2"','"NC_045512_2').replace("'",'"')
                    options_dict = json.loads(options)
                    
                new_tracks=[]
                for i in options_dict['tracks']:
                    if 'name' in i and ('Recombinant' in i['name'] or 'Parent' in i['name'] or 'Variants Mutations' in i['name']):
                        pass
                    elif 'name' in i and (i['name'] == "NC_045512.2 Alignment" or i['name'] == "EC-19 Alignment" or i['name'] == "Alignment" ):
                        new_tracks.append(i)
                        if argvs.igv_variants:
                            new_tracks.append(variants_track)
                        if reads_stats['recomb1_reads'] > 0:
                            new_tracks.append(recomb1_track)
                        if reads_stats['recomb2_reads'] > 0:
                            new_tracks.append(recomb2_track)
                        if reads_stats['recombx_reads'] > 0:
                            new_tracks.append(recombx_track)
                        if reads_stats['parent1_reads'] > 0:
                            new_tracks.append(parent1_track)
                        if reads_stats['parent2_reads'] > 0:
                            new_tracks.append(parent2_track)
                    else:
                        if 'showGenotypes' not in i:
                            new_tracks.append(i)
                options_dict['tracks'] = new_tracks
                of.write("var options = " + json.dumps(options_dict) + ";\n" + createBrowser + "\n")
            else:
                of.write(line)
    of.close()
    shutil.move(update_igv_html_file,igv_html_file)

def calculate_recombinant_rate(cr_file,reads_stats,filtered_nt_to_variants,parent_pos_count, argvs):
    recomb_rate=list()
    df_cr= pd.read_table(cr_file,sep="\t")
    df_cr = df_cr.sort_values(by="Cross_region",key = lambda col: [ int(x.split("-")[0]) for x in col] )
    filtered_nt_to_variants_list = list(filtered_nt_to_variants.keys())
    mut_nt_pos = [ i.split(":")[1] for i in filtered_nt_to_variants_list ]
    recomb_read_count_by_index=dict()
    recomb_read_count_by_pos=dict()
    recomb_rate_by_pos=dict()
    recomb_rate_by_index=dict()
    parent_read_count_by_index=dict()
    for index, row in df_cr.iterrows():
        region=row['Cross_region']
        reads= json.loads(row['Reads'])
        recomb_read_count = 0
        for type, read in reads.items():
            recomb_read_count = recomb_read_count + len(read)
        index1 = mut_nt_pos.index(region.split('-')[0])
        index2 = mut_nt_pos.index(region.split('-')[1])
        parent_read_count_in_cr=0
        region_size =0
        for i in range(int(region.split('-')[0])+1,int(region.split('-')[1])+1):
            recomb_read_count_by_pos[i] = recomb_read_count_by_pos[i] + recomb_read_count if i in recomb_read_count_by_pos else recomb_read_count
            parent_read_count_in_cr = parent_pos_count[i] + parent_read_count_in_cr if i in parent_pos_count else parent_read_count_in_cr + 0
            region_size += 1
        avg_parent_read_count_in_cr = int(parent_read_count_in_cr/region_size)
        for i in range(index2 - index1):
            recomb_read_count_by_index[index1 + i + 0.5] = recomb_read_count_by_index[index1 + i + 0.5] + recomb_read_count if (index1 + i + 0.5) in recomb_read_count_by_index else recomb_read_count
            parent_read_count_by_index[index1 + i + 0.5] = parent_read_count_by_index[index1 + i + 0.5] + avg_parent_read_count_in_cr if (index1 + i + 0.5) in parent_read_count_by_index else avg_parent_read_count_in_cr
    total_recomb_paraents_reads= reads_stats['parent1_reads'] + reads_stats['parent2_reads'] + reads_stats['recomb1_reads'] + reads_stats['recomb2_reads'] + reads_stats['recombx_reads'] 
    recomb_rate_by_index = { i+0.5 : recomb_read_count_by_index[i+0.5]/(parent_read_count_by_index[i+0.5]+recomb_read_count_by_index[i+0.5]) if (i + 0.5) in recomb_read_count_by_index and (parent_read_count_by_index[i+0.5]+recomb_read_count_by_index[i+0.5]) > 0 else 0 for i in range(len(mut_nt_pos)-1)}
    recomb_rate_by_pos = { i:recomb_read_count_by_pos[i]/total_recomb_paraents_reads for i in recomb_read_count_by_pos}
    #cr_coord = list({ i  for x in df_cr.Cross_region for i in x.split('-') })
    tmp_output = os.path.splitext(argvs.tsv)[0] + ".recombinant_rate_pos.tsv"
    with open(tmp_output, "w") as f:
        for key, value in recomb_read_count_by_pos.items():
            rate = value/(parent_pos_count[key] + value) if key in parent_pos_count else 1.0
            parent_depth = parent_pos_count[key] if key in parent_pos_count else 0
            f.write(f"{key}\t{value}\t{parent_depth}\t{rate}\n")

    for i in range(len(mut_nt_pos)-1):
        if (i + 0.5) in recomb_read_count_by_index:
            #logging.info(filtered_nt_to_variants_list[i] + "," + str(recomb_read_count_by_index[i + 0.5]) + "," + str(recomb_rate_by_index[i + 0.5]))
            if (recomb_read_count_by_index[i + 0.5] >= argvs.minRecombReadCount) and (recomb_rate_by_index[i + 0.5] >= argvs.minRecombRate):
                logging.info(f"There is hotspot of break point at region between {filtered_nt_to_variants_list[i]} and {filtered_nt_to_variants_list[i+1]} with {recomb_read_count_by_index[i + 0.5]} recombinant reads and {recomb_rate_by_index[i + 0.5]:.2%} recombinant rate.")
    return recomb_read_count_by_index, recomb_read_count_by_pos, parent_read_count_by_index

def parse_lineage(txt):
    with open(txt,'r') as f:
        content = f.readlines()
    result = content[1].split(',')
    return(result[1])

def parse_ec19_config(config):
    ec19_config=dict()
    with open(config,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.__contains__('='):
                k,v =line.strip().split('=')
                ec19_config[k]=v
    return(ec19_config)

def main():
    argvs = setup_argparse()
    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
    format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')

    (delta_uniq_nt, omicron_uniq_nt, nt_to_variant,
     nt_to_aa) = load_var_mutation(argvs.variantMutation)
    nt_to_lineage = load_lineage_mutation(argvs.lineageMutation)
    
    ec19_lineage =  'Unknown'
    ec19_config = dict()
    if argvs.ec19_projdir:
        ec19_config_file = os.path.join(argvs.ec19_projdir, "config.txt")
        ec19_lineage_file = os.path.join(argvs.ec19_projdir, "ReadsBasedAnalysis","readsMappingToRef","NC_045512.2_consensus_lineage.txt")
        ec19_lineage = parse_lineage(ec19_lineage_file) if ec19_lineage_file else 'Unknown'
        ec19_config =  parse_ec19_config(ec19_config_file)

    if argvs.vcf:
        parents_variant, filtered_nt_to_variants, filtered_nt_to_variants_af, filtered_nt_to_variants_dp=parse_vcf(argvs,nt_to_variant)
        mutation_reads, reads_coords, reads_stats = find_read_with_variants(
        filtered_nt_to_variants, argvs)
        two_parents_list = list(parents_variant.keys())[0:2]
    else:
        mutation_reads, reads_coords, reads_stats = find_read_with_variants(
        nt_to_variant, argvs)
        parents_variant = find_parents_by_mutation_reads(mutation_reads,argvs)
        two_parents_list = list(parents_variant.keys())[0:2]
        #two_parents_list = argvs.recombinant_variants

    if (len(two_parents_list) == 1):
        two_parents_list.append("Wuhan-Hu-1")
    if (len(two_parents_list) == 0):
        logging.error("No parents detected.")
        sys.exit(1)
    # if list(parents_v.keys())[0]['uniq'] > 2 and list(parents_v.keys())[1]['uniq'] > 2: #both parents need at least have two muations. recombinant?
    
    recomb_reads_df,  reads_stats, parent_pos_count = find_recomb(mutation_reads, reads_coords, two_parents_list, reads_stats, argvs)
    
    write_stats(reads_stats,argvs)

    write_recombinant_bam(recomb_reads_df, two_parents_list, argvs)

    if argvs.vcf and argvs.mutations_af_plot and  (len(two_parents_list) == 2):
        try:
            nt_to_aa_class = translate.NT_AA_POSITION_INTERCHANGE(importlib_resources.files(toolname).joinpath('data/SARS-CoV-2.json'))
        except:
            nt_to_aa_class = translate.NT_AA_POSITION_INTERCHANGE(os.path.join(bin_dir,"data", 'SARS-CoV-2.json'))
        cr_file  = os.path.splitext(argvs.tsv)[0] + "_by_cross_region.tsv"
        recomb_rate_by_index=dict()
        recomb_rate_by_pos=dict()
        if os.path.exists(cr_file):
            recomb_rate_by_index,  recomb_rate_by_pos, parent_read_count_by_index = calculate_recombinant_rate(cr_file,reads_stats,filtered_nt_to_variants,parent_pos_count, argvs)
        mutations_af_plot(two_parents_list, filtered_nt_to_variants, filtered_nt_to_variants_af, filtered_nt_to_variants_dp, nt_to_aa_class, reads_stats, recomb_rate_by_index, parent_read_count_by_index, parent_pos_count,ec19_lineage, ec19_config, nt_to_lineage, argvs)
        mutations_af_plot_genome(two_parents_list, filtered_nt_to_variants, filtered_nt_to_variants_af, filtered_nt_to_variants_dp, nt_to_aa_class, reads_stats, recomb_rate_by_pos, parent_pos_count, ec19_lineage, ec19_config, nt_to_lineage, argvs)
    if argvs.igv:
        update_igv_html(two_parents_list,reads_stats,argvs)
    # else:
        #logging.error(f"both parents need at least have two muations {parents_variant}.")
    

if __name__ == '__main__':
    main()
