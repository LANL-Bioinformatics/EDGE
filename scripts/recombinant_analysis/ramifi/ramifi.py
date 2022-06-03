#!/usr/bin/env python3
import argparse as ap
import json
import os
import sys
import re
import shutil
from collections import defaultdict

import pandas as pd
import pysam
# standardize the logging output
import logging

toolname = os.path.basename(__file__)
bin_dir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()
sys.path.append(bin_dir)
from __init__ import __version__

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
    parser.add_argument('--minMixed_n',metavar='[INT]',required=False, type=int, default=3, help="threadhold of mixed mutations count for vcf.")
    parser.add_argument('--minReadCount',metavar='[INT]',required=False, type=int, default=10, help="threadhold of read with variant count when no vcf provided.")
    parser.add_argument(
        '--lineageMutation', metavar='[FILE]', required=False, type=str, help=f"lineage mutation json file [default: variant_mutation.json]")
    parser.add_argument(
        '--variantMutation', metavar='[FILE]', required=False, type=str, help=f"variant mutation json file [default: lineage_mutation.json]")
    
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
                                     'readsMappingToRef', tsv_filename)
        else:
            argvs.tsv = tsv_filename
    if argvs.ec19_projdir and not argvs.igv:
        argvs.igv = os.path.join('..', '..', 'IGV', 'ref_tracks', 'igv.html')
    if not argvs.outbam:
        argvs.outbam = os.path.splitext(argvs.tsv)[0] + ".bam"
    return argvs


def load_var_mutation(file):
    with open(file, 'r') as f:
        v_data = json.load(f)
    delta_uniq_nt = list(
        set(v_data['Delta'].keys()) - set(v_data['Omicron'].keys()))
    omicron_uniq_nt = list(
        set(v_data['Omicron'].keys()) - set(v_data['Delta'].keys()))

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
    parents_v_uniq=dict()
    parents_v_all=dict()
    parents_v_in_af_range=dict()
    mix_count=0
    mix_count_in_known_variants=dict()
    mutations_count=0
    with open(argvs.vcf,'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            content = line.strip().split('\t')
            AFreq = 0
            if "ALT_FREQ" in line:
                content2 = content[-1].split(':')
                AFreq = float(content2[-1])
            if 'DP=' in line and 'DP4=' in line:
                m = re.search(r'DP=(\d+)', line)
                depth=m[1]
                m = re.search(r'DP4=(\d+),(\d+),(\d+),(\d+)', line)
                (ref_foward,ref_reverse,alt_forward,alt_reverse) = m.groups()
                AFreq = float((int(alt_forward) + int(alt_reverse)/int(depth)))
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
                        if ''.join(nt_to_variants[nt_v]) in parents_v_uniq:
                            parents_v_uniq[''.join(nt_to_variants[nt_v])] += 1 
                        else:
                            parents_v_uniq[''.join(nt_to_variants[nt_v])] = 1
                            parents_v_all[''.join(nt_to_variants[nt_v])] = 0
                            parents_v_in_af_range[''.join(nt_to_variants[nt_v])] = 0
                    filtered_nt_to_variants[nt_v]=nt_to_variants[nt_v]
                    filtered_nt_to_variants_af[nt_v]=AFreq
    # use the first scan, uniq list to count all mutataions with the variants and variants between AF range
    for v in parents_v_uniq:
        for nt_v in filtered_nt_to_variants:
            if v in filtered_nt_to_variants[nt_v]:
                parents_v_all[v] += 1
                if filtered_nt_to_variants_af[nt_v] > argvs.minMixAF and filtered_nt_to_variants_af[nt_v] < argvs.maxMixAF:
                    parents_v_in_af_range[v] += 1

    logging.debug(f"{mix_count}/{mutations_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF}.")
    logging.debug(f"{len(mix_count_in_known_variants)}/{mix_count} mutations (positions) have allelic frequency between {argvs.minMixAF} and {argvs.maxMixAF} and known variants.")
    ## sort by observerd variants counts
    parents_v_all = dict(sorted(parents_v_all.items(), key=lambda item: item[1], reverse=True))
    logging.debug(f"All probable parents, mutation count: {parents_v_all}")
    logging.debug(f"All probable parents, mutation count with AF {argvs.minMixAF}-{argvs.maxMixAF}: {parents_v_in_af_range}")
    if len(parents_v_all.keys())<2:
        logging.error(f'no parents variants detected. {parents_v_all}')
        sys.exit(1)
    if mix_count <= argvs.minMixed_n:
        logging.error(f'count of mixed mutations with AF between {argvs.minMixAF} and {argvs.maxMixAF} is less than {argvs.minMixed_n}.')
        sys.exit(1)
    return parents_v_all, filtered_nt_to_variants

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
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.
                        if ref == 'del' or ref == 'ins':
                            continue
                        if alt == pileupread.alignment.query_sequence[pileupread.query_position]:
                            for mut in nt_to_variant[i]:
                                if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                    mutation_reads[pileupread.alignment.query_name][pos].extend([mut])
                        elif ref == pileupread.alignment.query_sequence[pileupread.query_position]:
                            if 'ref' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(
                                    ['ref'] + nt_to_variant[i])
                        else:
                            other_mut = ref + ":" + \
                                str(pos) + ":" + \
                                pileupread.alignment.query_sequence[pileupread.query_position]

                            if other_mut in nt_to_variant:
                                for mut in nt_to_variant[other_mut]:
                                    if mut not in mutation_reads[pileupread.alignment.query_name][pos]:
                                        mutation_reads[pileupread.alignment.query_name][pos].extend([
                                                                                                    mut])
                            elif 'unknown' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['unknown'])
                    if pileupread.is_del and ref == 'del':
                        if pileupread.indel < 0:
                            if (int(alt) + int(pileupread.indel)) == 0:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(
                                    nt_to_variant[i])
                            else:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['unknown'])
                        elif pileupread.indel == 0:
                            if 'ref' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['ref'])
                        else:
                            if 'unknown' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['unknown'])

                    if pileupread.is_refskip and ref == 'ins':
                        if pileupread.indel > 0:
                            if (len(alt) - int(pileupread.indel)) == 0:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(
                                    nt_to_variant[i])
                            else:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['unknown'])
                        elif pileupread.indel == 0:
                            if 'ref' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['ref'])
                        else:
                            if 'unknown' not in mutation_reads[pileupread.alignment.query_name][pos]:
                                mutation_reads[pileupread.alignment.query_name][pos].extend(['unknown'])

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
            if 'ref' not in mutation_reads[name][pos] and 'unknown' not in mutation_reads[name][pos] and len(mutation_reads[name][pos]) == 1:
                if ''.join(mutation_reads[name][pos]) in parents_v:
                    parents_v[''.join(mutation_reads[name][pos])]['uniq'] += 1
                
                else:
                    parents_v[''.join(mutation_reads[name][pos])]['uniq'] = 0
                    parents_v[''.join(mutation_reads[name][pos])]['all'] = 0
                
    for v in parents_v:
        for name in mutation_reads:
            for pos in mutation_reads[name]:
                if 'ref' not in mutation_reads[name][pos] and v in mutation_reads[name][pos]:
                    parents_v[v]['all'] += 1

    parents_v = dict(sorted(parents_v.items(), key=lambda item: item[1]['all'], reverse=True))
    if (parents_v[list(parents_v.keys())[0]]['all'] < argvs.minReadCount or parents_v[list(parents_v.keys())[1]]['all'] < argvs.minReadCount):
        logging.error(f'Not have enough reads ({argvs.minReadCount}) to support both parents variants: {parents_v}')
        sys.exit(1)

    two_parents_list = list(parents_v.keys())[0:2]
    return two_parents_list

def find_recomb(mutation_reads, reads_coords, two_parents_list, reads_stats, argvs):
    logging.info("Finding recombinant reads")
    logging.info(f"Parents: {','.join(two_parents_list)} ")
    ## Can expand the list to check other recombinant???
    recomb_reads = []
    parent1_reads = []
    parent2_reads = []
    for read in mutation_reads:
        var1_count = 0
        var2_count = 0
        for k, v in sorted(mutation_reads[read].items()):
            var1 = two_parents_list[0]
            var2 = two_parents_list[1]
            if 'ref' in v:
                if var1 in v and var2 not in v:
                    var2_count += 1
                if var2 in v and var1 not in v:
                    var1_count += 1
            else:
                if var1 in v and var2 not in v:
                    var1_count += 1
                if var2 in v and var1 not in v:
                    var2_count += 1
        if var1_count > 0 and var2_count > 0 and read not in recomb_reads:
            #child
            recomb_reads.append(read)
        if var1_count == 0 and var2_count > 1 and read not in recomb_reads:
            parent2_reads.append(read)
        if var1_count > 1 and var2_count == 0 and read not in recomb_reads:
            parent1_reads.append(read)

    start_list = [reads_coords[i]['start'] for i in (recomb_reads + parent1_reads + parent2_reads)]
    end_list = [reads_coords[i]['end'] for i in (recomb_reads + parent1_reads + parent2_reads)]
    mutaions_list = []
    for read in (recomb_reads + parent1_reads + parent2_reads):
        tmp = dict()
        for k, v in sorted(mutation_reads[read].items()):
            if 'ref' in v:
                v = [" ".join(v).replace("ref", "ref of ")]

            tmp[k] = v
        mutaions_list.append(tmp)
    
    if argvs.ec19_projdir:
        rel_link = argvs.igv
        igv_list = [f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + '">recombinant</a>' for i in recomb_reads]
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + f'">parent {two_parents_list[0]}</a>' for i in parent1_reads])
        igv_list.extend([f'<a target="_new" href="{rel_link}?locus=NC_045512_2:' +
                    str(reads_coords[i]['start']) + "-" + str(reads_coords[i]['end']) + f'">parent {two_parents_list[1]}</a>' for i in parent2_reads])
        df = pd.DataFrame(list(zip((recomb_reads + parent1_reads + parent2_reads), start_list, end_list, mutaions_list, igv_list)), columns=[
            'read_name', 'start', 'end', 'mutaions_json', 'note'])
    else:
        note_list = [ 'recombinant' for i in recomb_reads]
        note_list.extend( [  'parent ' + two_parents_list[0]  for i in parent1_reads] )
        note_list.extend( [  'parent ' + two_parents_list[1]  for i in parent2_reads] )
        df = pd.DataFrame(list(zip((recomb_reads + parent1_reads + parent2_reads), start_list, end_list, mutaions_list, note_list)), columns=[
            'read_name', 'start', 'end', 'mutaions_json','note'])
    df.sort_values(by=['start'], inplace=True)
    df.to_csv(argvs.tsv, sep="\t",index=False)
    if argvs.ec19_projdir:
        html_file = os.path.splitext(argvs.tsv)[0] + ".html"
        df.to_html(html_file, index=False, escape=False)

    reads_stats['mutation_reads'] = len(mutation_reads)
    reads_stats['parents'] = ','.join(two_parents_list)
    reads_stats['recomb_reads'] = len(recomb_reads)
    reads_stats['parent1_reads'] = len(parent1_reads)
    reads_stats['parent2_reads'] = len(parent2_reads)
    total_recomb_paraents_reads = len(parent1_reads) + len(parent2_reads) + len(recomb_reads)
    reads_stats['recomb_perc'] = len(recomb_reads) /total_recomb_paraents_reads * 100 if total_recomb_paraents_reads > 0 else 0
    logging.info(f"Reads with Variants Mutations: {reads_stats['mutation_reads']}")
    logging.info(f"Recombinant Reads: {reads_stats['recomb_reads']} ({reads_stats['recomb_perc']:.2f}%)")
    logging.info(f"Parent 1 ({two_parents_list[0]}) Reads: {reads_stats['parent1_reads']}")
    logging.info(f"Parent 2 ({two_parents_list[1]}) Reads: {reads_stats['parent2_reads']}")

    return df, reads_stats

def write_recombinant_bam(recomb_reads, argvs):
    samfile = pysam.AlignmentFile(argvs.bam, "rb")
    logging.info("Index input Sam/BAM-file by query name. The index is kept in memory and can be substantial.")
    name_indexed = pysam.IndexedReads(samfile)
    name_indexed.build()
    tmp_bam = os.path.splitext(argvs.outbam)[0] +'.tmp.bam'
    if os.path.exists(tmp_bam):
        os.remove(tmp_bam)
    logging.info(f"Writting recombinant reads to {argvs.outbam}")
    outbam = pysam.AlignmentFile(tmp_bam, "wb", template=samfile)
    for name in recomb_reads['read_name']:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                outbam.write(x)
    outbam.close()
    samfile.close()
    pysam.sort("-o", argvs.outbam, tmp_bam)
    pysam.index(argvs.outbam)
    os.remove(tmp_bam)

def write_stats(stats, argvs):
    outstats = os.path.splitext(argvs.tsv)[0] + ".stats"
    logging.info(f"Writting recombinant stats to {outstats}")
    with open(outstats,"w") as of:
        of.write("\t".join(stats.keys()) + "\n")
        of.write(str("\t".join(str(x) for x in stats.values())) + "\n")


def update_igv_html(argvs):
    igv_html_file = os.path.join(os.path.dirname(argvs.bam), argvs.igv)
    if not os.path.exists(argvs.igv):
        return
    update_igv_html_file = os.path.splitext(igv_html_file)[0] + '.recombreads.html'
    logging.info(f"Writting recombinant reads as track in IGV view to {update_igv_html_file}")
    of = open(update_igv_html_file,'w')
    with open(igv_html_file, 'r') as f:
        for line in f:
            if 'recombinant_reads.bam' in line:
                of.close()
                if os.path.exists(update_igv_html_file):
                    os.remove(update_igv_html_file)
                return
            if 'options = {"tr' in line:
                options = [i.strip().replace(";","") for i in line.split("=")]
                options_dict = json.loads(options[1])
                new_tracks=[]
                recomb_track = {
                    'name': 'Recombinant', 
                    'type':'alignment', 
                    'format': 'bam', 
                    'colorBy': 'strand', 
                    'url': '../../ReadsBasedAnalysis/readsMappingToRef/recombinant_reads.bam', 
                    'indexURL': '../../ReadsBasedAnalysis/readsMappingToRef/recombinant_reads.bam.bai',
                    'squishedRowHeight': 10,
                    'height': 300,
                    'displayMode': 'SQUISHED' }
                for i in options_dict['tracks']:
                    if 'name' in i and i['name'] == "NC_045512.2 Alignment":
                        new_tracks.append(recomb_track)
                        new_tracks.append(i)
                    else:
                        if 'showGenotypes' not in i:
                            new_tracks.append(i)
                options_dict['tracks'] = new_tracks
                of.write("var options = " + json.dumps(options_dict) + ";")
            elif 'EC-19 Alignment' in line or 'NC_045512.2 Alignment' in line or "'Alignment'," in line:
                of.write("""
                    name: 'Recombinant', 
                    type:'alignment', 
                    format: 'bam', 
                    colorBy: 'strand', 
                    url: '../ReadsBasedAnalysis/readsMappingToRef/recombinant_reads.bam', 
                    indexURL: '../ReadsBasedAnalysis/readsMappingToRef/recombinant_reads.bam.bai',
                    squishedRowHeight: 10,
                    height: 300,
                    displayMode: 'SQUISHED' 
                    },
                    {
                    """ )
                of.write(line)
            else:
                of.write(line)
    of.close()
    shutil.move(update_igv_html_file,argvs.igv)

def main():
    argvs = setup_argparse()
    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
    format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')

    (delta_uniq_nt, omicron_uniq_nt, nt_to_variant,
     nt_to_aa) = load_var_mutation(argvs.variantMutation)
    nt_to_lineage = load_lineage_mutation(argvs.lineageMutation)
    if argvs.vcf:
        parents_variant, filtered_nt_to_variants=parse_vcf(argvs,nt_to_variant)
        mutation_reads, reads_coords, reads_stats = find_read_with_variants(
        filtered_nt_to_variants, argvs)
        two_parents_list = list(parents_variant.keys())[0:2]
    else:
        mutation_reads, reads_coords, reads_stats = find_read_with_variants(
        nt_to_variant, argvs)
        two_parents_list = find_parents_by_mutation_reads(mutation_reads,argvs)
        #two_parents_list = argvs.recombinant_variants

    recomb_reads_df,  reads_stats = find_recomb(mutation_reads, reads_coords, two_parents_list, reads_stats, argvs)
    
    write_stats(reads_stats,argvs)

    write_recombinant_bam(recomb_reads_df, argvs)
    if argvs.igv:
        update_igv_html(argvs)

    

if __name__ == '__main__':
    main()
