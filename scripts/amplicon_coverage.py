#!/usr/bin/env python3
import os, errno
import argparse as ap
import numpy as np
import csv
import subprocess
import plotly.graph_objects as go
from pathlib import Path
import pysam
from collections import defaultdict

class SmartFormatter(ap.HelpFormatter):
    def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            # this is the RawTextHelpFormatter._split_lines
            return ap.HelpFormatter._split_lines(self, text, width)

def setup_argparse():
    parser = ap.ArgumentParser(prog='amplicon_coverage.py',
        description = '''Script to parse amplicon region coverage and generate barplot in html''',
        formatter_class = SmartFormatter)

    inGrp = parser.add_argument_group('Amplicon Input (required, mutually exclusive)')
    inGrp_me = inGrp.add_mutually_exclusive_group(required=True)
    inGrp_me.add_argument('--bed', metavar='[FILE]', type=str,help='amplicon bed file')
    inGrp_me.add_argument('--bedpe', metavar='[FILE]', type=str,help='amplicon bedpe file')
   
    covGrp = parser.add_argument_group('Coverage Input (required, mutually exclusive)')
    covGrp_me = covGrp.add_mutually_exclusive_group(required=True)
    covGrp_me.add_argument('--bam', metavar='[FILE]', type=str,help='bam file')
    covGrp_me.add_argument('--cov', metavar='[FILE]', type=str,help='coverage file [position\tcoverage]')
    
    outGrp = parser.add_argument_group('Output')
    outGrp.add_argument('-o', '--outdir', metavar='[PATH]',type=str, default='.', help='output directory')
    outGrp.add_argument('-p', '--prefix', metavar='[STR]',type=str , help='output prefix')

    args_parsed = parser.parse_args()
    if not args_parsed.outdir:
        args_parsed.outdir = os.getcwd()

    return args_parsed

def mkdir_p(directory_name):
    try:
        os.makedirs(directory_name)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(directory_name):
            pass

def covert_bed_to_amplicon_dict(input):
    ## convert bed file to amplicon region dictionary
    input_bed = input
    amplicon=defaultdict(dict)
    cmd = 'grep -v alt %s | paste - - | awk \'{print $4"\t"$3"\t"$8}\' | sed -e "s/_LEFT//g" -e "s/_RIGHT//g" ' % (input_bed)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    for line in proc.stdout:
         id, start, end = line.decode().rstrip().split("\t")
         amplicon[id]['start']=start
         amplicon[id]['end']=end

    outs, errs = proc.communicate()
    if proc.returncode != 0:
        print("Failed %d %s %s" % (proc.returncode, outs, errs))

    return amplicon

def covert_bedpe_to_amplicon_dict(input):
    ## convert bed file to amplicon region dictionary 
    input_bedpe = input
    amplicon=defaultdict(dict)
    cmd = 'awk \'{print $7"\t"$3"\t"$5}\' %s ' % (input_bedpe)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    for line in proc.stdout:
         id, start, end = line.decode().rstrip().split("\t")
         amplicon[id]['start']=start
         amplicon[id]['end']=end
    outs, errs = proc.communicate()
    if proc.returncode != 0:
        print("Failed %d %s %s" % (proc.returncode, outs, errs))

    return amplicon

def parse_cov_file(input):
    ## read the coverage and save in a list
    cov_file = input
    cov_list = []
    with open(cov_file,'r') as cov_file:
        for line in cov_file:
            pos, cov = line.rstrip().split("\t")
            cov_list.append(int(cov))

    cov_array=np.array(cov_list)
    return cov_array

def parse_bam_file(bam):
    cov_list = []
    for line in pysam.samtools.depth("-aa","-d0", "NC_045512.2.sort_sorted.bam",split_lines=True):
        id, pos, cov = line.rstrip().split("\t")
        cov_list.append(int(cov))
    cov_array=np.array(cov_list)
    return cov_array

def calculate_mean_cov_per_amplicon(cov_np_array,amplicon_d):
    mean_dict=dict()
    for key in amplicon_d:
        start = amplicon_d[key]['start']
        end = amplicon_d[key]['end']
        mean_dict[key] = cov_np_array[int(start):int(end)].mean()
    return mean_dict

def barplot(mean_dict,input_bed,overall_mean,outdir,prefix):
    #plot bar chart
    x_name=Path(input_bed).stem
    x=list(mean_dict.keys())
    y=list(mean_dict.values())
    fig = go.Figure(data=[go.Bar(x=x, y=y)])

    updatemenus = list([
        dict(
             buttons=list([
                dict(label='Linear Scale',
                     method='relayout',
                     args=[{'title': '',
                           'yaxis.type': 'linear'}]),
                dict(label='Log Scale',
                     method='relayout',
                     args=[{'title': '',
                            'yaxis.type': 'log'}])
                ]),
             direction="down",
             x=0.1,
             xanchor="left",
             y=1.1,
             yanchor="top"
            )
        ])

    fig.update_layout(
        updatemenus=updatemenus,
        xaxis_title=x_name,
        yaxis_title="Mean Coverage(X)",
        font=dict(
            family="Courier New, monospace",
            size=10,
        ),
        annotations=[
            dict(text="Y-axis:", showarrow=False,
            x=0, y=1.085, yref="paper", align="left")
        ],
        shapes=[
    	    dict(type='line',
    		xref='paper',x0=0,x1=1,
    		yref='y',y0=overall_mean,y1=overall_mean,
    		line=dict(
    			color="red",
    			width=1,
    			dash='dashdot',
    		),
    	    ),
        ]
    )

    output_html = outdir + os.sep + prefix + '_amplicon_coverage.html'
    fig.write_html(output_html)

def run(argvs):
    if (argvs.bed):
        amplicon_dict = covert_bed_to_amplicon_dict(argvs.bed)
        bedfile = argvs.bed
    if (argvs.bedpe):
        amplicon_dict = covert_bedpe_to_amplicon_dict(argvs.bedpe)
        bedfile = argvs.bedpe
    if (argvs.cov):
        cov_array = parse_cov_file(argvs.cov)
        prefix = Path(argvs.cov).stem 
    if (argvs.bam):
        cov_array = parse_bam_file(argvs.bam)
        prefix = Path(argvs.bam).stem
    if not argvs.prefix:
        argvs.prefix = prefix
        
    mean_d = calculate_mean_cov_per_amplicon(cov_array,amplicon_dict)
    barplot(mean_d,bedfile,cov_array.mean(),argvs.outdir,argvs.prefix)

if __name__ == '__main__':
    argvs = setup_argparse()
    mkdir_p(argvs.outdir)
    run(argvs)