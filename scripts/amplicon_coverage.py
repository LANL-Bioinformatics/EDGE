#!/usr/bin/env python3
import os, errno, sys
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
    parser = ap.ArgumentParser(prog='amplicov',
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

    #optGrp = parser.add_argument_group('Options')
    parser.add_argument('--pp', action='store_true', help='process proper paired only reads from bam file (illumina)')
    
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.2')
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

def covert_bed_to_amplicon_dict(input,unique=False):
    ## convert bed file to amplicon region dictionary
    input_bed = input
    amplicon=defaultdict(dict)
    if unique:
       cmd = 'grep -v alt %s | paste - - | awk \'{print $4"\t"$2"\t"$9}\' | sed -e "s/_LEFT//g" -e "s/_RIGHT//g" ' % (input_bed)
    else:
       cmd = 'grep -v alt %s | paste - - | awk \'{print $4"\t"$3"\t"$8}\' | sed -e "s/_LEFT//g" -e "s/_RIGHT//g" ' % (input_bed)
 
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    previous_id=''
    for line in proc.stdout:
         id, start, end = line.decode().rstrip().split("\t")
         amplicon[id]=range(int(start),int(end)+1)
         if (previous_id and unique):
            set_new=set(amplicon[id])
            set_previous=set(amplicon[previous_id])
            set_new_update=sorted(set_new - set_previous)
            set_previous_update=sorted(set_previous - set_new)
            amplicon[id] = range(set_new_update[0],set_new_update[-1]+1)
            amplicon[previous_id] = range(set_previous_update[0],set_previous_update[-1]+1)
         previous_id = id

    outs, errs = proc.communicate()
    if proc.returncode != 0:
        print("Failed %d %s %s" % (proc.returncode, outs, errs))

    return amplicon

def covert_bedpe_to_amplicon_dict(input,unique=False):
    ## convert bed file to amplicon region dictionary 
    input_bedpe = input
    amplicon=defaultdict(dict)
    if unique:
        cmd = 'awk \'{print $7"\t"$2"\t"$6}\' %s ' % (input_bedpe)
    else:
        cmd = 'awk \'{print $7"\t"$3"\t"$5}\' %s ' % (input_bedpe)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    previous_id=''
    for line in proc.stdout:
         id, start, end = line.decode().rstrip().split("\t")
         amplicon[id]=range(int(start),int(end)+1)
         if (previous_id and unique):
            set_new=set(amplicon[id])
            set_previous=set(amplicon[previous_id])
            set_new_update=sorted(set_new - set_previous)
            set_previous_update=sorted(set_previous - set_new)
            amplicon[id] = range(set_new_update[0],set_new_update[-1]+1)
            amplicon[previous_id] = range(set_previous_update[0],set_previous_update[-1]+1)
         previous_id = id

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

def parse_bam_file(bam,pp,outdir):

    bam_input = bam
    if pp:
        bam_prefix = Path(bam).stem
        bam_input = outdir + os.sep + bam_prefix + ".pp.bam"
        samfile = pysam.AlignmentFile(bam, "rb")
        ppcount=0;
        allcount=0;
        proper_pairedreads = pysam.AlignmentFile(bam_input, "wb", template=samfile)
        for read in samfile.fetch(until_eof=True):
            allcount+=1
            if read.is_proper_pair:
                ppcount+=1;
                proper_pairedreads.write(read)

        proper_pairedreads.close()
        samfile.close()
        print("Total Reads: %d" % (allcount))
        print("Proper Paired Reads: %d" % (ppcount))
        if ppcount==0:
            print("Failed: There is no proper paired reads in your input file %s." % (bam))
            sys.exit()
        
    cov_list = []
    for line in pysam.samtools.depth("-aa","-d0", bam_input ,split_lines=True):
        id, pos, cov = line.rstrip().split("\t")
        cov_list.append(int(cov))
    cov_array=np.array(cov_list)
    return cov_array

def calculate_mean_cov_per_amplicon(cov_np_array,amplicon_d):
    mean_dict=dict()
    for key in amplicon_d:
        start = amplicon_d[key][0]
        end = amplicon_d[key][-1]
        mean_dict[key] = cov_np_array[start:end].mean()
    return mean_dict
    
def write_dict_to_file(mean_d,uniq_mean_d,outdir,prefix):
    output_amplicon_cov_txt = outdir + os.sep + prefix + '_amplicon_coverage.txt'
    try:
        os.remove(output_amplicon_cov_txt)
    except OSError:
        pass
    with open(output_amplicon_cov_txt,"w") as f:
        f.write("ID\tWhole_Amplicon\tUnique\n")
        for key in mean_d:
            f.write("%s\t%.2f\t%.2f\n" % (key,mean_d[key],uniq_mean_d[key]))
      

def barplot(mean_dict,uniq_mean_d,input_bed,overall_mean,outdir,prefix):
    #plot bar chart
    x_name=Path(input_bed).stem
    x=list(mean_dict.keys())
    y=list(mean_dict.values())
    uniq_x=list(uniq_mean_d.keys())
    uniq_y=list(uniq_mean_d.values())
    fig = go.Figure(data=[go.Bar(x=x, y=y)])
    fig.add_trace(go.Bar(x=uniq_x,y=uniq_y,marker_color='lightsalmon',visible=False))
    
    depthMean = [dict(type='line',
                    xref='paper',x0=0,x1=1,
                    yref='y',y0=overall_mean,y1=overall_mean,
                    line=dict(
                        color="red",
                        width=1,
                        dash='dot',
                    )
                )]

    depth5x = [dict(type='line',
                    xref='paper',x0=0,x1=1,
                    yref='y',y0=5,y1=5,
                    line=dict(
                        color="black",
                        width=1,
                        dash='dot',
                    )
                )]
    depth10x = [dict(type='line',
                    xref='paper',x0=0,x1=1,
                    yref='10',y0=10,y1=10,
                    line=dict(
                        color="black",
                        width=1,
                        dash='dot',
                    )
                )]
    depth20x = [dict(type='line',
                    xref='paper',x0=0,x1=1,
                    yref='y',y0=20,y1=20,
                    line=dict(
                        color="black",
                        width=1,
                        dash='dot',
                    )
                )]

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
             x=0,
             xanchor="left",
             y=1.03,
             yanchor="top"
        ),
        dict(
            buttons=list([
                dict(label='Whole Amplicon',
                     method='update',
                     args=[{"visible": [True, False]}]),
                dict(label='Unique',
                     method='update',
                     args=[{"visible": [False, True]}])
                
                ]),
             direction="down",
             x=0.15,
             xanchor="left",
             y=1.03,
             yanchor="top"
        ),
        dict(
            buttons=list([
                dict(label="Mean(" + str(int(overall_mean)) + 'X)',
                 method="relayout",
                 args=["shapes", depthMean]),
                dict(label="5X",
                 method="relayout",
                 args=["shapes", depth5x]),
                dict(label="10X",
                 method="relayout",
                 args=["shapes", depth10x]),
                dict(label="20X",
                 method="relayout",
                 args=["shapes", depth20x]),
                dict(label="All",
                 method="relayout",
                 args=["shapes", depthMean + depth5x + depth10x + depth20x])
                ]),
             direction="down",
             x=0.30,
             xanchor="left",
             y=1.03,
             yanchor="top"
        )
    ])

    fig.update_xaxes(
        tickfont=dict(family='Courier New, monospace', size=8),
        ticks="outside",
        tickwidth=0.5,
        ticklen=3
    )
    fig.update_yaxes(
        tickfont=dict(family='Courier New, monospace', size=8),
        ticks="outside",
        tickwidth=0.5,
        ticklen=3
    )

    fig.update_layout(
        updatemenus=updatemenus,
        xaxis_title=x_name,
        yaxis_title="Mean Coverage(X)",
        font=dict(
            family="Courier New, monospace",
            size=10,
        ),
        annotations=[
            dict(text="Y-axis:",x=0, xref="paper",
                                y=1.06, yref="paper", 
                                showarrow=False),
            dict(text="Region:",x=0.15, xref="paper",
                                y=1.06, yref="paper", 
                                showarrow=False),
            dict(text="DepthLine:", x=0.30, xref="paper",
                                y=1.06, yref="paper", 
                                showarrow=False),
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
        uniq_amplicon_dict = covert_bed_to_amplicon_dict(argvs.bed,True)
        bedfile = argvs.bed
    if (argvs.bedpe):
        amplicon_dict = covert_bedpe_to_amplicon_dict(argvs.bedpe)
        uniq_amplicon_dict = covert_bedpe_to_amplicon_dict(argvs.bedpe)
        bedfile = argvs.bedpe
    if (argvs.cov):
        cov_array = parse_cov_file(argvs.cov)
        prefix = Path(argvs.cov).stem 
    if (argvs.bam):
        cov_array = parse_bam_file(argvs.bam,argvs.pp, argvs.outdir)
        prefix = Path(argvs.bam).stem
    if not argvs.prefix:
        argvs.prefix = prefix
        
    amplicon_mean_d = calculate_mean_cov_per_amplicon(cov_array,amplicon_dict)
    uniq_amplicon_mean_d = calculate_mean_cov_per_amplicon(cov_array,uniq_amplicon_dict)
    write_dict_to_file(amplicon_mean_d,uniq_amplicon_mean_d,argvs.outdir,argvs.prefix)
    barplot(amplicon_mean_d,uniq_amplicon_mean_d,bedfile,cov_array.mean(),argvs.outdir,argvs.prefix)

if __name__ == '__main__':
    argvs = setup_argparse()
    mkdir_p(argvs.outdir)
    run(argvs)
