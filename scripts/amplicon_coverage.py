#!/usr/bin/env python3
import os, errno, sys
import argparse as ap
import numpy as np
import csv
import subprocess
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import pysam
from collections import defaultdict

# standardize the logging output
import logging
this_script=os.path.basename(__file__)
logging.basicConfig(format='['+this_script+'] %(levelname)s:%(message)s', level=logging.INFO)

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
    covGrp_me.add_argument('--bam', metavar='[FILE]', type=str,help='sorted bam file (ex: samtools sort input.bam -o sorted.bam)')
    covGrp_me.add_argument('--cov', metavar='[FILE]', type=str,help='coverage file [position\tcoverage]')
    
    outGrp = parser.add_argument_group('Output')
    outGrp.add_argument('-o', '--outdir', metavar='[PATH]',type=str, default='.', help='output directory')
    outGrp.add_argument('-p', '--prefix', metavar='[STR]',type=str , help='output prefix')

    #optGrp = parser.add_argument_group('Options')
    parser.add_argument('--pp', action='store_true', help='process proper paired only reads from bam file (illumina)')
    parser.add_argument('--mincov', metavar='[INT]', type=int, help='minimum coverage to count as ambiguous N site [default:10]', default=10)
    parser.add_argument('-r', '--refID', metavar='[STR]',type=str , help='reference accession (bed file first field')
    
    parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')
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

def covert_bed_to_amplicon_dict(input,cov_array,RefID="",unique=False):
    ## convert bed file to amplicon region dictionary
    input_bed = input
    cov_zero_array = np.zeros_like(cov_array)
    amplicon=defaultdict(dict)
    primers_pos=list()
    RefID = '""' if RefID is None else RefID
    cmd = 'grep -v alt %s | grep %s | paste - - | cut -f 2,3,4,8,9 | sed -e "s/_LEFT//g" -e "s/_RIGHT//g" ' % (input_bed, RefID)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    previous_id=''
    outs, errs = proc.communicate()
    if proc.returncode == 0:
    	outs = outs.splitlines()
    else:
    	logging.error("Failed %d %s %s" % (proc.returncode, outs, errs))
   
    for i in range(len(outs)):
        fstart, fend, id, rstart, rend = outs[i].decode().rstrip().split("\t")
        if id in amplicon:
            continue
        amplicon[id]=range(int(fend),int(rstart)+1)   
        primers_pos.extend(range(int(fstart),int(fend)+1))
        primers_pos.extend(range(int(rstart),int(rend)+1))
        for i in range(int(fend),int(rstart)+1):
            cov_zero_array[i] += 1

    if unique:
        unique_region_set = set(np.where(cov_zero_array == 1)[0]) - set(primers_pos)
        for i in range(len(outs)):
            fstart, fend, id, rstart, rend = outs[i].decode().rstrip().split("\t")
            unique_region = sorted(unique_region_set.intersection(set(range(int(fstart),int(rend)+1))))
            if unique_region:
                amplicon[id] = range(list(unique_region)[0],list(unique_region)[-1]+1)
            else:
                amplicon[id] = range(0)
         
    return amplicon

def covert_bedpe_to_amplicon_dict(input,cov_array,RefID="",unique=False):
    ## convert bed file to amplicon region dictionary 
    input_bedpe = input
    cov_zero_array = np.zeros_like(cov_array)
    amplicon=defaultdict(dict)
    primers_pos=list()
    RefID = '""' if RefID is None else RefID
    cmd = 'grep %s %s | cut -f 2,3,5,6,7 ' % (RefID, input_bedpe)
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE) 
    previous_id=''
    outs, errs = proc.communicate()
    if proc.returncode == 0:
        outs = outs.splitlines()
    else:
        logging.error("Failed %d %s %s" % (proc.returncode, outs, errs))

    for i in range(len(outs)):
        fstart, fend, rstart, rend, id = outs[i].decode().rstrip().split("\t")
        if id in amplicon:
            continue
        amplicon[id]=range(int(fend),int(rstart)+1)   
        primers_pos.extend(range(int(fstart),int(fend)+1))
        primers_pos.extend(range(int(rstart),int(rend)+1))
        for i in range(int(fend),int(rstart)+1):
            cov_zero_array[i] += 1
    
    if unique:
        unique_region_set = set(np.where(cov_zero_array == 1)[0]) - set(primers_pos)
        for i in range(len(outs)):
            fstart, fend, rstart, rend, id = outs[i].decode().rstrip().split("\t")
            unique_region = sorted(unique_region_set.intersection(set(range(int(fstart),int(rend)+1))))
            if unique_region:
                amplicon[id] = range(list(unique_region)[0],list(unique_region)[-1]+1)
            else:
                amplicon[id] = range(0)

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

def parse_bam_file(bam,pp,outdir,RefID):

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
        logging.info("Total Reads: %d" % (allcount))
        logging.info("Proper Paired Reads: %d" % (ppcount))
        if ppcount==0:
            logging.error("Failed: There is no proper paired reads in your input file %s." % (bam))
            sys.exit()
        
    cov_list = []
    for line in pysam.samtools.depth("-aa","-d0", bam_input ,split_lines=True):
        id, pos, cov = line.rstrip().split("\t")
        if RefID and RefID == id:
            cov_list.append(int(cov))
        else:
            cov_list.append(int(cov))
    cov_array=np.array(cov_list)
    return cov_array

def calculate_mean_cov_per_amplicon(cov_np_array,amplicon_d):
    mean_dict=dict()
    for key in amplicon_d:
        if amplicon_d[key]:
            start = amplicon_d[key][0]
            end = amplicon_d[key][-1]
            #print(key,start,end)
            mean_dict[key] = 0 if start > cov_np_array.size else cov_np_array[start:end].mean()
        else:
            # not have any range in this key region
            mean_dict[key] = 0 
    return mean_dict

def calculate_N_per_amplicon(cov_np_array,amplicon_d, mincov):
    num_low_depth=dict()
    for key in amplicon_d:
        if amplicon_d[key]:
            start = amplicon_d[key][0]
            end = amplicon_d[key][-1]
            if start > cov_np_array.size:
                num_low_depth[key] = 0
            else:
                arr = cov_np_array[start:end]
                num_low_depth[key] = np.size(np.where(arr < mincov))
        else:
            # not have any range in this key region
            num_low_depth[key] = 0
    return num_low_depth
    
def write_dict_to_file(mean_d,uniq_mean_d,ambiguity_d,uniq_ambiguity_d,outdir,prefix,mincov):
    output_amplicon_cov_txt = outdir + os.sep + prefix + '_amplicon_coverage.txt'
    try:
        os.remove(output_amplicon_cov_txt)
    except OSError:
        pass
    with open(output_amplicon_cov_txt,"w") as f:
        f.write(f"ID\tWhole_Amplicon\tUnique_Amplicon\tWhole_Amplicon_Ns(cov<{mincov})\tUnique_Amplicon_Ns(cov<{mincov})\n")
        for key in mean_d:
            f.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\n" %
              (
                key,
                mean_d[key],uniq_mean_d[key],
                ambiguity_d[key],uniq_ambiguity_d[key]
              )
            )
      

def barplot(mean_dict,uniq_mean_d,ambiguity_d,uniq_ambiguity_d,input_bed,overall_mean,outdir,prefix,mincov):
    #plot bar chart
    x_name=Path(input_bed).stem
    x=list(mean_dict.keys())
    y=list(mean_dict.values())
    uniq_x=list(uniq_mean_d.keys())
    uniq_y=list(uniq_mean_d.values())
    barcolor1 = ['lightsalmon' if i >= 20 else 'blue' if i >5 else 'black' for i in y]
    barcolor2 = ['lightsalmon' if i >= 20 else 'blue' if i >5 else 'black' for i in uniq_y]
    
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # coverage levels
    fig.add_trace(
      go.Bar(x=uniq_x,y=uniq_y,marker_color=barcolor2,visible=True,name="Coverage",showlegend=False),
      secondary_y=False
    )
    fig.add_trace(
      go.Bar(x=x, y=y, marker_color=barcolor1,visible=False,name="Coverage",showlegend=False),
      secondary_y=False
    )
    #fig = go.Figure(data=[go.Bar(x=x, y=y, marker_color=barcolor1,visible=False)])

    #fig.add_trace(go.Bar(x=uniq_x,y=uniq_y,marker_color=barcolor2,visible=True))

    # Add Ns
    y2=list(ambiguity_d.values())
    uniq_y2=list(uniq_ambiguity_d.values())
    fig.add_trace(
      go.Scatter(x=uniq_x,y=uniq_y2,marker_color='blueviolet',visible=True,name="Num of Ns (cov<" + str(mincov) +")" ),
      secondary_y=True
    )
    fig.add_trace(
      go.Scatter(x=x,y=y2,marker_color='blueviolet',visible=False,name="Num of Ns (cov<" + str(mincov) +")" ),
      secondary_y=True
    )
    
    depthMean = [dict(type='line',
                    xref='paper',x0=0,x1=0.94,
                    yref='y',y0=overall_mean,y1=overall_mean,
                    line=dict(
                        color="red",
                        width=1,
                        dash='dot',
                    )
                )]

    depth5x = [dict(type='line',
                    xref='paper',x0=0,x1=0.94,
                    yref='y',y0=5,y1=5,
                    line=dict(
                        color="black",
                        width=1,
                        dash='dot',
                    )
                )]
    depth10x = [dict(type='line',
                    xref='paper',x0=0,x1=0.94,
                    yref='10',y0=10,y1=10,
                    line=dict(
                        color="black",
                        width=1,
                        dash='dot',
                    )
                )]
    depth20x = [dict(type='line',
                    xref='paper',x0=0,x1=0.94,
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
                dict(label='Unique',
                     method='update',
                     args=[{"visible": [True, False, True, False]}]),
                dict(label='Whole Amplicon',
                     method='update',
                     args=[{"visible": [False, True, False, True]}])
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
    fig.update_yaxes(
        tickfont=dict(family='Courier New, monospace', size=8, color="blueviolet"),
        ticks="outside",
        tickwidth=0.5,
        ticklen=3,
        secondary_y=True,
        title_text="Num of Ns (cov<" + str(mincov) +")",
        titlefont=dict(
            color="blueviolet"
        ),
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
            xref='paper',x0=0,x1=0.94,
            yref='y',y0=overall_mean,y1=overall_mean,
            line=dict(
                color="red",
                width=1,
                dash='dashdot',
            ),
            ),
        ],
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=0.95
        )
    )

    output_html = outdir + os.sep + prefix + '_amplicon_coverage.html'
    fig.write_html(output_html)

def run(argvs):
    if (argvs.cov):
        cov_array = parse_cov_file(argvs.cov)
        prefix = Path(argvs.cov).stem
    if (argvs.bam):
        cov_array = parse_bam_file(argvs.bam,argvs.pp, argvs.outdir, argvs.refID)
        prefix = Path(argvs.bam).stem 
    if (argvs.bed):
        amplicon_dict = covert_bed_to_amplicon_dict(argvs.bed,cov_array,argvs.refID)
        uniq_amplicon_dict = covert_bed_to_amplicon_dict(argvs.bed,cov_array,argvs.refID,True)
        bedfile = argvs.bed
    if (argvs.bedpe):
        amplicon_dict = covert_bedpe_to_amplicon_dict(argvs.bedpe,cov_array,argvs.refID)
        uniq_amplicon_dict = covert_bedpe_to_amplicon_dict(argvs.bedpe,cov_array,argvs.refID,True)
        bedfile = argvs.bedpe
   
    
    if not argvs.prefix:
        argvs.prefix = prefix

    # Count the number of Ns < argvs.mincov per amplicon
    ambiguity_d = calculate_N_per_amplicon(cov_array,amplicon_dict, argvs.mincov)
    uniq_ambiguity_d = calculate_N_per_amplicon(cov_array,uniq_amplicon_dict,argvs.mincov)

    # Get the coverage per amplicon
    amplicon_mean_d = calculate_mean_cov_per_amplicon(cov_array,amplicon_dict)
    uniq_amplicon_mean_d = calculate_mean_cov_per_amplicon(cov_array,uniq_amplicon_dict)

    # Write results to a tab-delimited file
    write_dict_to_file(amplicon_mean_d,uniq_amplicon_mean_d,ambiguity_d,uniq_ambiguity_d,argvs.outdir,argvs.prefix,argvs.mincov)

    barplot(amplicon_mean_d,uniq_amplicon_mean_d,ambiguity_d,uniq_ambiguity_d,bedfile,cov_array.mean(),argvs.outdir,argvs.prefix,argvs.mincov)

if __name__ == '__main__':
    argvs = setup_argparse()
    mkdir_p(argvs.outdir)
    run(argvs)
