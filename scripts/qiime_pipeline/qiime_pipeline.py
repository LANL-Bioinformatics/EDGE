#!/usr/bin/env python
__version__ = "0.3.0"
__author__    = "Chienchi Lo, Bioscience Division, Los Alamos National Laboratory"
__date__      = "2018/09/20"
__copyright__ = "BSD 3-Clause"

import sys
import os
import errno
import argparse
import subprocess
import shutil
import io
import datetime
import time
import platform
import gzip
import glob
import re
import numpy as np
import pandas as pd
from collections import defaultdict, Counter


script_path = os.path.dirname(os.path.abspath(__file__))

target_path = script_path + '/data'

# ME
#print("The script path is: {}".format(script_path))
#print("The target path is: {}".format(target_path))

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
	    if text.startswith('R|'):
	        return text[2:].splitlines()  
	    # this is the RawTextHelpFormatter._split_lines
	    return argparse.HelpFormatter._split_lines(self, text, width)

def setup_argparse():
    parser = argparse.ArgumentParser(prog='qiime_pipeline.py',
        description = '''Script to run Qiime2 pipeline based on Moving Pictures tutorial''',
        formatter_class = SmartFormatter)

    metaGrp = parser.add_argument_group('Mapping File (required)')
    metaGrp.add_argument('-m', '--mappingFile',required=True, type=str,help='MAPPING file.  #SampleID BarcodeSequence LinkerPrimerSequence sampleType ... Description')   

    readGrp = parser.add_argument_group('Reads Input (required, mutually exclusive)')
    eg = readGrp.add_mutually_exclusive_group(required=True)
    eg.add_argument('-p', '--paired', metavar='[FASTQ]',nargs=2,type=str,help='Paired reads in two fastq files and separate by space')
    eg.add_argument('-s', '--single', metavar='[FASTQ]',nargs=1,type=str,help='Unpaired reads fastq')
    eg.add_argument('-d', '--dir', metavar='[PATH]',type=str,help='''Contains Demultiplexed fastq files. To use this option.
                            To use this option,the mapping file need a 'Files' column with filenames for each sampleID.''')
    bcGrp = parser.add_argument_group('Barcode Fastq (requried when -p or -s)')
    bcGrp.add_argument('-b', '--barcode', metavar='[FASTQ]',nargs=1,type=str,help='Barcodes fastq')
    
    outGrp = parser.add_argument_group('Output')
    outGrp.add_argument('-o', '--outdir', metavar='[PATH]',type=str, default='.', help='Output directory')
    outGrp.add_argument('--zip',action='store_true', help='Zip output files for visualization')

    optGrp = parser.add_argument_group('Parameters')
    optGrp.add_argument('--target', metavar='<STR>', default='Greengenes', type=str, help='Greengenes, SILVA, SILVA-V3-V4 or ITS. [default: Greengenes]\nGreengenes and SILVA-V3-V4 are for 16s. SILVA is for 16s/18s. ITS is from https://unite.ut.ee/ and for fungal rDNA ITS sequences.\n## the trained qza files in data/')
    optGrp.add_argument('--qcMethod', metavar='<STR>', default='dada2', type=str, help='Quality control method. dada2 or deblur or otus[default: dada2]')
    
    # ME
    # optGrp.add_argument('--singleOrPairedQC', metavar='<STR>', default='single', type=str, help = 'Choose whether to use single or paired QC. single or paired [defaule: single]')
    optGrp.add_argument('--trimLeftForward', metavar='<INT>', default=20, type=int, help = 'This is for Dada2 QC on paired end reads. The number of bases to trim from the left of the forward (R1) read to e.g. remove PCR primer sequences. [default: 20]')
    optGrp.add_argument('--trimLeftReverse', metavar='<INT>', default=20, type=int, help = 'This is for Dada2 QC on paired end reads. This is the number of bases to trim from the left of the reverse (R2) read to e.g. remove PCR primer sequences. [default: 20]')    
    optGrp.add_argument('--truncLenForward', metavar='<INT>', default=0, type=int, help = 'This is for Dada2 QC on paired end reads. This is the truncation length of the forward reads after any trimming. "0" is no truncation. [default: 0]')    
    optGrp.add_argument('--truncLenReverse', metavar='<INT>', default=0, type=int, help = 'This is for Dada2 QC on paired end reads. This is the truncation length of the reverse reads after any trimming. "0" is no truncation. [default: 0]')    
    
    
    optGrp.add_argument('--trimLen', metavar='<INT>', default=20, type=int, help = 'This is for Dada2 QC on single end reads. This is the number of bases to trim from the reads to e.g. remove PCR primer sequences.[default: 20]')    
    optGrp.add_argument('--truncLen', metavar='<INT>', default=0, type=str, help='This works for Dada2 and Deblur on single end reads to truncate sequences at position [default: 0] no truncate')
    optGrp.add_argument('--minQuality',metavar='<INT>',default=4, type=int,help='This is for Deblur and OTUs QC. The minimum acceptable PHRED score. All PHRED scores less that this value are considered to be low PHRED scores. [default: 4]' )
    optGrp.add_argument('--minLengthFraction',metavar='<FLOAT>',default=0.75,type=float, help='This is for Deblur and OTUs QC. The minimum length that a sequence read can be following truncation and still be retained. This length should be provided as a fraction of the input sequence length. [default: 0.75]')
    optGrp.add_argument('--maxAmbiguous',metavar='<INT>',type=int,default=0,help='This is for Deblur and OTUs QC. The maximum number of ambiguous (i.e., N) base calls. This is applied after trimming sequences based on `min_length_fraction. [default: 0]')
    optGrp.add_argument('--percIdentity',metavar='<FLOAT>',type=float,default=0.99,help='This is for OTUs clustering. The percent identity at which clustering should be performed. This parameter maps to vsearch\'s --id parameter.  [default: 0.99]')
    #optGrp.add_argument('--minOTUsize', metavar='<INT>', default=2, type=int, help='the minimum OTU size (in number of sequences) to retain the OTU')
    optGrp.add_argument('--samplingDepth', metavar='<INT>', default=1000, type=int, help='Filter sample less this amount of sequences.The minimium of sequenceing depth of samples after this filter will be  use for even sub-sampling and maximum rarefaction depth.')
    optGrp.add_argument('--autoDepth', action='store_true', help='Automatically adjust the sampling to the minimum sequences count of all samples. The minimum > samplingDepth option above.')
    optGrp.add_argument('--maxRarefactionDepth', metavar='<INT>', default='0', type=int, help='The maximum rarefaction depth.[default: same as samplingDepth]')
    optGrp.add_argument('--phred_offset', metavar='<INT>', default=33, type=int, help='The ascii offset to use when decoding phred scores')
    optGrp.add_argument('-c', '--cpus', metavar='<INT>', default=8, type=int, help='Number of CPUS')
    optGrp.add_argument('--title', metavar='<STR>', default='', type=str, help='Project Title')
    
    parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
    #parser.add_argument('--quiet', action='store_true', help='Keep messages in terminal minimal')
    return parser.parse_args()

def mkdir_p(directory_name):
    try:
        os.makedirs(directory_name)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(directory_name):
            pass

def symlink_force(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)
        else:
            raise e

def dependency_check(cmd):
	proc = subprocess.Popen("which " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	outs, errs = proc.communicate()
	return outs.decode().rstrip() if proc.returncode == 0 else False

def tool_version_check(cmd):
	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	outs, errs = proc.communicate()
	return outs.decode().rstrip() if proc.returncode == 0 else False

def get_runtime(start):
    end = time.time()
    hours, rem = divmod(float(end-start), 3600)
    minutes, seconds = divmod(rem, 60)
    runtime ="Running time: {:0>2}:{:0>2}:{:0>2}".format(int(hours),int(minutes), int(seconds))
    return runtime

def process_cmd(cmd, msg=''):
    if msg:
        print("");
        print("###########################################################################")
        print("Qiime [%s] %s" % ( datetime.datetime.now().strftime('%Y %b %d %H:%M:%S'), msg) )
        print("###########################################################################")
    
    print("Qiime CMD: %s" %(cmd) )
    
    start=time.time()

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    if proc.returncode != 0: 
        print("Failed %d %s %s" % (proc.returncode, outs, errs))

    if msg:
        print("Qiime %s" % get_runtime(start))

def qiime_export_html(input,outdir):
    html_cmd = ('qiime tools export --input-path '+ input + ' --output-path ' + outdir)
    process_cmd(html_cmd)

def auto_determine_sampling_depth(file,fixedDepth):
    f =  open (file, "r")
    depth_list=list()
    depth=0
    for line in f:
        if not line.strip():continue
        sid, num =line.split(',')
        depth_list.append(float(num))
    f.close
    
    depth_list.sort()
    if (depth_list[0] > fixedDepth):
        depth = depth_list[0]
    if (depth_list[-1] < fixedDepth):
        depth = depth_list[-1]

    median = np.median(depth_list)
    
    return depth, median, depth_list

def auto_determine_truncate_len_deblur(file):
    f =  open (file, "r")
    len=-1
    for line in f:
        if not line.strip():continue
        if '2%' in line: 
            line = next(f) 
            match = re.search(r"\d+",line)
            len = int(match.group())
    f.close
    
    return len

def fix_mappingFile(mappingFile,out_dir):
#name reserved for the ID column header. 
#Reserved ID column headers:\n\n  
#Case-insensitive: 'feature id', 'feature-id', 'featureid', 'id', 'sample id', 'sample-id', 'sampleid'\n\n  
#Case-sensitive: '#OTU ID', '#OTUID', '#Sample ID', '#SampleID', 'sample_name'\n\n
    fix_f = out_dir + '/sample_metadata.txt'
    f =  open (mappingFile, "r")
    ff = open (fix_f, 'w')
    category_list=list()
    fix_metadata_list=list()
    non_category_list=['#SampleID','Barcode','Linker','Day','Description','Files']
    num_sample = 0
    cat_dict=defaultdict(Counter)
    for line in f:
        if not line.strip():continue
        temp = [ x.strip(' "\'\n') for x in line.strip().split('\t')]
        temp[0] = temp[0].replace('.', '-').replace('_', '-').replace(' ','-')
        if line.lower().startswith('#'):
            line = line.strip().replace('sample_name', 'Sample_Name')
            temp = [ x.strip(' "\'\n').replace(' ','-') for x in line.split('\t')]
            header = temp
        else:
            if len(temp) > len(header):
                sys.exit("[Error] Metadata row contains more cells than are declared by the header. Row detail: %s" % line)
            if len(temp) < len(header):
                # fill short temp list with NaN
                temp.extend( ['NaN'] * (len(header) - len(temp)) )
            for i, e in enumerate(temp):
                if header[i] not in non_category_list:
                    cat_dict[header[i]][temp[i]]=1
                    
        fix_metadata_list.append("\t".join(temp))
        if not line.lower().startswith('#'):
            num_sample += 1
    for x in header:
        if x not in non_category_list:
            try:
                ## Assume all numeric is not category
                if (all( x.isnumeric() for x in cat_dict[x])):
                    continue
            except:
                pass
            ## for Beta diversity test, each group should have > 1 samples to get within sample distances ??
            ##  https://forum.qiime2.org/t/error-on-beta-group-significance/1789/6
            ## if len(cat_dict[x]) > 1 and len(cat_dict[x]) < num_sample and all(v > 1 for v in cat_dict[x].values()):
            if len(cat_dict[x]) > 1 and len(cat_dict[x]) < num_sample:
                category_list.append(x)        
    f.close    

    for line in fix_metadata_list:
        temp = line.split('\t')
        if line.lower().startswith('#'):
            header = temp
        else:
            for i, e in enumerate(temp):
                if header[i] in category_list:
                    temp[i]=re.sub('[^0-9a-zA-Z]+', '-',temp[i])
        new_line = "\t".join(temp)
        ff.write(new_line + "\n")
        
    ff.close
    os.remove(mappingFile)
    
    return fix_f , num_sample, category_list

def get_fastq_manifest(mappingFile,in_dir,out_dir):
    mf_file = out_dir + '/file_manifest.txt'
    f =  open (mappingFile, "r")
    mf = open (mf_file, "w")
    abs_inDir = os.path.abspath(in_dir)
    mf.write("sample-id,absolute-filepath,direction\n")
    f_index = 0
    for line in f:
        if not line.strip():continue
        temp = line.strip().split('\t')
        if line.lower().startswith('#'):
            if 'files' in line.lower():
                header = line.lower().strip().split('\t')
                f_index = header.index('files')
        else:
            if (f_index > 0):
                if (',' in temp[f_index]):
                    f_fq,r_fq = temp[f_index].split(',')
                    type='pe'
                else:
                    f_fq = temp[f_index]
                    type='se'
            else:
                sys.exit("[ERROR] 'Files' column not found in meta data mapping file.")
   
            f_fq = f_fq.replace('"', '').strip()
            mf.write('%s,%s/%s,%s\n' % (temp[0],abs_inDir,f_fq,'forward'))
            symlink_force(abs_inDir+'/'+ f_fq,out_dir + '/input/' + os.path.basename(f_fq))
            if (type == 'pe'):
                r_fq = r_fq.replace('"', '').strip()
                mf.write('%s,%s/%s,%s\n' % (temp[0],abs_inDir,r_fq,'reverse'))
                symlink_force(abs_inDir+'/'+ r_fq,out_dir + '/input/' + os.path.basename(r_fq))

    f.close
    mf.close
    return type , mf_file       
            

def add_html_table(analysis_type,category=list()):
    tab_list=list()
    src_list=list()
    beta_emperor_result = ['bray_curtis','jaccard','unweighted_unifrac','weighted_unifrac']
    beta_emperor_title = ['Bray-Curtis','Jaccard','Unweighted Unifrac','Weighted Unifrac']
    alpha_analysis_result = ['evenness_vector','faith_pd_vector','observed_otus_vector','shannon_vector']
    alpha_analysis_title = ['Evenness','Faith_pd','Observed OTUs','Shannon']
    activate_class_count = 0 
    if ( analysis_type == 'beta' ):
        frame_height = 600
        out_html = 'DiversityAnalysis/table.html'
        for i, val in enumerate(beta_emperor_result):
            if os.path.isfile("DiversityAnalysis/%s_emperor/emperor.html" % val ):
                activate_class_count += 1
                activate_class = "active" if activate_class_count < 2 else ""
                tab_list.append('<li id="%s_pcoa" class="%s"><a href="#">%s</a></li>' % (val, activate_class, beta_emperor_title[i]))
                src_list.append("'%s_pcoa':'./%s_emperor/emperor.html'" % (val ,val))
    elif ( analysis_type == 'beta-significance'):
        frame_height = 1400
        out_html = 'DiversityAnalysis/beta-table.html'
        if ( len(category)>0):
            for g in category:
                for i, val in enumerate(beta_emperor_result):
                    if os.path.isfile("DiversityAnalysis/%s_%s_boxplot/index.html" % (val,g)):
                        activate_class_count += 1 
                        activate_class = "active" if activate_class_count < 2 else ""
                        tab_list.append('<li id="%s_%s_boxplot" class="%s"><a href="#">%s-%s</a></li>' % (val,g,activate_class, g, beta_emperor_title[i]))
                        src_list.append("'%s_%s_boxplot':'./%s_%s_boxplot/index.html'" % (val ,g ,val, g))         
    elif ( analysis_type == 'alpha' ):
        frame_height = 1400
        out_html = 'DiversityAnalysis/alpha-table.html'
        for i, val in enumerate(alpha_analysis_result):
            if os.path.isfile("DiversityAnalysis/%s_boxplot/index.html" % val ):
                activate_class_count += 1
                activate_class = "active" if activate_class_count < 2 else ""
                tab_list.append('<li id="%s_boxplot" class="%s"><a href="#">%s</a></li>' % (val, activate_class, alpha_analysis_title[i]))
                src_list.append("'%s_boxplot':'./%s_boxplot/index.html'" % (val ,val))
    elif ( analysis_type == 'taxonomy-ancom' ):
        frame_height = 2000
        out_html = 'TaxonomyAnalysis/ancom.html'
        if ( len(category)>0):
            for g in category:
                if os.path.isfile("TaxonomyAnalysis/l6-ancom-%s/index.html" % (g)):
                    activate_class_count += 1
                    activate_class = "active" if activate_class_count < 2 else ""
                    tab_list.append('<li id="l6-ancom-%s" class="%s"><a href="#">%s</a></li>' % (g,activate_class, g))
                    src_list.append("'l6-ancom-%s':'./l6-ancom-%s/index.html'" % (g,g))
  
    first_src_file = src_list[0].split(':')[1] if len(src_list) > 0 else ""
    html="""<!DOCTYPE html>
<html>
  <head>
    <meta charset='utf-8'>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <title>q2 beta diversity</title>
    <link rel="stylesheet" href="../q2templateassets/css/bootstrap.min.css"/>
    <link rel="stylesheet" href="../q2templateassets/css/normalize.css"/>
    <link rel="stylesheet" href="../q2templateassets/css/tab-parent.css">
    <script src="../q2templateassets/js/jquery-3.2.0.min.js" charset="utf-8"></script>
  </head>
  <body>
    <div class="container-fluid">
        <div class="row">
        <ul class="nav nav-tabs">
          %s
        </ul>
        </div>
        <iframe id="tab-frame" src=%s height="%dpx"></iframe>
        <!--<div class="row">
            <div class="col-lg-12">
            <h1>Beta Diversity PCoA Plots</h1>
            <table class="table table-striped table-hover">
            <tr><td>Bray-Curtis</td><td><a href="bray_curtis_emperor/emperor.html" target="_new">html</a></td></tr>
            <tr><td>Jaccard</td><td><a href="jaccard_emperor/emperor.html" target="_new">html</a></td></tr>
            <tr><td>unweighted UniFrac</td><td><a href="unweighted_unifrac_emperor/emperor.html" target="_new">html</a></td></tr>
            <tr><td>weighted UniFrac</td><td><a href="weighted_unifrac_emperor/emperor.html" target="_new">html</a></td></tr>
            </table>
            </div>
        </div>-->
        
    </div>
    <script type="text/javascript">
      //if (window.frameElement) {
      //if (window.self != window.top) {
      //    document.getElementById('q2templatesheader').remove();
      //}
      var frameSrc = {
          %s
      };
      var frame;

      function frameLoad(event) {
        frame.height = `${event.data + 50}px`;
      }

      function toggleClass() {
        if (document.querySelector('.active').id != this.id) {
            document.querySelector('.active').className = '';
            this.className = 'active';
            frame.height = '%dpx';
            frame.contentWindow.location.replace(frameSrc[this.id]);
        }
      }

      function init() {
        for (const frame of Object.keys(frameSrc)) {
            document.getElementById(frame).addEventListener('click', toggleClass);
        }
        frame = document.getElementById('tab-frame');
      }

      document.addEventListener('DOMContentLoaded', init);
      //window.addEventListener('message', frameLoad);
    </script>
  </body>
</html>""" % ("\n".join(tab_list), first_src_file , frame_height, ",\n".join(src_list), frame_height)

    if len(tab_list) > 0 :
        with open( out_html, "w") as f:
            f.write(html)
            f.close()


def html_report(template):
    ## need add for each index.html if not exist<script src="./q2templateassets/js/child.js"></script>
    # javascript modified if (window.self != window.top) {
    if not os.path.exists('q2templateassets'):
        shutil.copytree(template,'q2templateassets')

    tab_list=list()
    src_list=list()
    if os.path.isfile('input/index.html'):
        tab_list.append('<li id="sampleMetadata"><a href="#">Sample Metadata</a></li>')
        src_list.append("'sampleMetadata': './input/index.html'")
    if os.path.isfile('demux/overview.html'):
        tab_list.append('<li id="overview" class="active"><a href="#">Input Overview</a></li>')
        src_list.append("'overview': './demux/overview.html'")
    if os.path.isfile('demux/quality-plot.html'):
        tab_list.append('<li id="quality-plot"><a href="#">Reads Quality Plot</a></li>')
        src_list.append("'quality-plot': './demux/quality-plot.html'")
    if os.path.isfile('demux-joined/overview.html'):
        tab_list.append('<li id="joined-overview"><a href="#">Joined Reads Overview</a></li>')
        src_list.append("'joined-overview': './demux-joined/overview.html'")
    if os.path.isfile('demux-joined/quality-plot.html'):
        tab_list.append('<li id="joined-quality-plot"><a href="#">Joined Reads Quality Plot</a></li>')
        src_list.append("'joined-quality-plot': './demux-joined/quality-plot.html'")
    if os.path.isfile('QCandFT/filter-stats/index.html'):
        tab_list.append('<li id="filter-stats"><a href="#">Filter QC Stats</a></li>')
        src_list.append("'filter-stats': './QCandFT/filter-stats/index.html'")
    #if os.path.isfile('QCandFT/chimera-stats/index.html'):
    #    tab_list.append('<li id="chimera-stats"><a href="#">Chimera QC Stats</a></li>')
    #    src_list.append("'chimera-stats': './QCandFT/chimera-stats/index.html'")
    if os.path.isfile('QCandFT/deblur-stats/index.html'):
        tab_list.append('<li id="deblur-stats"><a href="#">Deblur QC Stats</a></li>')
        src_list.append("'deblur-stats': './QCandFT/deblur-stats/index.html'")
    if os.path.isfile('QCandFT/table-summary/index.html'):
        tab_list.append('<li id="table-summary"><a href="#">FeatureTable Summary</a></li>')
        src_list.append("'table-summary': './QCandFT/table-summary/index.html'")
    if os.path.isfile('QCandFT/rep-seq-table/index.html'):
        tab_list.append('<li id="rep-seq-table"><a href="#">Rep-Seq Table</a></li>')
        src_list.append("'rep-seq-table': './QCandFT/rep-seq-table/index.html'")
    if os.path.isfile('DiversityAnalysis/alpha-rarefaction/index.html'):
        tab_list.append('<li id="alpha-rarefaction"><a href="#">Alpha Rarefaction</a></li>')
        src_list.append("'alpha-rarefaction': './DiversityAnalysis/alpha-rarefaction/index.html'")
    if os.path.isfile('DiversityAnalysis/alpha-table.html'):
        tab_list.append('<li id="alpha-boxplots"><a href="#">Alpha-Diversity Boxplots</a></li>')
        src_list.append("'alpha-boxplots': './DiversityAnalysis/alpha-table.html'")
    if os.path.isfile('DiversityAnalysis/table.html'):
        tab_list.append('<li id="PCoA-plots"><a href="#">PCoA plots</a></li>')
        src_list.append("'PCoA-plots': './DiversityAnalysis/table.html'")
    if os.path.isfile('DiversityAnalysis/beta-table.html'):
        tab_list.append('<li id="beta-boxplot"><a href="#">Beta-Diversity Analysis</a></li>')
        src_list.append("'beta-boxplot': './DiversityAnalysis/beta-table.html'")
    if os.path.isfile('TaxonomyAnalysis/Table/index.html'):
        tab_list.append('<li id="rep-seq-taxonomy"><a href="#">OTU Table summary</a></li>')
        src_list.append("'rep-seq-taxonomy': './TaxonomyAnalysis/Table/index.html'")
    if os.path.isfile('TaxonomyAnalysis/barplots/index.html'):
        tab_list.append('<li id="taxonomy-barplot"><a href="#">Taxonomy BarPlot</a></li>')
        src_list.append("'taxonomy-barplot': './TaxonomyAnalysis/barplots/index.html'")
    if os.path.isfile('TaxonomyAnalysis/ancom.html'):
        tab_list.append('<li id="taxonomy-ancom"><a href="#">Taxonomy Analysis</a></li>')
        src_list.append("'taxonomy-ancom': './TaxonomyAnalysis/ancom.html'")

      
    html="""<!DOCTYPE html>
<html>
  <head>
    <meta charset='utf-8'>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <title>Qiime2 report</title>
    <link rel="stylesheet" href="./q2templateassets/css/bootstrap.min.css"/>
    <link rel="stylesheet" href="./q2templateassets/css/normalize.css"/>
    <link rel="stylesheet" href="./q2templateassets/css/tab-parent.css">
    <script src="./q2templateassets/js/jquery-3.2.0.min.js" charset="utf-8"></script>
    <script src="./q2templateassets/js/parent.js"></script>
  </head>
  <body>
    <div class="container-fluid">
      <h1>%s</h1>
      <div class="row">
        <ul class="nav nav-tabs">
          %s
        </ul>
      </div>
      <iframe id="tab-frame" src="./demux/overview.html"></iframe>
    </div>
    <script type="text/javascript">
      //if (window.frameElement) {
      //if (window.self != window.top) {
      //    document.getElementById('q2templatesheader').remove();
      //}
      var frameSrc = {
          %s
      };
    </script>
  </body>
</html>""" % (argvs.title,"\n".join(tab_list),",\n".join(src_list))

    with open( 'index.html', "w") as f:
        f.write(html)
        f.close()


if __name__ == '__main__':
    argvs    = setup_argparse()
    begin_t  = time.time()
    
    if not dependency_check("qiime"):
	    sys.exit("[ERROR] Executable qiime not found.")

    mkdir_p(argvs.outdir)
    abs_output = os.path.abspath(argvs.outdir)
    os.environ["TMPDIR"] =  abs_output

    print("Pipeline v%s Start: [%s]" % (__version__, datetime.datetime.now().strftime('%Y %b %d %H:%M:%S')) )

    # arguments info
    print("Arguments:")
    print("    Paired-End Reads       : %s" % argvs.paired)
    print("    Single-End Reads       : %s" % argvs.single)
    print("    Barcode Fastq File     : %s" % argvs.barcode)
    print("    MetaData Mapping File  : %s" % argvs.mappingFile)
    print("    Output Path            : %s" % abs_output)
    print("    Quality Control Method : %s" % argvs.qcMethod)
    print("    Target Amplicon        : %s" % argvs.target)
    #print("    Minimum OTU size      : %s" % argvs.minOTUsize)
    # ME
    if argvs.qcMethod.lower() == 'dada2':
        print("    Num bases trim 5'end R1: %s" % argvs.trimLeftForward)
        print("    Num bases trim 5'end R2: %s" % argvs.trimLeftReverse)
        print("    Truncate R1 Seq At Pos : %s" % argvs.truncLenForward)
        print("    Truncate R2 Seq At Pos : %s" % argvs.truncLenReverse)
        print("    Num bases trim         : %s" % (argvs.trimLen if argvs.trimLen > 0 and not argvs.truncLenForward  else 'False'))    
    
    if argvs.qcMethod.lower() == 'deblur' or argvs.qcMethod.lower() == 'otus':
        print("    Minimum Quality:       : %s" % argvs.minQuality)
        print("    Minimum Len Fraction   : %s" % argvs.minLengthFraction)
        print("    Maximum # Ambiguous    : %s" % argvs.maxAmbiguous)
    if argvs.qcMethod.lower() == 'otus':
        print("    Percent Identity       : %s" % argvs.percIdentity)

    print("    Truncate Seq At Pos    : %s" % (argvs.truncLen if argvs.truncLen > 0 and not argvs.truncLenForward else 'False'))
    print("    Minimum Sampling Depth : %s" % ('Auto' if argvs.autoDepth else argvs.samplingDepth ))
    if argvs.autoDepth:
        print("    Max Rarefaction Depth  : %s" % ('Auto'))
    else:
        print("    Max Rarefaction Depth  : %s" % ( argvs.maxRarefactionDepth if argvs.maxRarefactionDepth > 0 else argvs.samplingDepth ) )
    print("    Auto Sampling Depth    : %s" % argvs.autoDepth)
    print("    Phred Score Offset     : %d" % argvs.phred_offset)
    print("    Project Title          : %s" % (argvs.title if argvs.title else "N/A"))
    print("    CPU Number             : %d" % argvs.cpus)
    print("    Qiime2 Path            : %s" % dependency_check("qiime"))
    print("    Qiime2 Version         : %s" % tool_version_check("qiime --version | grep q2"))

   # Convert excel format to tsv 
    mappingFile = os.path.abspath(argvs.mappingFile)
    if mappingFile.lower().endswith(('.xlsx','.xls')):
        if not dependency_check("xlsx2csv"):
            shutil.rmtree(abs_output,ignore_errors=True) 
            sys.exit("[ERROR] Executable xlsx2csv not found.")
        process_cmd("xlsx2csv -d tab %s %s/tmp_sample_metadata.txt" % (mappingFile,abs_output)) 
    else:
        shutil.copy(mappingFile, abs_output+'/tmp_sample_metadata.txt') 

    # Parse Mapping file
    mappingFile, num_sample, category_list = fix_mappingFile(abs_output+'/tmp_sample_metadata.txt',abs_output)
    print("    Number of Samples      : %d" % num_sample)
    import_type = 'EMPSingleEndSequences'
    input_path = abs_output + '/input'
    demux_type = 'emp-single'
    if os.path.exists(input_path):
        for f in glob.glob("*.gz"):
            os.remove(f)
        #shutil.rmtree(input_path)  
    
    mkdir_p(input_path)

    if argvs.barcode:      
        if argvs.barcode[0].lower().endswith('.gz'):
            symlink_force(os.path.abspath(argvs.barcode[0]), input_path + '/barcodes.fastq.gz')
        else:
            process_cmd('gzip -c ' + argvs.barcode[0] + ' > ' + input_path + '/barcodes.fastq.gz')
        if argvs.single:
            read_type='se'
            if argvs.single[0].lower().endswith('.gz'):
                symlink_force(os.path.abspath(argvs.single[0]), input_path + '/sequences.fastq.gz')
            else:
                process_cmd('gzip -c ' + argvs.single[0] + ' > ' + input_path + '/sequences.fastq.gz')
        elif argvs.paired:
            read_type='pe'
            demux_type = 'emp-paired'
            import_type = 'EMPPairedEndSequences'
            if argvs.paired[0].lower().endswith('.gz'):
                symlink_force(os.path.abspath(argvs.paired[0]), input_path + '/forward.fastq.gz')
                symlink_force(os.path.abspath(argvs.paired[1]), input_path + '/reverse.fastq.gz')
            else:
                process_cmd('gzip -c ' + argvs.paired[0] + ' > ' + input_path + '/forward.fastq.gz')
                process_cmd('gzip -c ' + argvs.paired[1] + ' > ' + input_path + '/reverse.fastq.gz')
        else:
            shutil.rmtree(abs_output,ignore_errors=True) 
            sys.exit( "ERROR: Expected barcode with emp single end or paired end data" )

        import_cmd = ("qiime tools import --type  %s --input-path %s --output-path %s/input.qza" )  % (import_type, input_path, input_path)
    
    if (argvs.single or argvs.paired) and (not argvs.barcode):
        shutil.rmtree(abs_output,ignore_errors=True) 
        sys.exit( "ERROR: Expected barcode with emp single end or paired end data" )

    if argvs.dir:
        # read metadata and generteat files fof and phred scroe
        read_type, mainfest = get_fastq_manifest(mappingFile,argvs.dir,abs_output)
        if (read_type == 'pe'):
            input_format = "PairedEndFastqManifestPhred%s" % argvs.phred_offset
            import_type = 'SampleData[PairedEndSequencesWithQuality]'
        if (read_type == 'se'):
            input_format = "SingleEndFastqManifestPhred%s" % argvs.phred_offset
            import_type = 'SampleData[SequencesWithQuality]'
        
        import_cmd = ("qiime tools import --type  %s --input-path %s --output-path %s/input.qza --input-format %s" )  % (import_type, mainfest, input_path, input_format)

    
    os.chdir(abs_output)
    

    # Import Reads Data
    if not os.path.isfile('input/input.qza'):
        process_cmd(import_cmd,"Importing reads")
    
    # smaple metadata to index.html
    if not os.path.isfile('input/index.html'):
    	mappingfileTable_cmd=('qiime metadata tabulate --m-input-file %s --o-visualization input/sample-metadata.qzv') % (mappingFile)
    	process_cmd(mappingfileTable_cmd,"Sample metadata summary")
    	qiime_export_html('input/sample-metadata.qzv','input')
    ## Demultiplex
    if not os.path.isfile('demux/index.html'):
        mkdir_p('demux')
        demux_cmd = ("qiime demux %s --i-seqs %s/input.qza --m-barcodes-file %s --m-barcodes-column BarcodeSequence --o-per-sample-sequences demux/demux.qza --o-error-correction-details demux/demux-details.qza") % (demux_type,input_path,mappingFile)
        if argvs.barcode:
            process_cmd(demux_cmd,"Demultiplexing")
        else:
            symlink_force(input_path+"/input.qza","demux/demux.qza")
        
        demux_sum_cmd = ('qiime demux summarize --i-data demux/demux.qza --o-visualization demux/demux.qzv')
        process_cmd(demux_sum_cmd,"Demultiplexing summary")

        qiime_export_html('demux/demux.qzv','demux')

    # Quality control and Feature Table contruction
    if not os.path.isfile('QCandFT/rep-seq-table/index.html'):
        mkdir_p('QCandFT')
        chimera_method = "consensus"       
        
        # ME added the 'single/paired' switch
        if argvs.qcMethod.lower() == 'dada2':
            if argvs.single or read_type == 'se':
                dada2_cmd = ("qiime dada2 denoise-single --i-demultiplexed-seqs demux/demux.qza "
                            "--p-chimera-method %s "
                            "--o-representative-sequences QCandFT/rep-seqs.qza "
                            "--o-table QCandFT/table.qza "
                            "--o-denoising-stats QCandFT/stats-dada2.qza "
                            "--p-trunc-len %d "
                            "--p-trim-left %d "
                            "--p-n-threads %d "  ) % (chimera_method, argvs.truncLen, argvs.trimLen, argvs.cpus)
            
            elif argvs.paired or read_type == 'pe':

                dada2_cmd = ("qiime dada2 denoise-paired --i-demultiplexed-seqs demux/demux.qza "
                            "--p-chimera-method %s "
                            "--o-representative-sequences QCandFT/rep-seqs.qza "
                            "--o-table QCandFT/table.qza "
                            "--o-denoising-stats QCandFT/stats-dada2.qza "
                            "--p-trim-left-f %d "
                            "--p-trim-left-r %d "
                            "--p-trunc-len-f %d "
                            "--p-trunc-len-r %d "
                            "--p-n-threads %d "  ) % (chimera_method, argvs.trimLeftForward, argvs.trimLeftReverse, argvs.truncLenForward, argvs.truncLenReverse, argvs.cpus)
                            
            process_cmd(dada2_cmd,"Dada2 QC and FeatureTable construction")
            stats_cmd = ('qiime metadata tabulate --m-input-file QCandFT/stats-dada2.qza --o-visualization QCandFT/stats-dada2.qzv')
            process_cmd(stats_cmd, 'Tabulate QC Stats')
            qiime_export_html('QCandFT/stats-dada2.qzv','QCandFT/filter-stats')
        elif argvs.qcMethod.lower() == 'deblur' or  argvs.qcMethod.lower() == 'otus':
            file_for_check_truncate_len = 'demux/quality-plot.html'
            if argvs.paired or read_type == 'pe':
                join_cmd = ('qiime vsearch join-pairs --i-demultiplexed-seqs demux/demux.qza --o-joined-sequences demux/demux-joined.qza')
                if not os.path.isfile('demux/demux-joined.qza'):
                    process_cmd(join_cmd,"Vsearch join-pairs")
                
                demux_sum_cmd = ('qiime demux summarize --i-data demux/demux-joined.qza --o-visualization demux/demux-joined.qzv')
                if not os.path.isfile('demux-joined/index.html'):
                    process_cmd(demux_sum_cmd,"Vsearch join-pairs summary")
                    qiime_export_html('demux/demux-joined.qzv','demux-joined')

                file_for_check_truncate_len = 'demux-joined/quality-plot.html'

                qf_cmd = ('qiime quality-filter q-score --i-demux demux/demux-joined.qza --o-filtered-sequences QCandFT/filtered.qza '
                        '--p-min-quality %s '
                        '--p-min-length-fraction %s '
                        '--p-max-ambiguous %s '
                        '--o-filter-stats QCandFT/filter-stats.qza') % (argvs.minQuality, argvs.minLengthFraction, argvs.maxAmbiguous)
            elif argvs.single or read_type == 'se':
                qf_cmd = ('qiime quality-filter q-score --i-demux demux/demux.qza --o-filtered-sequences QCandFT/filtered.qza '
                        '--p-min-quality %s '
                        '--p-min-length-fraction %s '
                        '--p-max-ambiguous %s '
                        '--o-filter-stats QCandFT/filter-stats.qza') % (argvs.minQuality, argvs.minLengthFraction, argvs.maxAmbiguous)
            
            if not os.path.isfile('QCandFT/filtered.qza'):
                process_cmd(qf_cmd,"Quality-filter q-score")

            if not os.path.isfile('QCandFT/filter-stats/index.html'):
                stats_cmd = ('qiime metadata tabulate --m-input-file QCandFT/filter-stats.qza --o-visualization QCandFT/filter-stats.qzv')
                process_cmd(stats_cmd, 'Tabulate QC Stats')
                qiime_export_html('QCandFT/filter-stats.qzv','QCandFT/filter-stats')
            
            # Deblur
            # need use denoise-other for other db and provide reference sequence
            # data/sh_refs_qiime_ver8_dynamic_02.02.2019.fasta
            if argvs.qcMethod.lower() == 'deblur':
                if (argvs.target.lower() == 'its'):
                    reference_seqs_opts = '--i-reference-seqs ' + target_path + '/ITS.qza'
                    deblur_method = 'denoise-other'
                else:
                    reference_seqs_opts = ""
                    deblur_method = 'denoise-16S'
           
                trim_length = argvs.truncLen if argvs.truncLen > 0 else auto_determine_truncate_len_deblur(file_for_check_truncate_len)
           
                deblur_cmd = ("qiime deblur %s %s --i-demultiplexed-seqs QCandFT/filtered.qza "
                          "--p-trim-length %d --o-representative-sequences QCandFT/rep-seqs.qza "
                          "--o-table QCandFT/table.qza --p-sample-stats --o-stats QCandFT/deblur-stats.qza "
                          "--p-jobs-to-start %d " )  % (deblur_method,reference_seqs_opts,trim_length,argvs.cpus)
                process_cmd(deblur_cmd,"Deblur QC and FeatureTable construction")
                deblur_vis_cmd = ('qiime deblur visualize-stats --i-deblur-stats QCandFT/deblur-stats.qza --o-visualization QCandFT/deblur-stats.qzv')
                process_cmd(deblur_vis_cmd, 'Deblur vis-Stats')
                qiime_export_html('QCandFT/deblur-stats.qzv','QCandFT/deblur-stats')
            # OTU clustering
            elif argvs.qcMethod.lower() == 'otus':
                reference_seqs = target_path + '/gg_2022_10.qza'
                if (argvs.target.lower() == 'silva'):
                    reference_seqs = target_path + '/silva-138-99.qza'
                elif (argvs.target.lower() == 'silva-v3-v4'):
                    reference_seqs = target_path + '/silva_138_99PercClust_16SOnly_unaligned_F341Mod_R806Mod_primerMatchPortionOnly_forTaxClassTraining.qza'
                elif (argvs.target.lower() == 'its'):
                    reference_seqs = target_path + '/ITS.qza'

                if not os.path.isfile('QCandFT/filtered-dedup-seq.qza'):
                    dedup_cmd = ("qiime vsearch dereplicate-sequences --i-sequences QCandFT/filtered.qza --o-dereplicated-table QCandFT/filtered-dedup-table.qza "
                                "--o-dereplicated-sequences QCandFT/filtered-dedup-seq.qza")
                    process_cmd(dedup_cmd,"Dereplicate sequences")
        
                if not os.path.isfile('QCandFT/clustered-seqs.qza'):
                    clustering_cmd = ("qiime vsearch cluster-features-open-reference --i-sequences QCandFT/filtered-dedup-seq.qza "
                                "--i-table QCandFT/filtered-dedup-table.qza "
                                "--o-clustered-table QCandFT/clustered-table.qza "
                                "--o-clustered-sequences QCandFT/clustered-seqs.qza "
                                "--o-new-reference-sequences QCandFT/new-ref-seqs.qza "
                                "--i-reference-sequences %s "
                                "--p-perc-identity %f "
                                "--p-threads %d ") % (reference_seqs, argvs.percIdentity, argvs.cpus)
                    process_cmd(clustering_cmd,"Open Reference OTU clustering")

                if not os.path.isfile('QCandFT/non-chimeras.qza'):
                    chimera_cmd = ("qiime vsearch uchime-ref --i-table QCandFT/clustered-table.qza --i-sequences QCandFT/clustered-seqs.qza "
                                "--o-chimeras QCandFT/chimeras.qza "
                                "--o-nonchimeras QCandFT/non-chimeras.qza "
                                "--o-stats QCandFT/chimera-stats.qza "
                                "--i-reference-sequences %s "
                                "--p-threads %d ") % (reference_seqs, argvs.cpus)
                    process_cmd(chimera_cmd, "Chimera Detection")

                if not os.path.isfile('QCandFT/table.qza'):
                    chimera_rmt_cmd = ("qiime feature-table filter-features --i-table QCandFT/clustered-table.qza --m-metadata-file QCandFT/non-chimeras.qza --o-filtered-table QCandFT/table.qza")
                    process_cmd(chimera_rmt_cmd, "Chimera Remove from Features Table")
                if not os.path.isfile('QCandFT/rep-seqs.qza'):
                    chimera_rms_cmd = ("qiime feature-table filter-seqs --i-table QCandFT/table.qza --i-data QCandFT/clustered-seqs.qza --o-filtered-data QCandFT/rep-seqs.qza")
                    process_cmd(chimera_rms_cmd, "Chimera Sequences Remove")

                if not os.path.isfile('QCandFT/chimera-stats/index.html'):
                    stats_cmd = ('qiime metadata tabulate --m-input-file QCandFT/chimera-stats.qza --o-visualization QCandFT/chimera-stats.qzv')
                    process_cmd(stats_cmd, 'Tabulate Chimera Detection Stats')
                    qiime_export_html('QCandFT/chimera-stats.qzv','QCandFT/chimera-stats')

        else:
            sys.exit( "ERROR: Expected QC methods are dada2 and deblur. Your input is %s" ) % (argvs.qcMethod)

        ft_sum_cmd = ('qiime feature-table summarize --i-table QCandFT/table.qza --o-visualization QCandFT/table.qzv'
                      ' --m-sample-metadata-file %s ')  % (mappingFile)
        process_cmd(ft_sum_cmd, 'Tabulate FeatureTable')
        qiime_export_html('QCandFT/table.qzv','QCandFT/table-summary')
        qiime_export_html('QCandFT/table.qza','QCandFT/table-summary')
        biom_otu_cmd=('biom convert -i %s -o %s --to-tsv') % ("QCandFT/table-summary/feature-table.biom","QCandFT/table-summary/feature-table.tsv")
        process_cmd(biom_otu_cmd, 'Generate Biom OTU FeatureTable')
        ft_tseq_cmd = ('qiime feature-table tabulate-seqs --i-data QCandFT/rep-seqs.qza --o-visualization QCandFT/rep-seqs.qzv')
        process_cmd(ft_tseq_cmd, 'Tabulate rep seqs')
        qiime_export_html('QCandFT/rep-seqs.qzv','QCandFT/rep-seq-table')
    
    # Automatically adjust the sampling depth
    samplingDepth = argvs.samplingDepth
    median_depth = argvs.samplingDepth
    num_sample_over_samplingDepth = num_sample
    depth_list=list()
    sample_freq_file = "QCandFT/table-summary/sample-frequency-detail.csv"
    if os.path.isfile(sample_freq_file):
        autoSamplingDepth, median_depth, depth_list = auto_determine_sampling_depth(sample_freq_file,samplingDepth)

        if argvs.autoDepth and autoSamplingDepth:
            samplingDepth = autoSamplingDepth

    num_sample_over_samplingDepth=sum(1 for i in depth_list if i >= samplingDepth)
    #get_depth_cutoff_from_biom_summary_table

    # Phylogenetic diversity analyses
    if not os.path.isfile('PhyloAnalysis/rooted-tree/tree.nwk'):
        mkdir_p('PhyloAnalysis')
        phy_cmd = ("qiime phylogeny align-to-tree-mafft-fasttree --i-sequences QCandFT/rep-seqs.qza --o-alignment PhyloAnalysis/aligned-rep-seqs.qza "
                   "--o-masked-alignment PhyloAnalysis/masked-aligned-rep-seqs.qza --o-tree PhyloAnalysis/unrooted-tree.qza "
                   "--o-rooted-tree PhyloAnalysis/rooted-tree.qza "
                   "--p-n-threads %d " ) % (argvs.cpus)
        process_cmd(phy_cmd, 'Phylogenetic diversity analyses')
        qiime_export_html('PhyloAnalysis/rooted-tree.qza','PhyloAnalysis/rooted-tree')

    # Alpha and beta diversity analysis
    
    if not os.path.isfile('DiversityAnalysis/unweighted_unifrac_emperor.qzv') and num_sample_over_samplingDepth > 2 :
        if os.path.exists('DiversityAnalysis'):
            shutil.rmtree('DiversityAnalysis',ignore_errors=True)
        diversity_cmd = ("qiime diversity core-metrics-phylogenetic --i-phylogeny PhyloAnalysis/rooted-tree.qza --i-table QCandFT/table.qza "
                         "--p-sampling-depth %d --m-metadata-file %s --output-dir DiversityAnalysis " 
                         "--p-n-jobs-or-threads %d " ) % (samplingDepth,mappingFile, int(argvs.cpus/2) if argvs.cpus > 1 else 1 )
        process_cmd(diversity_cmd, 'Alpha and Beta diversity analyses')
        beta_analysis_results = ['unweighted_unifrac','jaccard','bray_curtis','weighted_unifrac']
        for x in beta_analysis_results:
            if os.path.isfile("DiversityAnalysis/%s_emperor.qzv" % x ):
                qiime_export_html("DiversityAnalysis/%s_emperor.qzv" % x ,"DiversityAnalysis/%s_emperor" % x)
            if os.path.isfile("DiversityAnalysis/%s_distance_matrix.qza" % x ):
                qiime_export_html("DiversityAnalysis/%s_distance_matrix.qza" % x ,"DiversityAnalysis/%s_distance_matrix" % x)
                if len(category_list) > 0:
            	    for g in category_list:
                        ## disable --p-pairwise to avoid error in https://forum.qiime2.org/t/error-on-beta-group-significance/1789/6 ??
                        cmd = ("qiime diversity beta-group-significance --i-distance-matrix DiversityAnalysis/%s_distance_matrix.qza "
                                "--m-metadata-file %s --m-metadata-column '%s' --o-visualization DiversityAnalysis/%s-%s-significance.qzv --p-pairwise") % (x, mappingFile, g, x, g)
                        process_cmd(cmd, '%s: Beta diversity PERMANOVA test' % (g) )
                        if os.path.isfile("DiversityAnalysis/%s-%s-significance.qzv" % (x, g)):
                            qiime_export_html("DiversityAnalysis/%s-%s-significance.qzv" % (x,g),"DiversityAnalysis/%s_%s_boxplot" % (x,g))

        add_html_table('beta')
   
        alpha_analysis_result = ['evenness_vector','faith_pd_vector','observed_otus_vector','shannon_vector']
        for x in alpha_analysis_result:
            if os.path.isfile("DiversityAnalysis/%s.qza" % x) :
                cmd = ("qiime diversity alpha-group-significance --i-alpha-diversity DiversityAnalysis/%s.qza "
                       "--m-metadata-file %s --o-visualization DiversityAnalysis/%s-group-significance.qzv ") % (x, mappingFile,x)
                process_cmd(cmd)

                if os.path.isfile("DiversityAnalysis/%s-group-significance.qzv" % x) :
                    qiime_export_html("DiversityAnalysis/%s-group-significance.qzv" % x,"DiversityAnalysis/%s_boxplot" % x)

        add_html_table('alpha')

    add_html_table('beta-significance',category_list)

    
   
    # Alpha rarefaction plotting
    max_rarefaction_depth = argvs.maxRarefactionDepth if argvs.maxRarefactionDepth > 0 else samplingDepth
    #max_rarefaction_depth = argvs.maxRarefactionDepth if argvs.maxRarefactionDepth > 0 else median_depth
    if not os.path.isfile('DiversityAnalysis/alpha-rarefaction.qzv'):
        mkdir_p('DiversityAnalysis')
        diversity_cmd = ("qiime diversity alpha-rarefaction --i-table QCandFT/table.qza --i-phylogeny PhyloAnalysis/rooted-tree.qza "
                         "--p-max-depth %d --m-metadata-file %s "
                         "--o-visualization DiversityAnalysis/alpha-rarefaction.qzv")  % (max_rarefaction_depth,mappingFile)
        process_cmd(diversity_cmd, 'Alpha rarefaction analyses')
        qiime_export_html('DiversityAnalysis/alpha-rarefaction.qzv','DiversityAnalysis/alpha-rarefaction')

    # Taxonomic analysis
    nb_classifier = target_path + '/gg_2022_10_backbone_full_length.nb.qza'
    if (argvs.target.lower() == 'silva'):
        nb_classifier = target_path + '/silva-138-99-nb-classifier.qza'
    elif (argvs.target.lower() == 'silva-v3-v4'):
        nb_classifier = target_path + '/silva_138_99PercClust_16SOnly_unaligned_F341Mod_R806Mod_primerMatchPortionOnly_7LevelTaxonomy_naiveBayesClassifier.qza'
    elif (argvs.target.lower() == 'its'):
        nb_classifier = target_path + '/ITS.classifier.qza'

    #  QCandFT/table.qza is the full data 
    #  use rarefied_table for sub-sampling taxanomy plots  
    raw_table = 'QCandFT/table.qza'
    data_table = 'DiversityAnalysis/rarefied_table.qza' if os.path.isfile('DiversityAnalysis/rarefied_table.qza') else 'QCandFT/table.qza'
    if not os.path.isfile('TaxonomyAnalysis/taxa-bar-plots.qzv'):
        mkdir_p('TaxonomyAnalysis')
        classifier_cpu = 8 if argvs.cpus >8 else argvs.cpus
        taxa_cmd = ("qiime feature-classifier classify-sklearn --i-classifier %s " 
                    "--i-reads QCandFT/rep-seqs.qza --o-classification TaxonomyAnalysis/taxonomy.qza "
                    "--p-n-jobs %d " ) % (nb_classifier,classifier_cpu) 
        process_cmd(taxa_cmd, 'Taxonomic analyses')
        taxa_table_cmd = ('qiime metadata tabulate --m-input-file QCandFT/rep-seqs.qza --m-input-file TaxonomyAnalysis/taxonomy.qza --o-visualization TaxonomyAnalysis/taxonomy.qzv')
        process_cmd(taxa_table_cmd, 'Taxonomic tabulate')
        qiime_export_html('TaxonomyAnalysis/taxonomy.qzv','TaxonomyAnalysis/Table')

        taxa_barplot_cmd = ("qiime taxa barplot --i-table %s --i-taxonomy TaxonomyAnalysis/taxonomy.qza "
                            "--m-metadata-file %s " 
                            "--o-visualization TaxonomyAnalysis/taxa-bar-plots.qzv") % (data_table, mappingFile)
        process_cmd(taxa_barplot_cmd, 'Taxonomic barplots')
        qiime_export_html('TaxonomyAnalysis/taxa-bar-plots.qzv','TaxonomyAnalysis/barplots')

    ## Taxonomic ANCOM analysis with category metadata field
    if (len(category_list)>0 and not os.path.isfile('TaxonomyAnalysis/ancom.html')):
        cmd = ("qiime taxa collapse --i-table %s --i-taxonomy TaxonomyAnalysis/taxonomy.qza "
              " --p-level 6 --o-collapsed-table TaxonomyAnalysis/table-l6.qza") % (raw_table)
        process_cmd(cmd)
        cmd = ("qiime composition add-pseudocount --i-table TaxonomyAnalysis/table-l6.qza --o-composition-table TaxonomyAnalysis/comp-table-l6.qza");
        process_cmd(cmd)
        for g in category_list:
            cmd = ("qiime composition ancom --i-table TaxonomyAnalysis/comp-table-l6.qza --m-metadata-file %s "
                   "--m-metadata-column '%s' --o-visualization TaxonomyAnalysis/l6-ancom-%s.qzv ") % (mappingFile, g, g)
            process_cmd(cmd,'%s: Taxonomic abundance testing with ANCOM' % (g))
            qiime_export_html('TaxonomyAnalysis//l6-ancom-%s.qzv' % (g),'TaxonomyAnalysis/l6-ancom-%s' % (g) )

        add_html_table('taxonomy-ancom',category_list)
    
    ## Creating a OTU table with taxonomy annotations use rarefied_table
    if not os.path.isfile('DiversityAnalysis/rarefied-table-summary/feature-table.tsv'):
        qiime_export_html(data_table,'DiversityAnalysis/rarefied-table-summary')  
        biom_otu_cmd=('biom convert -i %s -o %s --to-tsv') % ("DiversityAnalysis/rarefied-table-summary/feature-table.biom","DiversityAnalysis/rarefied-table-summary/feature-table.tsv")
        process_cmd(biom_otu_cmd, 'Generate Biom OTU Rarefied FeatureTable')
    if os.path.isfile('DiversityAnalysis/rarefied-table-summary/feature-table.tsv') and os.path.isfile('TaxonomyAnalysis/Table/metadata.tsv'):
        ft=pd.read_csv('DiversityAnalysis/rarefied-table-summary/feature-table.tsv',delimiter='\t',encoding='utf-8',header=1)
        tax=pd.read_csv('TaxonomyAnalysis/Table/metadata.tsv',delimiter='\t',encoding='utf-8',skiprows=[1])
        if not os.path.isfile('TaxonomyAnalysis/Table/feature-table-taxanomy.tsv'):
            ft.merge(tax,left_on='#OTU ID',right_on='id').drop(['id'],axis=1).to_csv("TaxonomyAnalysis/Table/feature-table-taxanomy.tsv",sep="\t",index=False)
        if not os.path.isfile('TaxonomyAnalysis/Table/feature-table-taxanomy.qzv'):
            mappingfileTable_cmd=('qiime metadata tabulate --m-input-file TaxonomyAnalysis/Table/feature-table-taxanomy.tsv --o-visualization TaxonomyAnalysis/Table/feature-table-taxanomy.qzv')
            process_cmd(mappingfileTable_cmd,"OTUs like table summary")
            qiime_export_html('TaxonomyAnalysis/Table/feature-table-taxanomy.qzv','TaxonomyAnalysis/Table')
    
    # HTML REPORT
    html_report(script_path+'/q2templateassets')
    
    # zip all report html
    if argvs.zip:
        os.chdir("..")
        output_dirname = os.path.basename(abs_output)
        if argvs.title:
            zip_file = 'qiime2_' + argvs.title + '.zip'
        else:
            zip_file = 'qiime2_' + output_dirname + '.zip'
        zip_cmd = ("zip --symlinks -r %s %s -x *qza ") % (zip_file,output_dirname)
        process_cmd(zip_cmd, 'Compressing output')
        os.rename(zip_file,abs_output+'/'+zip_file)

    print("\nTotal %s" % get_runtime(begin_t))
    
    
    
