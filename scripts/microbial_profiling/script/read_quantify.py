#!/usr/bin/env python3
''''''''''''''''''''''''''''' 
.'.         .'.                                                  
|  \       /  |                                          
'.  \  |  /  .'                                       
  '. \\|// .'      Hey!                                              
    '-- --'       Listen!      
    .'/|\'.            
   '..'|'..'        
'''''''''''''''''''''''''''''
 ###################################################################
#       This script is currently being reworked.                    #
#       It is now able to take one input fastq file for analysis    #
#          with multiple marker files. These marker files need      #
#          to be named by the number of target sequences that were  #
#          hit during inclusivity testing. That's all this does     #
#          for now. Well, it prints out the table to a file, too.   #
#                                                                   #
#       It will loop over all marker sequences fed to it and expects#
#           that these have been somewhat curated already.          #
#       If you want to analyze multiple files, loop this over them. #
#                                                                   #
#       ShortBRED quantify is at the core of this script and expects#
#           fasta files. A sed command is used to convert fastqs to #
#           fasta format.                                           #
#           It does require fastq format to start.                  #
 ###################################################################

 #####################################
#                                     #
#    It's dangerous to go alone!      #
#           Take these!               #
#                                     #
 #####################################
import matplotlib
matplotlib.use('Agg')
import shutil
import argparse as ap, pandas as pd, sys, os, fnmatch, glob, numpy as np, re
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib.pyplot import xticks

version ='beta'

def parse_params(ver): #Are there too many parameters?
    class SmartFormatter(ap.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()
            return ap.HelpFormatter._split_lines(self, text, width)

    p = ap.ArgumentParser(prog='quantify_and_tabulate', description="""Script for reporting pathogen detection by input target pathogen %s""" % ver, formatter_class=SmartFormatter)

    p.add_argument('-fq', '--FastQ',
        metavar='[FQ]', type=str,
            help="""Complete path to the input fastq. """)

    p.add_argument('-mp', '--MarkerPath',
        metavar='[MP]', type=str,
            help="""The complete path to SB Marker Files. """)

    p.add_argument('-rp', '--ResultsPath',
        metavar='[RP]', type=str,
            help="""The complete path to Results Files. """)

    p.add_argument('-t', '--Threads',
        metavar='[THREADS]', type=str,
            help="""Number of threads for SB to run wtih. """)

    p.add_argument('-sp', '--SearchProgram', default ="rapsearch2", 
        metavar='[SP]', type=str,
            help="""Choose program for serach. Default is \"rapsearch2\". Support: rapsearch2, diamond and usearch in PATH""")

    args_parsed = p.parse_args()

    """ Check the options """

    return args_parsed

#Use these to human-sort marker files that begin with an integer
def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    l.sort(key=alphanum_key)
#done

def quantify(marker_path, fastq, result_dir, threads, search_program):


    fq = os.path.basename(fastq) #get the filename
    fq = fq.split('.fastq')[0] #remove the extension

    reads_path = '{0}/reads'.format(result_dir)
    if not os.path.exists(reads_path):
        os.makedirs(reads_path)

    fasta = '{0}/{1}.fna'.format(reads_path,fq)
    fasta

    os.system('sed -n \'1~4s/^@/>/p;2~4p\' {0} > {1}'.format(fastq,fasta))
#            SeqIO.convert(combined_fastq,'fastq',combined_fasta,'fasta')
            #ShortBRED Quantify starts here. Note this is only invoked if the parent taxID is found. Hence it still being in the if structure.
    markerfiles = glob.glob('{0}/*markers.txt'.format(marker_path))

    sort_nicely(markerfiles) #sort them from least specific to most specific

    for marker in markerfiles:
        prefix = os.path.basename(marker)
        prefix = prefix.split('.')[0] #for naming results
        print(prefix)
        tmp = '{0}/testing.{1}.tmp'.format(result_dir, prefix) #SB quantify parameter for temp files
        family_results = '{0}/testing.{1}.family.results.txt'.format(result_dir,prefix) #SB quantify parameter for marker-family-level results
        marker_results = '{0}/testing.{1}.marker.results.txt'.format(result_dir,prefix) #SB  parameter for marker-level results
        quant_cmd = ("shortbred_quantify.py --markers {0} --wgs {1} --results {2} --tmp {3} --marker_results {4} --threads {5} --search_program {6}".format(marker, fasta, family_results, tmp, marker_results, threads, search_program))
        os.system(quant_cmd) #Run SB  quantifywith parameters. Not that it's relevant here, but note that you can't do this from within a python interpreter-for loop directly. However, you can do this in a python interpreter and not in a for loop. Riddle me that.
        table = pd.read_csv('{0}/testing.{1}.marker.results.txt'.format(result_dir,prefix),sep='\t',usecols=['Hits']) #Load results into a table so we can see if there's a reason to continue
 #       if sum(table.Hits) == 0: # Assuming the marker files are iterated through from least specific to most specific markers, once a set of reads doesn't hit to any of them, the more specific ones don't need to be run. It would be redundant.
  #          break


#HEY!
#THE REST STILL NEEDS FIXING!!
def tabulate(results_dir):

    resultsfiles = glob.glob('{0}/*marker.results.txt'.format(results_dir)) #grab resultsfiles
    columns = []
    sums = []
    for resultfile in resultsfiles:
        file = pd.read_csv(resultfile,sep='\t') #read in each file
        name = os.path.basename(resultfile)
        columns.append(resultfile.split('.')[1]) #name the column with the first indicator per marker
#CHANGE THE ABOVE FROM 1 TO 0 AFTER WE STOP "TESTING" BECAUSE CURRENTLY THE WORD "TESTING" IS THE FIRST SPLIT.
        sums.append(file.Hits.sum()) #get sum of Hits column for each marker set.

    df = pd.DataFrame(columns=columns) #build dataframe
    df.loc[0] = sums #add sums to df

    sum_table = '{0}/marker.sum.table.tsv'.format(results_dir) #this is where they're going

    df.to_csv(sum_table,sep='\t') #and there it is



''' Replaceable
#Step 1: sum marker file columns - get one row for each minimum marker hit value by file
#Step 2: include the RNR values by file 
#Step 3: transpose so that rows are files and columns are marker thresholds and RNR?
#Step 4: graph it up!

#Step 0 - read in combined marker file

   #GRAPH IT UP!
    graph_handle = (tax_path + '/' + target + '.RNR_markers.png')
    results = pd.read_csv(out_handle,sep='\t',index_col=0)
    ax = results.plot.bar(x=results.index, y=['Sum','RNR'], title='Comparison of Characteristic Markers and PanGIA RNR for {0}'.format(targ_name), legend=True)
    plt.yscale('symlog')
    locs, labels = xticks()
#    xticks(np.arange(len(results.index)), results.index, rotation=90)
    plt.xlabel('Sample')
    plt.ylabel('Read Abundance')
#    plt.grid(True)
    plt.tight_layout()
    plt.savefig(graph_handle)
'''  

if __name__ == '__main__':
    argvs = parse_params(version)

    # Add EDGE ShortBRED path
    edge_home = os.environ["EDGE_HOME"]
    if os.path.exists(edge_home):
        os.environ["PATH"] = edge_home + "/bin/ShortBRED" + os.pathsep + os.environ["PATH"]

    #dependency check
    if sys.version_info < (3,0):
        sys.exit( "[ERROR] Python 3.0 or above is required.")
    if not shutil.which("shortbred_quantify.py"):
        sys.exit("[ERROR]: Cannot find EDGE ShortBred scirpts")

    result_dir = '{0}'.format(argvs.ResultsPath)

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    quantify(argvs.MarkerPath, argvs.FastQ, result_dir, argvs.Threads, argvs.SearchProgram)

    tabulate(result_dir)
#This will need to change - it shouldn't be going into each marker path and tabulating - each one is only one column. It should be glomming each combined results file together.
#    tabulate(argvs.RNRFile, argvs.TaxID, tax_path)
