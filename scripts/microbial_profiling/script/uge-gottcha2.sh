#!/bin/bash
#$ -l h_vmem=2.6G
#$ -j y
#$ -cwd

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -o <OUTDIR> -p <PREFIX> -l <LEVEL> [OPTIONS]

ARGUMENTS:
   -i      Input a FASTQ file or pair-ended FASTAQ files sperated by comma
   -o      Output directory
   -p      Output prefix
   -l      Level [genus|species|strain]
   -k      Database type - bacteria or viruses [bact|viru|both]. Default is "bact".

OPTIONS:
   -t      Number of threads. [default is 4]
   -a      BWA alignment parameters. [default is "mem -k 30 -T 0 -B 100 -O 100 -E 100"]
   -q      minQ: minimum quality of any single base [default is 20]
   -f      fixL: chunk all fragments of reads into smaller fragments of length [default is 30]
   -m      minL (*** disabled, don't use ***)
   -s      The FULL PATH and prefix of pre-splitrimmed sequences and stats file.
   -h      help
EOF
}

###########################
#
# Default values
#
###########################

FASTQ=
PREFIX=
OUTPATH=
LVL="species"
THREADS=4
PRE_SPLITRIM=
DB=$EDGE_HOME/database/GOTTCHA2/

while getopts "i:o:p:l:d:t:a:q:f:m:s:h" OPTION
do
     case $OPTION in
        i) FASTQ=$OPTARG
           ;;
        o) OUTPATH=$OPTARG
           ;;
        p) PREFIX=$OPTARG
           ;;
        l) LVL=$OPTARG
           ;;
        d) DB+=$OPTARG
           ;;
        t) THREADS=$OPTARG
           ;;
        a) BWAMETHOD=$OPTARG
           ;;
        q) TRIM_MINQ=$OPTARG
           ;;
        f) TRIM_FIXL=$OPTARG
           ;;
        m) TRIM_MINL=$OPTARG
           ;;
        s) PRE_SPLITRIM=$OPTARG
           ;;
        h) usage
           exit
           ;;
     esac
done

## path
export PATH=$EDGE_HOME/thirdParty/gottcha2:$EDGE_HOME/bin:$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;

mkdir -p $OUTPATH

set -xe;

gottcha.py -r READ_COUNT -i $FASTQ -t $THREADS --outdir $OUTPATH -p $PREFIX --database $DB

awk -F\\t '{if(NR==1){out=$1"\t"$2"\tROLLUP\tASSIGNED"; { for(i=3;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.summary.tsv > $OUTPATH/$PREFIX.out.list
awk -F\\t '{if(NR>1){out=$1"\t"$2"\t"$4"\t"; { for(i=3;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.summary.tsv >> $OUTPATH/$PREFIX.out.list

gottcha.py -r READ_COUNT -s $OUTPATH/$PREFIX.gottcha_*.sam -m lineage -c > $OUTPATH/$PREFIX.out.tab_tree

#generate KRONA chart
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

set +xe;
echo "";
echo "[END] $OUTPATH $PREFIX";
