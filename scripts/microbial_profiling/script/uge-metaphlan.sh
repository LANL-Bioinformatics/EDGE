#!/bin/bash
#$ -cwd
#$ -l h_vmem=2.6G
#$ -m abe
#$ -j y
set -e;

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -o <OUTDIR> -p <PREFIX> [OPTIONS]

OPTIONS:
   -i      Input a FASTQ file
   -o      Output directory
   -p      Output prefix
   -t      Number of threads
   -h      help
EOF
}

FASTQ=
PREFIX=
OUTPATH=
THREADS=24

while getopts "i:o:p:t:h" OPTION
do
     case $OPTION in
        i) FASTQ=$OPTARG
           ;;
        o) OUTPATH=$OPTARG
           ;;
        p) PREFIX=$OPTARG
           ;;
        t) THREADS=$OPTARG
           ;;
        h) usage
           exit
           ;;
     esac
done

if [[ -z "$FASTQ" || -z "$OUTPATH" || -z "$PREFIX" ]]
then
     usage;
     exit 1;
fi

export PYTHONPATH=$EDGE_HOME/bin/python/lib
export PATH=$EDGE_HOME/bin:$EDGE_HOME/scripts:$PATH;
mkdir -p $OUTPATH

set -x;
#run metaphlan
time metaphlan.py --bowtie2db $EDGE_HOME/database/metaphlan/mpa --bowtie2out ${OUTPATH}/${PREFIX}.bt2out --nproc ${THREADS} ${FASTQ} ${OUTPATH}/${PREFIX}.out.mpln

#parse mpln
convert_metaphlan2tabTree.pl < $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.tab_tree
mpln2krona.pl -l s -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.krona

ktImportText -a $OUTPATH/$PREFIX.out.krona -o $OUTPATH/$PREFIX.krona.html

mpln2list.pl  -l s -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.list
mpln2megan.pl -l s -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.megan

set +ex;
echo "";
echo "[END] $OUTPATH $PREFIX";
