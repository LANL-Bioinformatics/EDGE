#!/bin/bash
#$ -cwd
#$ -l h_vmem=10G
#$ -m abe
#$ -j y

set -e;

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -d <DATABASE> -o <OUTDIR> -p <PREFIX> -f <PLATFORM> [OPTIONS]

OPTIONS:
   -i      Input a FASTQ file (singleton or pair-end sequences separated by comma)
   -d      Database
   -o      Output directory
   -p      Output prefix
   -a      BWA alignment parameters.
   -t      Number of threads (default: 24)
   -f      platform pacbio
   -h      help
EOF
}

FASTQ=
REFDB=$EDGE_HOME/database/bwa_index/NCBI-Bacteria-Virus.fna
PREFIX=
OUTPATH=
THREADS=24
INPUT_OPTION=
PLATFORM=
PLATFORM_OPTION=
BWASCORECUT=
BWAMETHOD=

while getopts "i:d:o:p:t:f:a:b:h" OPTION
do
     case $OPTION in
        i) FASTQ=$OPTARG
           ;;
        d) REFDB=$OPTARG
           ;;
        o) OUTPATH=$OPTARG
           ;;
        p) PREFIX=$OPTARG
           ;;
        t) THREADS=$OPTARG
           ;;
        f) PLATFORM=$OPTARG
           ;;
        b) BWASCORECUT=$OPTARG
           ;;
	a) BWAMETHOD=$OPTARG
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

export PATH=$EDGE_HOME/bin:$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;

mkdir -p $OUTPATH

echo "[BEGIN]"

set -x;

bwa mem -t $THREADS $BWAMETHOD -T $BWASCORECUT $REFDB $FASTQ > $OUTPATH/$PREFIX.sam
# bwa_sam2read_taxa.pl species preload < $OUTPATH/$PREFIX.sam > $OUTPATH/$PREFIX.out.read_classification &
bwa_sam2giReadCount.pl < $OUTPATH/$PREFIX.sam > $OUTPATH/${PREFIX}.csv &

wait

convert_gi2list.pl < $OUTPATH/$PREFIX.csv > $OUTPATH/$PREFIX.out.list &
convert_gi2tabTree.pl < $OUTPATH/$PREFIX.csv > $OUTPATH/$PREFIX.out.tab_tree
ktImportText $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

wait

set +ex;
echo "";
echo "[END] $OUTPATH $PREFIX";
