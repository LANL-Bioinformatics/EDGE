#!/bin/bash
#$ -cwd
#$ -l h_vmem=2.6G
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

while getopts "i:d:o:p:t:f:h" OPTION
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

if [[ "$FASTQ" == *,* ]]
then
	IFS=',' read -a file <<< "$FASTQ"
	if [[ -z "${file[2]}" ]]
	then
		INPUT_OPTION="-p $FASTQ"
	else
		INPUT_OPTION="-u ${file[0]} -p '${file[1]} ${file[2]}'"
	fi
else
	INPUT_OPTION="-u $FASTQ"
fi


if [[ "$PLATFORM" == *pacbio* ]]
then 
	PLATFORM_OPTION="--pacbio"
fi

export PATH=$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;

mkdir -p $OUTPATH

echo "[BEGIN]"

set -x;
time runReadsToContig.pl $INPUT_OPTION $PLATFORM_OPTION -pre $PREFIX -bwa_options "-t $THREADS" -ref $REFDB -d $OUTPATH

time id_mapping_w_gi.pl $OUTPATH/${PREFIX}_coverage.table > $OUTPATH/$PREFIX.csv

#convert CSV to out.list format
convert_bwa2list.pl $OUTPATH/$PREFIX.csv > $OUTPATH/$PREFIX.out.list
convert_bwa2tabTree.pl < $OUTPATH/$PREFIX.csv > $OUTPATH/$PREFIX.out.tab_tree

ktImportBWA  $OUTPATH/$PREFIX.csv -o $OUTPATH/$PREFIX.krona.html
/bin/grep strain $OUTPATH/$PREFIX.out.list | awk -F "\t" '{print $2"\t"$4}' > $OUTPATH/$PREFIX.out.megan
/bin/grep species $OUTPATH/$PREFIX.out.list | awk -F "\t" '{print $2"\t"$4}' >> $OUTPATH/$PREFIX.out.megan
/bin/grep genus $OUTPATH/$PREFIX.out.list | awk -F "\t" '{print $2"\t"$4}' >> $OUTPATH/$PREFIX.out.megan

#convert BWA SAM file to read classification result
samtools view $OUTPATH/$PREFIX.sort.bam | bwa_sam2read_taxa.pl > $OUTPATH/$PREFIX.out.read_classification

set +ex;
echo "";
echo "[END] $OUTPATH $PREFIX";
