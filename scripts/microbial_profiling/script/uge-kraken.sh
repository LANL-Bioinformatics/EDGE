#!/bin/bash
#$ -cwd
#$ -l h_vmem=2.6G
#$ -m abe
#$ -j y

set -e;

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -d <DATABASE> -o <OUTDIR> -p <PREFIX> [OPTIONS]

OPTIONS:
   -i      Input a FASTQ file (singleton or pair-end sequences separated by comma)
   -d      Database
   -n      Additional options for kraken
   -o      Output directory
   -p      Output prefix
   -t      Number of threads (default: 24)
   -h      help
EOF
}

FASTQ=
REFDB=$EDGE_HOME/database/minikraken
PREFIX=
OUTPATH=
THREADS=24
OPTIONS=

while getopts "i:d:o:p:t:n:h" OPTION
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
        n) OPTIONS=$OPTARG
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
time kraken --threads $THREADS --db $REFDB $OPTIONS --output $OUTPATH/$PREFIX.classification.csv --fastq-input $FASTQ
time kraken-report --db $REFDB $OUTPATH/$PREFIX.classification.csv > $OUTPATH/$PREFIX.report.csv
set +e;

#generate out.list
convert_krakenRep2list.pl < $OUTPATH/$PREFIX.report.csv > $OUTPATH/$PREFIX.out.list
convert_krakenRep2tabTree.pl < $OUTPATH/$PREFIX.report.csv > $OUTPATH/$PREFIX.out.tab_tree

# Make Krona plot
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

# This is old 101817. See above for current Krona plot creation (consistent with the method used for other tools).
#generate krona plot
#cat $OUTPATH/$PREFIX.out.list | awk -F\\t "{if(\$4>0 && \$4!=\"ASSIGNED\") print \$5\"\t\"\$4}" > $OUTPATH/$PREFIX.out.krona
#ktImportTaxonomy -t 1 -s 2 -o $OUTPATH/$PREFIX.krona.html $OUTPATH/$PREFIX.out.krona

#preparing megan CSV file
cat $OUTPATH/$PREFIX.out.list | awk -F\\t "{if(\$4>0 && \$4!=\"ASSIGNED\") print \$2\"\t\"\$4}" > $OUTPATH/$PREFIX.out.megan

set +x;
echo "";
echo "[END] $OUTPATH $PREFIX";
