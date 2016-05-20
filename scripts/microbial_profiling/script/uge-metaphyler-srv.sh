#!/bin/bash
#$ -cwd
#$ -l h_vmem=5G
#$ -m abe
#$ -j y
set -e;

usage(){
cat << EOF
USAGE: $0 -i <FASTA> -o <OUTDIR> -p <PREFIX> [OPTIONS]

OPTIONS:
   -i      Input a FASTA file
   -o      Output directory
   -p      Output prefix
   -t      Number of threads
   -c      cut-off score (default 0.9)
   -h      help
EOF
}

FASTA=
PREFIX=
OUTPATH=
THREADS=24

while getopts "i:o:p:t:c:h" OPTION
do
     case $OPTION in
        i) FASTA=$OPTARG
           ;;
        o) OUTPATH=$OPTARG
           ;;
        p) PREFIX=$OPTARG
           ;;
        t) THREADS=$OPTARG
           ;;
        c) CUTOFF=$OPTARG
           ;;
        h) usage
           exit
           ;;
     esac
done

if [[ -z "$FASTA" || -z "$OUTPATH" || -z "$PREFIX" ]]
then
     usage;
     exit 1;
fi

mkdir -p $OUTPATH

export PATH=$EDGE_HOME/bin/microbial_profiling/script:$PATH;

set -x;
#run metaphyler.pl
time metaphyler.pl ${THREADS} ${FASTA} ${OUTPATH}/${PREFIX}
set +e;

#parse mpln
convert_metaphyler_srv2metaphlan.pl $OUTPATH/$PREFIX.classify.tab > $OUTPATH/$PREFIX.out.mpln
convert_metaphylerSrv2tabTree.pl < $OUTPATH/$PREFIX.classify.tab > $OUTPATH/$PREFIX.out.tab_tree

mpln2krona.pl -l g -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.krona
ktImportText $OUTPATH/$PREFIX.out.krona -o $OUTPATH/$PREFIX.krona.html
krona_portable.pl -inhtml $OUTPATH/$PREFIX.kronaInternet.html  -outhtml $OUTPATH/$PREFIX.krona.html
rm $OUTPATH/$PREFIX.kronaInternet.html

mpln2list.pl  -l g -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.list
mpln2megan.pl -l g -i $OUTPATH/$PREFIX.out.mpln > $OUTPATH/$PREFIX.out.megan

set +x;
echo "";
echo "[END] $OUTPATH $PREFIX";
