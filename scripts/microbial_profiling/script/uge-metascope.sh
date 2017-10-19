#!/bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe smp 60
#$ -m abe
#$ -M po-e@lanl.gov
#$ -j y

set -e;

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -d <DATABASE> -o <OUTDIR> -p <PREFIX> [OPTIONS]

OPTIONS:
   -i      Input a FASTQ file (singleton or pair-end sequences separated by comma)
   -d      Database (not avaliable)
   -o      Output directory
   -p      Output prefix
   -t      Number of threads (default: 24)
   -e      Extra-option: long or pacbio
   -s      Sequencer (pacbio, roche, illumina or ion)
   -h      help
EOF
}

FASTQ=
PREFIX=
OUTPATH=
THREADS=24
SEQUENCER="illumina"
INPUT_OPTION=
EXTRAOPT=

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
        e) EXTRAOPT=$OPTARG
           ;;
        s) SEQUENCER=$OPTARG
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

export PATH=$EDGE_HOME/thirdParty/metascope_plus:$EDGE_HOME/thirdParty/metascope_plus/metascope/bin:$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;
mkdir -p $OUTPATH

echo "[BEGIN]"

set -x;
cd $EDGE_HOME/thirdParty/metascope_plus;
#Usage: metascope_plus -in <input FastQ> -o <output file name> -s <sequencer, one of pacbio,roche,illumina,ion>
metascope_plus.sh $EXTRAOPT -in $FASTQ -o $OUTPATH/$PREFIX -t $THREADS -s $SEQUENCER
set +e;

cd -;

#parse mpln
convert_metascope2list.pl < $OUTPATH/$PREFIX.metascope.log > $OUTPATH/$PREFIX.out.list
convert_list2tabTree.pl < $OUTPATH/$PREFIX.out.list > $OUTPATH/$PREFIX.out.tab_tree

# Make Krona plot
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

# Not used 101817.
# cat $OUTPATH/$PREFIX.out.list | awk -F\\t "{if(\$4>0 && \$4!=\"ASSIGNED\") print \$5\"\t\"\$4}" > $OUTPATH/$PREFIX.out.krona
# This is the old way to make Krona plots. See above for current way (consistent with other methods)
#ktImportTaxonomy -t 1 -s 2 -o $OUTPATH/$PREFIX.krona.html $OUTPATH/$PREFIX.out.krona

set +x;
echo "";
echo "[END] $OUTPATH $PREFIX";

