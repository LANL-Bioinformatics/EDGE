#!/bin/bash
#$ -l h_vmem=2.6G
#$ -j y
#$ -cwd

usage(){
cat << EOF
USAGE: $0 -i <FASTQ> -o <OUTDIR> -p <PREFIX> -d <DB> [OPTIONS]

ARGUMENTS:
   -i      Input a FASTQ file or pair-ended FASTAQ files sperated by comma
   -o      Output directory
   -p      Output prefix
   -d      Database

OPTIONS:
   -b      background JSON file
   -t      Number of threads. [default is 4]
   -a      additional options
   -r      The field will be used to calculate relative 
           abundance. You can specify one of the following
           fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED",
           "READ_COUNT" and "LINEAR_DOC". [default: LINEAR_DOC]      
   -S      minScore
   -R      minReads
   -B      minRsnb
   -L      minLen
   -h      help
EOF
}

###########################
#
# Default values
#
###########################

FASTQ=
OUTPATH=
PREFIX=
DB=$EDGE_HOME/database/PanGIA/NCBI_genomes_refseq86*.fa
BG_JSON=
BG_JSON_OPT=
READ_COND="R_MAT"
THREADS=4
OPTIONS="-sb -se"
FIELD="DEPTH_COV"
MIN_SCORE=0
MIN_READ=3
MIN_RSNB=1
MIN_LEN=50

while getopts "i:o:p:d:b:c:t:a:r:S:R:B:L:h" OPTION
do
     case $OPTION in
        i) FASTQ=$OPTARG
           ;;
        o) OUTPATH=$OPTARG
           ;;
        p) PREFIX=$OPTARG
           ;;
        d) DB=$OPTARG
           ;;
        b) BG_JSON=$OPTARG
           ;;
        c) READ_COND=$OPTARG
           ;;
        t) THREADS=$OPTARG
           ;;
        a) OPTIONS=$OPTARG
           ;;
        r) FIELD=$OPTARG
           ;;
        S) MIN_SCORE=$OPTARG
           ;;
        R) MIN_READ=$OPTARG
           ;;
        B) MIN_RSNB=$OPTARG
           ;;
        L) MIN_LEN=$OPTARG
           ;;
        h) usage
           exit
           ;;
     esac
done

## path
export PATH=$EDGE_HOME/thirdParty/pangia:$EDGE_HOME/bin:$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;
mkdir -p $OUTPATH
set -xe;

# Background JSON
if [[ -s $BG_JSON ]]
then
    BG_JSON_OPT="-lb $BG_JSON"
fi

if [[ -z "$HOME" ]]
then
	#for parallel in pangia 
	HOME=$EDGE_HOME
fi

pangia.py --debug -r $FIELD -i $FASTQ -t $THREADS -o $OUTPATH -p $PREFIX -d $DB $BG_JSON_OPT -ms $MIN_SCORE -mr $MIN_READ -mb $MIN_RSNB -ml $MIN_LEN $OPTIONS

# prepare pangia-vis
#DIRMD5=`echo ${OUTPATH/\/ReadsBasedAnalysis*/} | rev | cut -d'/' -f1 | rev`
MERGED_SAM_DIR=$OUTPATH/${PREFIX}_tmp/merged_sam
LOCAL_PVIA_DATA=$OUTPATH/pangia-vis

mkdir -p $LOCAL_PVIA_DATA/report

cp $OUTPATH/$PREFIX.report.tsv $LOCAL_PVIA_DATA/report.tsv
cp $OUTPATH/$PREFIX.pangia.log $LOCAL_PVIA_DATA/report.pangia.log

# generate scaled-down genome depth files for PanGIA-VIS, save to $OUTPATH/pangia-vis then symlink to global pangia-vis
find $MERGED_SAM_DIR -type f -name '*.depth' | cut -d"|" -f3 | sort | uniq > $OUTPATH/${PREFIX}_tmp/unique_taxid.txt
parallel "$EDGE_HOME/thirdParty/pangia/pangia-vis/scripts/depth_scale_down.py $MERGED_SAM_DIR/*\|{}\|*.depth > $LOCAL_PVIA_DATA/report/{}.depth.scaledown" :::: $OUTPATH/${PREFIX}_tmp/unique_taxid.txt

# remove temporary files
find $OUTPATH/${PREFIX}_tmp -type f -exec rm {} +
rm -rf $OUTPATH/${PREFIX}_tmp


awk -F\\t '{if(NR==1){out=$1"\t"$2"\tROLLUP\tASSIGNED"; { for(i=3;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.report.tsv > $OUTPATH/$PREFIX.out.list
awk -F\\t '{if(NR>1&&$16==""){out=$1"\t"$2"\t"$14"\t"; { for(i=3;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.report.tsv >> $OUTPATH/$PREFIX.out.list

# prepare lineage and score files
awk -F\\t '{if($1=="species"){print $2"\t"$12}}' $OUTPATH/$PREFIX.report.tsv > $OUTPATH/$PREFIX.out.tab_tree.score

cat $OUTPATH/$PREFIX.report.tsv | pangia_report2lineage.py -dp ${DB%/*} > $OUTPATH/$PREFIX.out.tab_tree

#generate KRONA chart
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

set +xe;
echo "";
echo "[END] $OUTPATH $PREFIX";
