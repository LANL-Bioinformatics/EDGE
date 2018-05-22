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
   -d      Database

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
BWAMETHOD=" -k 30 -T 0 -B 100 -O 100 -E 100"

TRIM_MINL=
TRIM_FIXL=30
TRIM_MINQ=20
TRIM_ASCI=33
DB=$EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species

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
        d) DB=$OPTARG
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
export PATH=$EDGE_HOME/bin:$EDGE_HOME/scripts/microbial_profiling/script:$EDGE_HOME/scripts:$PATH;

if [[ -z "$OUTPATH" || -z "$PREFIX" ]]
then
     usage;
     exit 1;
fi

set -ex;

mkdir -p $OUTPATH

#############################
#
# splitrim FASTQ
#
#############################

# splitrim 1st FASTQ FILE
if [[ -z "$PRE_SPLITRIM" && ! -s "$OUTPATH/../../sequence_processed/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}/${PREFIX}_splitrim.fastq" ]]
then
	$EDGE_HOME/thirdParty/gottcha/bin/splitrim --inFile=$FASTQ --ascii=$TRIM_ASCI --fixL=$TRIM_FIXL --recycle --minQ=$TRIM_MINQ --prefix=$PREFIX --outPath=$OUTPATH/../../sequence_processed/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}
else
	echo "[TRIM] Skip splitrim step...";
fi

# pre-splitrim
if [[ -z "$PRE_SPLITRIM" ]]
then
	SPLITRIM="$OUTPATH/../../sequence_processed/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}"
else
	# checking the existence of splitrimmed and stats files
	if [[ -e "${PRE_SPLITRIM}/${PREFIX}_splitrim.fastq" && -e "${PRE_SPLITRIM}/${PREFIX}_splitrim.stats.txt" ]]
	then
		SPLITRIM="${PRE_SPLITRIM}"
	else
		echo "[TRIM] FATAL: missing ${PRE_SPLITRIM}/${PREFIX}_splitrim.fastq or ${PRE_SPLITRIM}/${PREFIX}_splitrim.stats.txt";
		exit 1;
	fi
fi

gottcha.pl --mode all --bwaOpt "$BWAMETHOD" -i $FASTQ -t $THREADS -stDir $SPLITRIM --outdir $OUTPATH -p $PREFIX --database $DB --dumpSam 

set +e;

#convert gottchaFilterRollupAbundance result to list
awk -F\\t '{if(NR==1){out=$1"\t"$2"\tROLLUP\tASSIGNED"; { for(i=4;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.gottcha.tsv > $OUTPATH/$PREFIX.out.list
awk -F\\t '{if(NR>1){out=$1"\t"$2"\t"$3"\t"; { for(i=4;i<=NF;i++){out=out"\t"$i}}; print out;}}' $OUTPATH/$PREFIX.gottcha.tsv >> $OUTPATH/$PREFIX.out.list

#generate tab_tree file
gottcha_sum_lineage.pl $LVL < $OUTPATH/${PREFIX}_temp/$PREFIX.lineage.tsv > $OUTPATH/$PREFIX.out.tab_tree

#generate KRONA chart
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

set +x;
echo "";
echo "[END] $OUTPATH $PREFIX";

