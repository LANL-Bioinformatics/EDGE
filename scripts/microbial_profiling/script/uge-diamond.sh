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
   -n      Additional options for diamond
   -o      Output directory
   -p      Output prefix
   -t      Number of threads (default: 16)
   -h      help
EOF
}

FASTQ=

# This is the database without many viral sequences
#REFDB=$EDGE_HOME/database/diamond/RefSeq_Release83.nr_protein.faa.dmnd

# This is the database with the NCBI viral database sequences added
REFDB=$EDGE_HOME/database/diamond/RefSeq_Release83.nr_protein_withRefSeq_viral_102317.protein.faa.dmnd

PREFIX=
OUTPATH=
THREADS=16
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

echo "[BEGIN DIAMOND]"


# -f 102 activates the taxonomy classification module of diamond and is required to make this work.
time diamond blastx -p $THREADS -q $FASTQ -d $REFDB --taxonmap $EDGE_HOME/database/diamond/taxonomy/prot.accession2taxid.gz --taxonnodes $EDGE_HOME/database/diamond/taxonomy/nodes.dmp -f 102 -o $OUTPATH/$PREFIX.diamondRawOutput.txt

# Make the taxonomy list file.
python $EDGE_HOME/scripts/microbial_profiling/script/convert_diamond2list.py -tp $EDGE_HOME/database/diamond/taxonomy < $OUTPATH/$PREFIX.diamondRawOutput.txt > $OUTPATH/$PREFIX.out.list

# Convert for krona plot and generate taxon-based krona plot. Use -s 2 flag for ignoring abundance data and splitting wedge size evenly among all listed taxa (if taxa are listed >1 time,
# then they will get a larger wedge. This is the same as ktImportText using the tab_tree file.
# Use flag -m 2 to scale the wedge size to each taxa's relative abundance (better matches the table output in EDGE and the dendrogram.
#python $EDGE_HOME/scripts/microbial_profiling/script/convert_diamond2krona.py -i $OUTPATH/$PREFIX.out.list > $OUTPATH/$PREFIX.out.krona

# Not used. See below for Krona plot making. An alternative way to make the krona plots.
#ktImportTaxonomy -t 1 -m 2 -o $OUTPATH/$PREFIX.krona.html $OUTPATH/$PREFIX.out.krona

# Convert into tab_tree file
python $EDGE_HOME/scripts/microbial_profiling/script/convert_diamond2tabTree.py -i $OUTPATH/$PREFIX.out.list -dp $EDGE_HOME/database/diamond/taxonomy > $OUTPATH/$PREFIX.out.tab_tree
# Make Krona plot
ktImportText  $OUTPATH/$PREFIX.out.tab_tree -o $OUTPATH/$PREFIX.krona.html

# Convert into a megan file (if necessary for EDGE)
python $EDGE_HOME/scripts/microbial_profiling/script/convert_diamond2megan.py -i $OUTPATH/$PREFIX.out.list > $OUTPATH/$PREFIX.out.megan

