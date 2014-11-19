User Guide
##########

Command-line
============
The command line usage is as followings.

 Usage: perl runPipeline.pl [options] -c config.txt -p 'reads1.fastq reads2.fastq' -o out_directory
     Version 1.0
     Input File:
            -u            Unpaired reads, Single end reads in fastq
            
            -p            Paired reads in two fastq files and separate by space in quote

            -c            Config File
     Output:
            -o            Output directory.
  
     Options:
            -ref          Reference genome file in fasta            

            -primer       A pair of Primers sequences in strict fasta format   

            -cpu          number of CPUs (default: 8)
 
            -version      print verison


A config file (example in the below section, the GUI will generate config automatically), reads Files in fastq format, and out_directory are required when run by command line. The pipeline will run following steps. Each steps contains at least one command line scripts/programs.

Config File 
-----------
    The config file is a text file with following information. If you are going to 
do host removal, you need to build the host genome index for it and change the 
fasta file path in the config file.

[Count Fastq]
DoCountFastq=auto

[Quality Trim and Filter]
## boolean, 1=yes, 0=no
DoQC=1
##Targets quality level for trimming
q=5
##Trimmed sequence length will have at least minimum length 
min_L=50
##Average quality cutoff
avg_q=0
##"N" base cutoff.  Trimmed read has more than this number of continuous base "N" will be discarded. 
n=1
##Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence
lc=0.85
## Trim reads with adapters or contamination sequences
adapter=/PATH/adapter.fasta
## phiX filter, boolean, 1=yes, 0=no
phiX=0
## Cut # bp from 5 end before quality trimming/filtering 
5end=0
## Cut # bp from 3 end before quality trimming/filtering 
3end=0

[Host Removal]
## boolean, 1=yes, 0=no
DoHostRemoval=1
## Use more Host=  to remove multiple host reads
Host=/PATH/all_chromosome.fasta
bwaMemOptions=

[IDBA Assembly]
## boolean, 1=yes, 0=no
DoAssembly=1
idbaOptions="--pre_correction  --mink 31"

[Contigs Annotation]
## boolean, 1=yes, 0=no
DoAnnotation=1
# kingdom: Archaea Bacteria Mitochondria Viruses
kingdom=Bacteria

[Reads Mapping To Reference]
# Reads mapping to reference
DoReadsMappingReference=auto
bowtieOptions=

[Reads Mapping To Contigs]
# Reads mapping to contigs
DoReadsMappingContigs=auto

[Contigs Mapping To Reference]
# Contig mapping to reference
DoContigMapping=auto
## identity cutoff
identity=85

[Variant Analysis]
DoVariantAnalysis=auto

[ProPhage Detection]
DoProPhageDetection=1

[Reads Taxonomy Classification]
## boolean, 1=yes, 0=no
DoTaxonomy=1
## If reference genome exists, only use unmapped reads to do Taxonomy Classification. Turn on AllReads=1 will use all reads instead.
AllReads=0
enabledTools=gottcha-genDB-b,gottcha-speDB-b,gottcha-strDB-b,gottcha-genDB-v,gottcha-speDB-v,gottcha-strDB-v,metaphlan,metaphyler-srv,bwa,kraken_mini,metascope

[SNP Phylogeny]
DoSNPtree=1
## Availabe choices are Ecoli, Yersinia, Francisella, Brucella, Bacillus
SNPdbName=Ecoli

[Primer Validation]
DoPrimerValidation=1
maxMismatch=1

[Primer Adjudication]
## boolean, 1=yes, 0=no
DoPrimerDesign=0
## desired primer tm
tm=59
## reject primer having Tm < tm_diff difference with background Tm
tm_diff=5
## display # top results for each target
top=5

[Contig Blast]
DoBlast=0
BLAST_nr_DB=/PATH/nr/
BLAST_nt_DB=/PATH/nt/

[Generate JBrowse Tracks]
DoJBrowse=1

[HTML Report]
DoHTMLReport=1

Test Run
--------
    The example data set is a Ecoli MiSeq dataset which has been subsample to around 10x fold coverage reads.

In the EDGE home directory, 
    cd testData
    sh runTest.sh


Graphic User Interface
======================