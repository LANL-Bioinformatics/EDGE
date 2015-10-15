SNPphy 
Version 1.0 07/15/2014

-----------------------------
---- GENERAL USAGE NOTES ----
-----------------------------

Given a reference extracts SNPs from contigs/draft genome and reads.
Uses SNP multiple sequence alignment to construct a phylogenetic tree.
Provides evolutionary analyses (genes under psoitive selection) using CDS SNPS.

-----------------------------
------- REQUIREMENTS --------
-----------------------------

(1) MUMmer Package
  - Pairwise alignment using NUCmer 

(2) Bowtie
  - Mapping of reads

(3) SAMtools and vcftools
  -  Convert BAM files created by Bowtie

(4) FastTree
  - Contruction of phylogenetic tree

(5) RAxML
  - Maximum likelihood construction of phylogenetic tree

(6) mafft
  - For optional evolutionary analyses

(7) pal2nal
  - For optional evolutionary analyses

(8) paml
  - For optional evolutionary analyses

(9) HyPhy
  - For optional evolutionary analyses

(10) FigTree
  - Visualization of phylogenetic tree

(11) Perl 
  - Version 5.0

-----------------------------
----- OVERVIEW OF STEPS -----
-----------------------------

INPUT files required
  - A reference directory with the following file suffixes 
    - *.fasta
    - *.fna
    - *.fa
  - A working directory 
    - Contig files with the following file suffixes
       - *.contig
       - *.ctg
       - *.contigs
    - Reads files with the following file suffixes
       - *_R1/2.fastq
       - *_R1/2.fq
    - A control file

Run runSNPphylogeny.pl from the working directory that contains the control file to extract whole genome SNPs

OUTPUT files generated
  - A SNP alignment file
  - A newick tree file
  - Text file containing the size of gaps, the core genome size, the average genome size, number of whole genome SNPs, and coding region SNPs.
  - A pairwise list of all SNPs and coordinates between references and samples
  - A matrix file that lists the number of SNPs present between genomes
  - 

DIRECTORIES generated
  - working directory/files
    - permutated versions of references (concatenated chromosomes)
  - working directory/reseults
    - All output files generated
  - working directory/reseults/snps
    - SNP coordinate files generated from NUCmer and Bowtie
  - working directory/results/gaps
    - Gap coordinate files generated from NUCmer and bowtie
  - working directory/results/stats
    - Intermediate stat files generated when parsing NUCmer and Bowtie results
  - working directory/results/temp
    - Temporary files generated
  - working directory/results/PSgenes
    - All gene fasta files  that contain at least 1 SNP
  - working directory/results/paml
    - PAML results
  - working directory/results/hyphy
    - HyPhy results

-----------------------------
------ STEPS IN DETAIL ------
-----------------------------

----- Step 1: Modify the Control file -----

