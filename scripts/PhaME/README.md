# Phylogenetic and Molecular Evolution (PhaME) analysis tool

Given a reference, PhaME extracts SNPs from complete genomes, draft genomes and/or reads. 
Uses SNP multiple sequence alignment to construct a phylogenetic tree. 
Provides evolutionary analyses (genes under positive selection) using CDS SNPs.


--------------------------------------------------------------
## Version
1.0

--------------------------------------------------------------
## SYSTEM REQUIREMENTS

* Perl version >5.8

* MUMmer version 3.23 - Pairwise alignment using NUCmer 

* Bowtie version >=2.1.0 - Mapping of reads

* SAMtools version 0.1.19 and vcftools - Convert BAM files created by Bowtie

* FastTree version >=2.1.8- Contruction of phylogenetic tree

* RAxML version >=8.0.26 - Maximum likelihood construction of phylogenetic tree

* mafft version >=7.0 - For optional evolutionary analyses

* pal2nal version >=14 - For optional evolutionary analyses

* paml version >=4.8 - For optional evolutionary analyses

* HyPhy version >=2.2 - For optional evolutionary analyses

The C/C++ compiling enviroment might be required for installing dependencies. Systems may vary. Please assure that your system has the essential software building packages (e.g. build-essential for Ubuntu, XCODE for Mac...etc) installed properly before running the installing script.

PhaME was tested successfully on our Linux servers (Ubuntu 14.04.3 LTS).

--------------------------------------------------------------
### Obtaining PhaME

You can use "git" to obtain the package:

    $ git clone https://github.com/LANL-Bioinformatics/PhaME.git

### Installing PhaME

    $ cd PhaME
    $ ./INSTALL.sh
    
--------------------------------------------------------------
### Running PhaME

#### Prepare the input files
INPUT files required

* A directory with reference files which have the following file suffixes
  - *.fasta
  - *.fna
  - *.fa
  - *.gff  (optional: to analyze Coding region SNPs of a selected reference file)
    
* A working directory 
  - Contig files with the following file suffixes
     - *.contig
     - *.ctg
     - *.contigs
     
  - Reads files with the following file suffixes
     - *_R1.fastq *_R2.fastq
     - *_R1.fq *_R2.fq
     
  - A control file (e.g. [phame.ctl](https://raw.githubusercontent.com/LANL-Bioinformatics/PhaME/master/phame.ctl))

#### Test run

* Please modify the values of 'refdir' and 'workdir' in the test/phame.ctl file to corresponding absolute PhaME installed path.

* From the PhaME directory 

    $ bin/runPhaME.pl test/phame.ctl

#### OUTPUT files

* Summary
  - A SNP alignment file
  - A newick tree file
  - Text file containing the size of gaps, the core genome size, the average genome size, number of whole genome SNPs, and coding region SNPs.
  - A pairwise list of all SNPs and coordinates between references and samples
  - A matrix file that lists the number of SNPs present between genomes
  
* Directories structure
  - working directory/files
      - permutated versions of references (concatenated chromosomes)
  - working directory/results
      - All output files generated
  - working directory/results/snps
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




--------------------------------------------------------------
## CITATION

From raw reads to trees: Whole genome SNP phylogenetics across the tree of life.

Sanaa Afroz Ahmed, Chien-Chi Lo, Po-E Li, Karen W Davenport, Patrick S.G. Chain

bioRxiv doi: http://dx.doi.org/10.1101/032250

--------------------------------------------------------------
## Contact

Ahmed, Sanaa Afroz <sahmed at lanl.gov>

--------------------------------------------------------------
## ACKNOWLEDGEMENTS
This project is funded by U.S. Defense Threat Reduction Agency [R-00059-12-0 and R-00332-13-0 to P.S.G.C.].

--------------------------------------------------------------
## Copyright

Los Alamos National Security, LLC (LANS) owns the copyright to PhaME, which it identifies internally as LA-CC-xx-xxxx.  The license is GPLv3.  See [LICENSE](https://github.com/losalamos/PhaME/blob/master/LICENSE) for the full text.

