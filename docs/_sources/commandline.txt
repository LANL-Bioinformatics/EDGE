Command Line Interface (CLI)
############################

The command line usage is as followings:: 

    Usage: perl runPipeline.pl [options] -c config.txt -p 'reads1.fastq reads2.fastq' -o out_directory
    Version 1.1
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


A config file (example in the below section, the :doc:`Graphic User Interface (GUI) <gui>` will generate config automatically), reads Files in fastq format, and a output directory are required when run by command line. Based on the configuration file, if all modules are turned on, EDGE will run the following steps. Each step contains at least one command line scripts/programs.

1. Data QC 
2. Host Removal QC 
3. *De novo* Assembling
4. Reads Mapping To Contig
5. Reads Mapping To Reference Genomes 
6. Taxonomy Classification on All Reads or unMapped to Reference Reads
7. Map Contigs To Reference Genomes 
8. Variant Analysis
9. Contigs Taxonomy Classification
10. Contigs Annotation
11. ProPhage detection
12. PCR Assay Validation 
13. PCR Assay Adjudication
14. Phylogenetic Analysis
15. Generate JBrowse Tracks
16. HTML report 


.. _config_example:

Configuration File 
==================

The config file is a text file with the following information. If you are going to do host removal, you need to :ref:`build host index<build-host-index>` for it and change the fasta file path in the config file. ::

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
	similarity=90

	[Assembly]
	## boolean, 1=yes, 0=no
	DoAssembly=1
	##Bypass assembly and use pre-assembled contigs
	assembledContigs=
	minContigSize=200
	## spades or idba_ud
	assembler=idba_ud
	idbaOptions="--pre_correction  --mink 31"
	## for spades
	singleCellMode=
	pacbioFile=
	nanoporeFile=

	[Reads Mapping To Contigs]
	# Reads mapping to contigs
	DoReadsMappingContigs=auto

	[Reads Mapping To Reference]
	# Reads mapping to reference
	DoReadsMappingReference=0
	bowtieOptions=
	# reference genbank or fasta file
	reference=
	MapUnmappedReads=0

	[Reads Taxonomy Classification]
	## boolean, 1=yes, 0=no
	DoReadsTaxonomy=1
	## If reference genome exists, only use unmapped reads to do Taxonomy Classification. Turn on AllReads=1 will use all reads instead.
	AllReads=0
	enabledTools=gottcha-genDB-b,gottcha-speDB-b,gottcha-strDB-b,gottcha-genDB-v,gottcha-speDB-v,gottcha-strDB-v,metaphlan,bwa,kraken_mini

	[Contigs Mapping To Reference]
	# Contig mapping to reference
	DoContigMapping=auto
	## identity cutoff
	identity=85
	MapUnmappedContigs=0
	
	[Variant Analysis]
	DoVariantAnalysis=auto

	[Contigs Taxonomy Classification]
	DoContigsTaxonomy=1

	[Contigs Annotation]
	## boolean, 1=yes, 0=no
	DoAnnotation=1
	# kingdom: Archaea Bacteria Mitochondria Viruses
	kingdom=Bacteria
	contig_size_cut_for_annotation=700
	## support tools: Prokka or RATT
	annotateProgram=Prokka
	annotateSourceGBK=

	[ProPhage Detection]
	DoProPhageDetection=1

	[Phylogenetic Analysis]
	DoSNPtree=1
	## Availabe choices are Ecoli, Yersinia, Francisella, Brucella, Bacillus
	SNPdbName=Ecoli
	## FastTree or RAxML
	treeMaker=FastTree
	## SRA accessions ByrRun, ByExp, BySample, ByStudy
	SNP_SRA_ids=

	[Primer Validation]
	DoPrimerValidation=1
	maxMismatch=1
	primer=

	[Primer Adjudication]
	## boolean, 1=yes, 0=no
	DoPrimerDesign=0
	## desired primer tm
	tm_opt=59
	tm_min=57
	tm_max=63
	## desired primer length
	len_opt=18
	len_min=20
	len_max=27
	## reject primer having Tm < tm_diff difference with background Tm
	tm_diff=5
	## display # top results for each target
	top=5

	[Generate JBrowse Tracks]
	DoJBrowse=1

	[HTML Report]
	DoHTMLReport=1

Test Run
========
    EDGE provides an example data set which is an E. coli MiSeq dataset and has been subsampled to ~10x fold coverage reads.

In the EDGE home directory, ::

    cd testData
    sh runTest.sh
    
.. figure:: img/commandline_screen.png
   :width: 400px
   :alt: Snapshot from the terminal.
   
   Snapshot from the terminal.

See :doc:`output`

Descriptions of each module
===========================

Each module comes with default parameters and user can see the optional parameters by entering the program name with –h or -help flag without any other arguments.

1. **Data QC**

  * Required step? **No**

  * Command example ::
       
       perl $EDGE_HOME/scripts/illumina_fastq_QC.pl  -p 'Ecoli_10x.1.fastq Ecoli_10x.2.fastq'  -q 5 -min_L 50 -avg_q 5 -n 0 -lc 0.85 –d QcReads -t 10
    
  * What it does
  
    * Quality control
    * Read filtering
    * Read trimming  
    
  * Expected input
  
    * Paired-end/Single-end reads in FASTQ format  
    
  * Expected output
  
    * QC.1.trimmed.fastq           
    * QC.2.trimmed.fastq  
    * QC.unpaired.trimmed.fastq
    * QC.stats.txt
    * QC_qc_report.pdf    

2. **Host Removal QC** 

  * Required step? **No** 

  * Command example ::
  
       perl $EDGE_HOME/scripts/host_reads_removal_by_mapping.pl  -p 'QC.1.trimmed.fastq QC.2.trimmed.fastq' -u QC.unpaired.trimmed.fastq -ref human_chromosomes.fasta -o QcReads -cpu 10       
  
  * What it does
  
    * Read filtering
       
  * Expected input
  
    * Paired-end/Single-end reads in FASTQ format

  * Expected output
  
    * host_clean.1.fastq  
    * host_clean.2.fastq  
    * host_clean.mapping.log  
    * host_clean.unpaired.fastq  
    * host_clean.stats.txt    

3. **IDBA Assembling**

  * Required step? **No**

  * Command example ::
  
      fq2fa --merge host_clean.1.fastq  host_clean.2.fastq  pairedForAssembly.fasta             
      idba_ud  --num_threads 10 -o AssemblyBasedAnalysis/idba --pre_correction pairedForAssembly.fasta       

  * What it does
  
    * Iterative kmers de novo Assembly, it performs well on isolates as well as metagenomes.  It may not work well on very large genomes.
    
  * Expected input
  
    * Paired-end/Single-end reads in FASTA format

  * Expected output
  
    * contig.fa
    * scaffold.fa  (input paired end)    

4. **Reads Mapping To Contig**

  * Required step? **No**

  * Command example ::

      perl $EDGE_HOME/scripts/runReadsToContig.pl  -p 'host_clean.1.fastq host_clean.2.fastq' -d AssemblyBasedAnalysis/readsMappingToContig -pre readsToContigs  -ref AssemblyBasedAnalysis/contigs.fa 

  * What it does
  
    * Mapping reads to assembled contigs
  
  * Expected input
  
    * Paired-end/Single-end reads in FASTQ format
    * Assembled Contigs in Fasta format
    * Output Directory
    * Output prefix

  * Expected output
    
    * readsToContigs.alnstats.txt
    * readsToContigs_coverage.table
    * readsToContigs_plots.pdf
    * readsToContigs.sort.bam
    * readsToContigs.sort.bam.bai

5. **Reads Mapping To Reference Genomes** 

  * Required step? **No**
  
  * Command example::

      perl $EDGE_HOME/scripts/runReadsToGenome.pl  -p 'host_clean.1.fastq host_clean.2.fastq'  -d ReadsBasedAnalysis -pre readsToRef -ref Reference.fna  
      
  * What it does

    * Mapping reads to reference genomes
    * SNPs/Indels calling     

  * Expected input
    
    * Paired-end/Single-end reads in FASTQ format
    * Reference genomes in Fasta format
    * Output Directory
    * Output prefix

  * Expected output
    
    * readsToRef.alnstats.txt
    * readsToRef_plots.pdf
    * readsToRef_refID.coverage        
    * readsToRef_refID.gap.coords
    * readsToRef_refID.window_size_coverage  
    * readsToRef.ref_windows_gc.txt
    * readsToRef.raw.bcf
    * readsToRef.sort.bam
    * readsToRef.sort.bam.bai
    * readsToRef.vcf

6. **Taxonomy Classification on All Reads or unMapped to Reference Reads**

  * Required step? **No**

  * Command example::
  
      perl $EDGE_HOME/scripts/microbial_profiling/microbial_profiling_configure.pl $EDGE_HOME/scripts/microbial_profiling/microbial_profiling.settings.tmpl gottcha-speDB-b > microbial_profiling.settings.ini
      perl $EDGE_HOME/scripts/microbial_profiling/microbial_profiling.pl -o  Taxonomy -s microbial_profiling.settings.ini -c 10 UnmappedReads.fastq 
      
  * What it does
    
    * Taxonomy Classification using multiple tools, including BWA mapping to NCBI Refseq, metaphlan, kraken, GOTTCHA.
    * Unify varies output format and generate reports

  * Expected input
    
    * Reads in FASTQ format
    * Configuration text file (generated by microbial_profiling_configure.pl)

  * Expected output
    
    * Summary EXCEL and text files.
    * Heatmaps tools comparison
    * Radarchart tools comparison
    * Krona and tree-style plots for each tool.

7. **Map Contigs To Reference Genomes**

  * Required step? **No**
  
  * Command example::
      
      perl $EDGE_HOME/scripts/nucmer_genome_coverage.pl  -e 1 -i 85 –p contigsToRef Reference.fna contigs.fa       

  * What it does
    
    * Mapping assembled contigs to reference genomes
    * SNPs/Indels calling     
 
  * Expected input
 
    * Reference genome in Fasta Format
    * Assembled contigs in Fasta Format
    * Output prefix

  * Expected output
    
    * contigsToRef_avg_coverage.table
    * contigsToRef.delta
    * contigsToRef_query_unUsed.fasta
    * contigsToRef.snps
    * contigsToRef.coords
    * contigsToRef.log
    * contigsToRef_query_novel_region_coord.txt
    * contigsToRef_ref_zero_cov_coord.txt

8. **Variant Analysis**

  * Required step?	**No** 

  * Command example:: 
  
      perl $EDGE_HOME/scripts/SNP_analysis.pl -genbank Reference.gbk -SNP contigsToRef.snps -format nucmer
      perl $EDGE_HOME/scripts/gap_analysis.pl -genbank Reference.gbk -gap  contigsToRef_ref_zero_cov_coord.txt

  * What it does
    
    * Analyze variants and gaps regions using annotation file.
    
  * Expected input
    
    * Reference in GenBank format
    * SNPs/INDELs/Gaps files from “Map Contigs To Reference Genomes“ 

  * Expected output
  
    * contigsToRef.SNPs_report.txt
    * contigsToRef.Indels_report.txt
    * GapVSReference.report.txt

9. **Contigs Taxonomy Classification**

  * Required step? **No**

  * Command example::

      perl $EDGE_HOME/scripts/contig_classifier_by_bwa/contig_classifier_by_bwa.pl --db $EDGE_HOME/database/bwa_index/NCBI-Bacteria-Virus.fna --threads 10 --prefix OuputCT --input contigs.fa

  * What it does
  
    * Taxonomy Classification on contigs using BWA mapping to NCBI Refseq
  
  * Expected input
    
    * Contigs in Fasta format
    * NCBI Refseq genomes bwa index
    * Output prefix
    
  * Expected output
    
    * prefix.assembly_class.csv
    * prefix.assembly_class.top.csv
    * prefix.ctg_class.csv
    * prefix.ctg_class.LCA.csv
    * prefix.ctg_class.top.csv
    * prefix.unclassified.fasta

10. **Contig Annotation**

  * Required step? **No**
  
  * Command example::
  
      prokka --force --prefix PROKKA --outdir Annotation contigs.fa 
      
  * What it does

    * The rapid annotation of prokaryotic genomes.

  * Expected input
    
    * Assembled Contigs in Fasta format
    * Output Directory
    * Output prefix

  * Expected output

    * It produces GFF3, GBK and SQN files that are ready for editing in Sequin and ultimately submitted to Genbank/DDJB/ENA.

11. **ProPhage detection**

  * Required step? **No**

  * Command example:: 
  
      perl $EDGE_HOME/scripts/phageFinder_prepare.pl -o Prophage –p Assembly Annotation/PROKKA.gff Annotation/PROKKA.fna
      $EDGE_HOME/thirdParty/phage_finder_v2.1/bin/phage_finder_v2.1.sh Assembly 

  * What it does

    * Identify and classify prophages within prokaryotic genomes.
       
  * Expected input
  
    * Annotated Contigs GenBank file
    * Output Directory
    * Output prefix

  * Expected output
  
    * phageFinder_summary.txt

12. **PCR Assay Validation**
      
  * Required step? **No**
  
  * Command example::
  
      perl $EDGE_HOME/scripts/pcrValidation/validate_primers.pl -ref contigs.fa -primer primers.fa -mismatch 1 -output AssayCheck 

  * What it does

    * In silico PCR primer validation by sequence alignment. 

  * Expected input

    * Assembled Contigs/Reference in Fasta format
    * Output Directory
    * Output prefix

  * Expected output

    * pcrContigValidation.log
    * pcrContigValidation.bam

13. **PCR Assay Adjudication**  
      
  * Required step? **No**
  
  * Command example::
  
      perl $EDGE_HOME/scripts/pcrAdjudication/pcrUniquePrimer.pl --input contigs.fa  --gff3 PCR.Adjudication.primers.gff3

  * What it does

    * Design unique primer pairs for input contigs.

  * Expected input
  
    * Assembled Contigs in Fasta format
    * Output gff3 file name

  * Expected output

    * PCR.Adjudication.primers.gff3
    * PCR.Adjudication.primers.txt

14.  **Phylogenetic Analysis**
      
  * Required step? **No**
   
  * Command example::

      perl $EDGE_HOME/scripts/prepare_SNP_phylogeny.pl -o output/SNP_Phylogeny/Ecoli -tree FastTree -db Ecoli -n output -cpu 10 -p QC.1.trimmed.fastq QC.2.trimmed.fastq -c contigs.fa -s QC.unpaired.trimmed.fastq
      perl $EDGE_HOME/scripts/SNPphy/runSNPphylogeny.pl output/SNP_Phylogeny/Ecoli/SNPphy.ctrl

  * What it does

    * Perform SNP identification against selected pre-built SNPdb or selected genomes
    * Build SNP based multiple sequence alignment for all and CDS regions
    * Generate Tree file in newick/PhyloXML format
    
  * Expected input
    
    * SNPdb path or genomesList
    * Fastq reads files
    * Contig files

  * Expected output
  
    * SNP based phylogentic multiple sequence alignment
    * SNP based phylogentic tree in newick/PhyloXML format.
    * SNP information table 

15. **Generate JBrowse Tracks**
      
  * Required step? **No**
    
  * Command example::

      perl $EDGE_HOME/scripts/edge2jbrowse_converter.pl --in-ref-fa Reference.fna --in-ref-gff3 Reference.gff --proj_outdir EDGE_project_dir

  * What it does
    
    * Convert several EDGE outputs into JBrowse tracks for visualization for contigs and reference, respectively.
  
  * Expected input
  
    * EDGE project output Directory

  * Expected output
    
    * EDGE post-processed files for JBrowse tracks in the JBrowse directory.
    * Tracks configuration files in the JBrowse directory.

16. **HTML Report**

  * Required step? **No**
   
  * Command example::
  
      perl $EDGE_HOME/scripts/munger/outputMunger_w_temp.pl EDGE_project_dir

  * What it does

    * Generate statistical numbers and plots in an interactive html report page. 

  * Expected input
    
    * EDGE project output Directory
  
  * Expected output
  
    * report.html


Other command-line utility scripts
==================================

1. To extract certain taxa fasta from contig classification result::

    cd /home/edge_install/edge_ui/EDGE_output/41/AssemblyBasedAnalysis/Taxonomy
    perl /home/edge_install/scripts/contig_classifier_by_bwa/extract_fasta_by_taxa.pl -fasta ../contigs.fa -csv ProjectName.ctg_class.top.csv -taxa "Enterobacter cloacae” > Ecloacae.contigs.fa

2. To extract unmapped/mapped reads fastq from the bam file::

    cd /home/edge_install/edge_ui/EDGE_output/41/AssemblyBasedAnalysis/readsMappingToContig
    # extract unmapped reads
    perl /home/edge_install/scripts/bam_to_fastq.pl -unmapped readsToContigs.sort.bam
    # extract mapped reads
    perl /home/edge_install/scripts/bam_to_fastq.pl -mapped readsToContigs.sort.bam

3. To extract mapped reads fastq of a specific contig/reference from the bam file:: 

    cd /home/edge_install/edge_ui/EDGE_output/41/AssemblyBasedAnalysis/readsMappingToContig
    perl /home/edge_install/scripts/bam_to_fastq.pl -id ProjectName_00001 -mapped readsToContigs.sort.bam 
    