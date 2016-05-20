This is an amplicon analysis pipeline based on Qiime v1.9.1 developed by Genome Science group of Los Alamos National Laboratory.

# REQUIREMENTS

* qiime version 1.9.1 (http://qiime.org/install/install.html)
* fasts-join from ea-utils 1.1.2-537 (https://code.google.com/p/ea-utils/)
* usearch version > 7 (optional for UPARSE, http://www.drive5.com/usearch/)

# USAGE 

Usage: perl qiime_pipeline.pl [options] -p reads1.fastq reads2.fastq -m mapping.txt -o out_directory 

Version 0.4

    Input File: (can use more than once)
            -p or -u      -u     <File> Unpaired reads fastq
                          -p     <File> Paired reads in two fastq files and separate by space
               or -d      -d     <PATH/directory> Contains Demultiplexed fastq files. To use this option, 
                                 the mapping file need a 'Files' column with filenames for each sampleID.
            
            -m            <File> MAPPING file.  #SampleID BarcodeSequence LinkerPrimerSequence sampleType ... Description
			              comma separated for multiple mapping file
    Output:
            -o            <PATH> Output directory.
    Options:
            -target       Greengenes, SILVA, or ITS. [default: Greengenes]
                          Greengenes is for 16s.
                          SILVA is for 16s/18s.
                          ITS is from https://unite.ut.ee/ and for fungal rDNA ITS sequences. 

            -b            <File> Barcodes fastq
            
            -barcode_len  Length of Barcode [default: 6]

            -demultiplex_fa  DE-MULTIPLEX FASTA file 
                          provide multiple previous demultiplex fasta by comma
			              separated files and will bypass the split fastq step. 
			              The fasta file will be concatenated as input for OTU clustering.
			               e.g: /path/run1/seqs.fna,/path/run2/seqs.fna

            -UPARSE	      <bool> use UPARSE pipeline clusters NGS amplicon reads into OTUs
			              Edgar, R.C. (2013) Nature Methods Pubmed:23955772
 
            -q            PHRED_QUALITY_THRESHOLD
                          the maximum unacceptable Phred quality score (e.g.,
                          for Q20 and better, specify -q 19) [default: 3]

            -phred_offset The ascii offset to use when decoding phred scores (33 or 64)
                          [default: determined automatically]
	
            -n            SEQUENCE_MAX_N
                          maximum number of N characters allowed in a sequence
                          to retain it -- this is applied after quality
                          trimming, and is total over combined paired end reads
                          if applicable [default: 1]

            -min_per_read_length_fraction  MIN_PER_READ_LENGTH_FRACTION
                          min number of consecutive high quality base calls to
                          include a read (per single end read) as a fraction of
                          the input read length [default: 0.5]
  
            -min_otu_size the minimum otu size (in number of sequences) to
                          retain the otu [default: 2]

            -similarity   sequence similarity threshold (for blast, cdhit,
                          uclust, uclust_ref, usearch, usearch_ref, usearch61,
                          or usearch61_ref) [default: 0.94]


            -filter_taxa   Comma-separated list of taxa to discard from OTU table
                           e.g: p__Bacteroidetes,p__Firmicutes

            -substract_NTC   A LIST OF SAMPLE ID separated by comma
                             substarct observation count from No Template Control (NTC)
                             and remove NTC samples from OTU table

            -e            Sequencing_depth_min_cutoff
                          Filter sample less this amount of sequences.
                          The minimium of sequenceing depth of samples after
                          this filter will be  use for even sub-sampling and
                          maximum rarefaction depth. You should review the
                          output of the 'biom summarize-table' command to decide
                          on this value.[default: 1000]

            -c            <INT> # of CPUs to run the script (default:4 )

            -t            <STRING>  Project title

            -debug        <bool> keep intermediate files

# OUTPUT
The major output file is out_directory/analysis/index.html file. It can be open by any browsers. 

# Reference

* UPARSE: http://drive5.com/uparse/
* QIIME: 
    Based on following tutorials
    
    http://qiime.org/tutorials/processing_illumina_data.html
    
    http://nbviewer.jupyter.org/github/biocore/qiime/blob/1.9.1/examples/ipynb/illumina_overview_tutorial.ipynb

    http://nbviewer.jupyter.org/github/biocore/qiime/blob/1.9.1/examples/ipynb/Fungal-ITS-analysis.ipynb

    http://qiime.org/1.4.0/tutorials/processing_18S_data.html
    
    http://qiime.org/tutorials/open_reference_illumina_processing.html
    
    Mapping File format 
    
    http://qiime.org/documentation/file_formats.html

    Data Source

    http://qiime.org/home_static/dataFiles.html
