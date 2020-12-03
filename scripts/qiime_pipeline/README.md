This is an amplicon analysis pipeline based on Qiime2 developed by Genome Science group of Los Alamos National Laboratory.

# REQUIREMENTS

* qiime2 version 2019.1 (http://qiime2.org/)
* python xlsx2csv (optional for EXCEL xlsx format mapping file support, https://github.com/dilshod/xlsx2csv) 

# USAGE 

```
usage: qiime_pipeline.py [-h] -m MAPPINGFILE
                         (-p [FASTQ] [FASTQ] | -s [FASTQ] | -d [PATH])
                         [-b [FASTQ]] [-o [PATH]] [--zip] [--target <STR>]
                         [--qcMethod <STR>] [--trimLeftForward <INT>]
                         [--trimLeftReverse <INT>] [--truncLenForward <INT>]
                         [--truncLenReverse <INT>] [--trimLen <INT>]
                         [--truncLen <INT>] [--minQuality <INT>]
                         [--minLengthFraction <FLOAT>] [--maxAmbiguous <INT>]
                         [--samplingDepth <INT>] [--autoDepth]
                         [--maxRarefactionDepth <INT>] [--phred_offset <INT>]
                         [-c <INT>] [--title <STR>] [--version]

Script to run Qiime2 pipeline based on Moving Pictures tutorial

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Mapping File (required):
  -m MAPPINGFILE, --mappingFile MAPPINGFILE
                        MAPPING file. #SampleID BarcodeSequence
                        LinkerPrimerSequence sampleType ... Description

Reads Input (required, mutually exclusive):
  -p [FASTQ] [FASTQ], --paired [FASTQ] [FASTQ]
                        Paired reads in two fastq files and separate by space
  -s [FASTQ], --single [FASTQ]
                        Unpaired reads fastq
  -d [PATH], --dir [PATH]
                        Contains Demultiplexed fastq files. To use this
                        option. To use this option,the mapping file need a
                        'Files' column with filenames for each sampleID.

Barcode Fastq (requried when -p or -s):
  -b [FASTQ], --barcode [FASTQ]
                        Barcodes fastq

Output:
  -o [PATH], --outdir [PATH]
                        Output directory
  --zip                 Zip output files for visualization

Parameters:
  --target <STR>        Greengenes, SILVA, SILVA-V3-V4 or ITS. [default:
                        Greengenes] Greengenes and SILVA-V3-V4 are for 16s.
                        SILVA is for 16s/18s. ITS is from https://unite.ut.ee/
                        and for fungal rDNA ITS sequences. ## the trained qza
                        files in data/
  --qcMethod <STR>      Quality control method. dada2 or deblur [default:
                        dada2]
  --trimLeftForward <INT>
                        This is for Dada2 QC on paired end reads. The number
                        of bases to trim from the left of the forward (R1)
                        read to e.g. remove PCR primer sequences. [default:
                        20]
  --trimLeftReverse <INT>
                        This is for Dada2 QC on paired end reads. This is the
                        number of bases to trim from the left of the reverse
                        (R2) read to e.g. remove PCR primer sequences.
                        [default: 20]
  --truncLenForward <INT>
                        This is for Dada2 QC on paired end reads. This is the
                        truncation length of the forward reads after any
                        trimming. "0" is no truncation. [default: 0]
  --truncLenReverse <INT>
                        This is for Dada2 QC on paired end reads. This is the
                        truncation length of the reverse reads after any
                        trimming. "0" is no truncation. [default: 0]
  --trimLen <INT>       This is for Dada2 QC on single end reads. This is the
                        number of bases to trim from the reads to e.g. remove
                        PCR primer sequences.[default: 20]
  --truncLen <INT>      This works for Dada2 and Deblur on single end reads to
                        truncate sequences at position [default: 0] no
                        truncate
  --minQuality <INT>    This is for Deblur QC. The minimum acceptable PHRED
                        score. All PHRED scores less that this value are
                        considered to be low PHRED scores. [default: 4]
  --minLengthFraction <FLOAT>
                        This is for Deblur QC. The minimum length that a
                        sequence read can be following truncation and still be
                        retained. This length should be provided as a fraction
                        of the input sequence length. [default: 0.75]
  --maxAmbiguous <INT>  This is for Deblur QC. The maximum number of ambiguous
                        (i.e., N) base calls. This is applied after trimming
                        sequences based on `min_length_fraction. [default: 0]
  --samplingDepth <INT>
                        Filter sample less this amount of sequences.The
                        minimium of sequenceing depth of samples after this
                        filter will be use for even sub-sampling and maximum
                        rarefaction depth.
  --autoDepth           Automatically adjust the sampling to the minimum
                        sequences count of all samples. The minimum >
                        samplingDepth option above.
  --maxRarefactionDepth <INT>
                        The maximum rarefaction depth.[default: same as
                        samplingDepth]
  --phred_offset <INT>  The ascii offset to use when decoding phred scores
  -c <INT>, --cpus <INT>
                        Number of CPUS
  --title <STR>         Project Title

```

# OUTPUT
The major output file is out_directory/index.html file. It can be open by any browsers. 

# Reference

* Based on following tutorials
    
    https://docs.qiime2.org/2019.1/tutorials/moving-pictures/
    
    https://docs.qiime2.org/2019.1/tutorials/overview/

    https://docs.qiime2.org/2019.1/tutorials/

    https://docs.qiime2.org/2019.1/tutorials/feature-classifier/

* Mapping File format 
    
    https://docs.qiime2.org/2019.1/tutorials/metadata/

* Data Source

    https://docs.qiime2.org/2019.1/data-resources/
