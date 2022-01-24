# SARS-CoV-2 variant abundance estimation in wastewater

This repository contains all code used for the analysis presented in this
manuscript. We show that by sequencing SARS-CoV-2 RNA in wastewater and applying
computational techniques initially used for RNA-Seq quantification, we can
estimate the abundance of variants in wastewater samples.

This computational pipeline can easily be reused for other wastewater samples.
The exact scripts and commands used to analyze the data discussed in our
manuscript can be found in the manuscript folder.


## Building a reference set
We build a reference set by selecting representative genomes per lineage from
the GISAID database. This selection can be based on country or region; for
example, in our analysis we only consider sequences from US origin. Because the
GISAID database is very large, we apply a quality filter and subsequently
call variants and compute allele frequencies from 1000 sequences per lineage.
Required input data for this step is a fasta file with GISAID sequences and a
matching metadata tsv file. Both the sequences (.fasta) and the metadata (.tsv)
are publicly available on GISAID but *only after registration*.

First we apply a quality filter while preprocessing GISAID sequences:

    python pipeline/preprocess_references.py -m <GISAID_metadata.tsv> -f <GISAID_sequences.fasta> -k 1000 --seed 0 --country USA --o reference_set   

Note that this preprocessing script also allows you to specify a specific state
to restrict the reference set to, by using the `--state` option.
Similarly, you can restrict the sequences used by collection date using
`--startdate` and `--enddate`.
These dates must be provided in ISO format (Y-M-D).

Then we call variants compared to the original SARS-CoV-2 reference (here we use
[NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512))
and compute allele frequencies per lineage. Note that this requires
`pyvcf`, `vcftools`, `bcftools`, `samtools` and `minimap2` to be installed, which
can all be installed through [bioconda](http://bioconda.github.io).
We run the following script to call variants:

    pipeline/call_variants.sh reference_set <full_path_to_main_ref_fasta>

Please make sure to provide the full path to the main reference sequence
(for example `/Users/username/wastewater_analysis/references/ref.fa`).

Based on the resulting allele frequencies, we select sequences per lineage such that all
mutations with an allele frequency of at least 50% were captured at least once.

    python pipeline/select_samples.py -m <GISAID_metadata.tsv> -f <GISAID_sequences.fasta> -o reference_set --vcf reference_set/*_merged.vcf.gz --freq reference_set/*_merged.frq

The resulting reference set is stored in `reference_set/sequences.fasta` and
the corresponding metadata in `reference_set/metadata.tsv`.

### Pre-selected reference set
GISAID sequence identifiers for the reference set used in our manuscript are
provided in `auxiliary_data/reference_set_03_2021.txt`.


## Preprocessing sequencing data
Before processing with kallisto, we recommend removing adapter sequences from
the reads using Trimmomatic, identifying primer sequences with iVar and then
trimming the primers off with [jvarkit](http://lindenb.github.io/jvarkit/Biostar84452).
Trimmomatic and iVar can be installed through [bioconda](http://bioconda.github.io)
and jvarkit can be downloaded [here](http://lindenb.github.io/jvarkit/Biostar84452).

All preprocessing commands are listed in `pipeline/trim_reads.sh` and can be
executed at once by running this script. The commands should be adjusted to
match the file locations locally and require adapter sequences (.fasta) and
primer sequences (.bed) for trimming.

## Predicting variant abundance
Now we can predict abundance per lineage using
[kallisto](https://pachterlab.github.io/kallisto/about), which can be
installed through [bioconda](http://bioconda.github.io).
First we build the kallisto index for our reference set:

    kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta

Then we calculate abundance per reference sequence:

    kallisto quant -t 20 -b <num_bootstraps> -i reference_set/sequences.kallisto_idx -o <outdir> -t 20 <forward.fastq> <reverse.fastq>

Finally, the kallisto output can be processed to obtain variant abundance estimates:

    python pipeline/output_abundances.py -m <min_ab> -o <outdir>/predictions.tsv --metadata reference_set/metadata.tsv --voc B.1.1.7,B.1.351,B.1.427,B.1.429,B.1.526,P.1 <outdir>/abundance.tsv

In this command, lineages of interest can be provided as a comma-separated list with
`--voc` (see above).
Alternatively, variants can be defined as combinations of lineages in a json file
(see auxiliary_data/WHO_variants_2021-10-06.json for an example) by using `--voc_file`.
If neither `--voc` nor `--voc_file` are provided, abundance predictions for ALL lineages
in the reference set are written to the output file.

The predictions are written to `<outdir>/predictions.tsv`, unless specified otherwise with `-o`.


## Example

The `example` directory contains a small example to test the prediction step.
A toy example of reference sequences is given in `sequences.fa` with metadata in
`metadata.tsv`.

    cd example
    python ../pipeline/preprocess_references.py -m metadata.tsv -f sequences.fa -k 10 --seed 0 -o reference_set
    bash ../pipeline/call_variants.sh reference_set "${PWD}/SARS-CoV-2-NC_045513.fasta"
    python ../pipeline/select_samples.py -m metadata.tsv -f sequences.fa -o reference_set --vcf reference_set/*_merged.vcf.gz --freq reference_set/*_merged.frq
    kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta
    kallisto quant -t 20 -b 0 -i reference_set/sequences.kallisto_idx -o kallisto -t 1 input_1.fastq input_2.fastq
    python ../pipeline/output_abundances.py --metadata metadata.tsv -o kallisto/predictions.tsv kallisto/abundance.tsv

The predictions can be found in `kallisto/predictions.tsv` and should show that
100% of this sample is SARS-CoV-2.

*Note that this is just a toy example, do NOT use this reference set for any of
datasets that you actually want to analyze. Download the GISAID data and build
a good reference set first.*
