# Manuscript

We show by sequencing samples from wastewater and clinical
isolates in Connecticut between January and April 2021 that the temporal
dynamics of variant strains broadly correspond. We further show that this
technique can be used with other wastewater sequencing techniques by expanding
to samples taken across the United States in a similar timeframe. Below,
we describe step by step how this analysis was done.

*Note that these scripts are slightly different from the most recent scripts in the pipeline folder, because the GISAID metadata formatting has changed since. The scripts in this folder (manuscript) are applicable only to the GISAID metadata formatting from early 2021.*

## US-specific reference set
Since all samples are taken from locations in the US, we work with a US-specific
reference set. The reference set is constructed from the GISAID database
(downloaded 9 March 2021) using the following commands.

We begin by counting the number of non-ambiguous
nucleotides per sequence:

    sed '/^>/d' sequences_2021-03-04_08-34.fasta | tr -d 'N' | awk '{ print length; }' > sequences_2021-03-04_08-34.nonN_chars.txt

Next, we preprocess references, call variants, select references and build a kallisto index:

    python manuscript/preprocess_references_v1.py -m GISAID/metadata_2021-03-04_10-31.tsv -f GISAID/sequences_2021-03-04_08-34.fasta -k 1000 --seed 0 --country USA --min_len 29500 -o reference_set -n GISAID/sequences_2021-03-04_08-34.nonN_chars_per_id.txt
    sbatch pipeline/call_variants.sh reference_set
    python manuscript/select_samples_v1.py -m GISAID/metadata_2021-03-04_10-31.tsv -f GISAID/sequences_2021-03-04_08-34.fasta -o reference_set -n GISAID/sequences_2021-03-04_08-34.nonN_chars_per_id.txt --vcf reference_set/*_merged.vcf.gz --freq reference_set/*_merged.frq
    kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta


## Benchmarking experiments
To evaluate the accuracy of the predictions obtained through our pipeline, we
created a collection of benchmarking datasets based on GISAID sequences in
Connecticut on 2021-03-04.

    mkdir benchmarks
    python benchmarking/create_benchmarks.py --voc_perc 0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 -m GISAID/metadata_2021-03-04_10-31.tsv -fr GISAID/sequences_2021-03-04_08-34.fasta -s Connecticut -d 2021-02-11 -fv GISAID/B.1.1.7_EPI_ISL_1064784.fasta,GISAID/P.1_EPI_ISL_1239974.fasta,GISAID/B.1.351_EPI_ISL_1038809.fasta,GISAID/B.1.427_EPI_ISL_755182.fasta,GISAID/B.1.429_EPI_ISL_1063907.fasta -o benchmarks/Connecticut-WG-2021-02-11-100x --total_cov 100
    python benchmarking/create_benchmarks.py --voc_perc 0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 -m GISAID/metadata_2021-03-04_10-31.tsv -fr GISAID/sequences_2021-03-04_08-34.fasta -s Connecticut -d 2021-02-11 -fv GISAID/B.1.1.7_EPI_ISL_1064784.fasta,GISAID/P.1_EPI_ISL_1239974.fasta,GISAID/B.1.351_EPI_ISL_1038809.fasta,GISAID/B.1.427_EPI_ISL_755182.fasta,GISAID/B.1.429_EPI_ISL_1063907.fasta -o benchmarks/Connecticut-WG-2021-02-11-1000x --total_cov 1000
    python benchmarking/create_benchmarks.py --voc_perc 0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 -m GISAID/metadata_2021-03-04_10-31.tsv -fr GISAID/sequences_2021-03-04_08-34.fasta -s Connecticut -d 2021-02-11 -fv GISAID/B.1.1.7_EPI_ISL_1064784.fasta,GISAID/P.1_EPI_ISL_1239974.fasta,GISAID/B.1.351_EPI_ISL_1038809.fasta,GISAID/B.1.427_EPI_ISL_755182.fasta,GISAID/B.1.429_EPI_ISL_1063907.fasta -o benchmarks/Connecticut-2021-02-11-100x --total_cov 100 --spike_only
    python benchmarking/create_benchmarks.py --voc_perc 0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 -m GISAID/metadata_2021-03-04_10-31.tsv -fr GISAID/sequences_2021-03-04_08-34.fasta -s Connecticut -d 2021-02-11 -fv GISAID/B.1.1.7_EPI_ISL_1064784.fasta,GISAID/P.1_EPI_ISL_1239974.fasta,GISAID/B.1.351_EPI_ISL_1038809.fasta,GISAID/B.1.427_EPI_ISL_755182.fasta,GISAID/B.1.429_EPI_ISL_1063907.fasta -o benchmarks/Connecticut-2021-02-11-1000x --total_cov 1000 --spike_only
    python benchmarking/create_benchmarks.py --voc_perc 0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100 -m GISAID/metadata_2021-03-04_10-31.tsv -fr GISAID/sequences_2021-03-04_08-34.fasta -s Connecticut -d 2021-02-11 -fv GISAID/B.1.1.7_EPI_ISL_1064784.fasta,GISAID/P.1_EPI_ISL_1239974.fasta,GISAID/B.1.351_EPI_ISL_1038809.fasta,GISAID/B.1.427_EPI_ISL_755182.fasta,GISAID/B.1.429_EPI_ISL_1063907.fasta -o benchmarks/Connecticut-2021-02-11-10000x --total_cov 10000 --spike_only

Then, we run kallisto on each of these benchmarks:

    bootstraps=100
    min_ab=0.1
    ref_dir=reference_set
    for dataset in Connecticut-WG-2021-02-11-100x Connecticut-WG-2021-02-11-1000x Connecticut-2021-02-11-100x Connecticut-2021-02-11-1000x Connecticut-2021-02-11-10000x; do \
        sbatch manuscript/run_kallisto_ref_sets.sh ${dataset} ${bootstraps} ${ref_dir} ${min_ab};
    done


## Experiments on real sequencing data
We preprocess sequencing data to remove adapters and primers from the reads:

    sbatch pipeline/trim_reads.sh ginkgo_1_fastqs analysis_july_5 ginkgo_1_ids.txt auxiliary_data/adapters.fa auxiliary_data/primers.bed MN908947.3.DNA.fasta
    sbatch pipeline/trim_reads.sh ww_fastqs_2021-05-25 yale/yale-batch3 yale/yale-batch3/sample_ids.txt auxiliary_data/adapters.fa auxiliary_data/primers.bed MN908947.3.DNA.fasta

Then, we predict variant abundance from the wastewater sequencing data:

    sbatch manuscript/run_kallisto.sh new_haven_ids.txt reference_set yale/yale-batch3 kallisto_refs_USA 100
    sbatch manuscript/run_kallisto.sh ginkgo_1_ids.txt reference_set biobot/analysis_july_5 kallisto_refs_USA 100


## Analysis
Finally, we analyze predictions by comparing to clinical frequencies or GISAID
frequencies. All figures presented in the manuscript are generated by running
`make_figures.sh <OUTDIR> <DATADIR> <HOMEDIR>` where OUTDIR is the directory
where all plots go, DATADIR is the directory with the sequencing data, and
HOMEDIR is the directory with the benchmarking data.
