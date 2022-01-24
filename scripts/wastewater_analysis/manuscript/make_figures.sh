#!/bin/bash

OUTDIR=$1
DATADIR=$2
HOMEDIR=$3

mkdir -p $OUTDIR

################### MAIN ###################

# figure 1 ---> designed in Inkscape

# figure 2
for dataset in Connecticut-2021-02-11-10000x Connecticut-WG-2021-02-11-1000x; do \
    for m in 0.1 1; do \
        python benchmarking/evaluate_abundances.py \
            -s _${dataset}_m${m} \
            -m $m \
            --voc B.1.1.7,B.1.351,B.1.427,B.1.429,P.1 \
            -o ${OUTDIR} \
            -v \
            --output_format png,svg \
            ${HOMEDIR}/kallisto/benchmarks/${dataset}_random1000_USA/*/predictions_m${m}.tsv;
    done;
done

# figure 3 panel a,b
python analysis/plot_RNA_levels.py \
        --case_rate_info ${DATADIR}/yale/raw_case_and_SARS_data_New_Haven_jan-april2021.csv \
        --run_info ${DATADIR}/yale/yale-batch3/kallisto_refs_USA/*/run_info.json \
        --metadata ${DATADIR}/yale/yale-batch3/metadata_full.csv \
        -o $OUTDIR \
        --aln_stats ${DATADIR}/yale/yale-batch3/ivar/*.stats\
        --outsuffix _yale \
        --outformat png,svg

# figure 3 panel c
python analysis/plot_results_cov_exp.py \
        --predictions ${DATADIR}/yale/cov_exp_ER2/p*/S*/predictions_m0.1.tsv \
        --voc B.1.1.7 \
        --repeats 100 \
        -o ${OUTDIR}/fig1c.png

# figure 4
for voc in B.1.1.7 B.1.351 B.1.427 B.1.429 B.1.526 P.1; do \
    python analysis/plot_results_clinical.py \
        --voc $voc \
        --predictions ${DATADIR}/yale/yale-batch3/plots/predictions_${voc}_m0.1_a500000.csv \
        --clinical_data ${DATADIR}/yale/NHV_lineages_by_week_2021Jan_23May.csv \
        --outdir $OUTDIR \
        --fig_width=15 \
        --outprefix yale_predictions_clinical \
        --output_format png,svg;
done

# figure 5
for voc in B.1.1.7 B.1.351 B.1.427 B.1.429 B.1.526 P.1; do \
    python analysis/plot_results_biobot.py \
        --run_info ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/run_info.json \
        --confidence 95 \
        --metadata ${DATADIR}/biobot/metadata_full.csv \
        --metadata_ref ${HOMEDIR}/GISAID/downloads/random1000_USA/metadata.tsv \
        --voc $voc \
        --min_ab 0.1 \
        --outdir $OUTDIR \
        --outprefix biobot_predictions_m0.1_a500000 \
        --min_aln 500000 \
        --output_format png,svg \
        --fig_width 12 \
        --gisaid_freqs ${DATADIR}/biobot/frequencies_gisaid.tsv \
        --sep_states ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/predictions_m0.1.tsv;
done

################ SUPPLEMENT ################

# figure S1
# "Within-lineage diversity observed in SARS-CoV-2 genomes on GISAID"
python analysis/plot_sequence_diversity.py \
        --outdir $OUTDIR \
        --voc_names B.1.1.7 B.1.351 B.1.427 B.1.429 B.1.526 P.1 \
        --site_pi_files ${HOMEDIR}/GISAID/downloads/random1000_USA/{B.1.1.7,B.1.351,B.1.427,B.1.429,B.1.526,P.1}_merged.sites.pi \
        --ref_size 29903 \
        --allele_freq_files ${HOMEDIR}/GISAID/downloads/random1000_USA/{B.1.1.7,B.1.351,B.1.427,B.1.429,B.1.526,P.1}_merged.frq \
        --min_af 0.5

# figure S2
# background noise in predictions
python analysis/plot_raw_predictions.py \
        ${HOMEDIR}/kallisto/benchmarks/Connecticut-WG-2021-02-11-100x_random1000_USA/B.1.1.7_EPI_ISL_1064784_ab100/abundance.tsv \
        --outprefix B.1.1.7_ab100_ \
        --outformat png,svg \
        --outdir $OUTDIR

# figure S3
# "additional benchmarking results"
for dataset in Connecticut-2021-02-11-100x Connecticut-2021-02-11-1000x Connecticut-WG-2021-02-11-100x; do \
    for m in 0.1 1; do \
        python benchmarking/evaluate_abundances.py \
            -s _${dataset}_m${m} \
            -m $m \
            --voc B.1.1.7,B.1.351,B.1.427,B.1.429,P.1 \
            -o ${OUTDIR} \
            -v \
            --output_format png,svg \
            ${HOMEDIR}/kallisto/benchmarks/${dataset}_random1000_USA/*/predictions_m${m}.tsv;
    done;
done

# figure S4
# "salmon versus kallisto abundance estimates per VOC"
python analysis/plot_kallisto_vs_salmon.py \
        --kallisto ${HOMEDIR}/kallisto/benchmarks/Connecticut-WG-2021-02-11-1000x/*/predictions_m0.1.tsv \
        --salmon ${HOMEDIR}/salmon/benchmarks/Connecticut-WG-2021-02-11-1000x/*/predictions_m0.1.tsv \
        --voc B.1.1.7,B.1.351,B.1.427,B.1.429,P.1 \
        -o $OUTDIR \
        -m 0.1 \
        --output_format png,svg

# figure S5
# "sequencing depth along the genome obtained from sludge samples in New Haven, CT."
python analysis/plot_coverage.py \
        --samples "FH1 (20%),EX1 (63%),ER2 (99%)" \
        -o ${OUTDIR}/FH1_EX1_ER2_depth.png \
        ${DATADIR}/yale/yale-batch3/ivar/FH1.trimmed.depth \
        ${DATADIR}/yale/yale-batch3/ivar/EX1.trimmed.depth \
        ${DATADIR}/yale/yale-batch3/ivar/ER2.trimmed.depth;
python analysis/plot_coverage.py \
        --samples "FH1 (20%),EX1 (63%),ER2 (99%)" \
        -o ${OUTDIR}/FH1_EX1_ER2_depth.svg \
        ${DATADIR}/yale/yale-batch3/ivar/FH1.trimmed.depth \
        ${DATADIR}/yale/yale-batch3/ivar/EX1.trimmed.depth \
        ${DATADIR}/yale/yale-batch3/ivar/ER2.trimmed.depth;
# python analysis/plot_coverage.py -o ${OUTDIR}/depth_all.svg ${DATADIR}/yale/yale-batch3/ivar/*.trimmed.depth

# figure S6 (subplots)
# "predictions per VOC with confidence intervals based on bootstrap analysis for New Haven samples"
for voc in B.1.1.7 B.1.351 B.1.427 B.1.429 B.1.526 P.1; do \
    python analysis/plot_bootstrap_results.py \
        --bootstraps ${DATADIR}/yale/yale-batch3/kallisto_bootstrap_refs_USA/*/abundance.h5 \
        --confidence 95 \
        --metadata ${DATADIR}/yale/yale-batch3/metadata_full.csv \
        --metadata_ref ${HOMEDIR}/GISAID/downloads/random1000_USA/metadata.tsv \
        --voc $voc \
        --min_ab 0 \
        --outdir $OUTDIR \
        --outprefix ${DATADIR}/yale_bootstrap_m0_a500000 \
        --output_format png \
        --run_info ${DATADIR}/yale/yale-batch3/kallisto_bootstrap_refs_USA/*/run_info.json \
        --min_aln 500000 \
        --fig_width 15 \
        ${DATADIR}/yale/yale-batch3/kallisto_bootstrap_refs_USA/*/predictions_m0.tsv;
done

# figure S7
# "percent genome coverage versus Ct values for samples across the US"
python analysis/plot_RNA_levels.py \
        --run_info ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/run_info.json \
        --metadata ${DATADIR}/biobot/metadata_full.csv \
        --aln_stats ${DATADIR}/biobot/analysis_july_5/ivar/*.stats \
        -o $OUTDIR \
        --outsuffix _biobot \
        --outformat png,svg

# figure S8 (subplots)
# "predictions per VOC with confidence intervals based on bootstrap analysis for samples across the US"
for voc in B.1.1.7 B.1.351 B.1.427 B.1.429 B.1.526 P.1; do \
    python analysis/plot_bootstrap_results.py \
        --bootstraps ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/abundance.h5 \
        --run_info ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/run_info.json \
        --confidence 95 \
        --metadata ${DATADIR}/biobot/metadata_full.csv \
        --metadata_ref ${HOMEDIR}/GISAID/downloads/random1000_USA/metadata.tsv \
        --voc $voc \
        --min_ab 0 \
        --outdir $OUTDIR \
        --outprefix biobot_bootstrap_m0_a500000 \
        --min_aln 500000 \
        --output_format svg,png \
        --fig_width 15 \
        --biobot \
        --sep_states \
        ${DATADIR}/biobot/analysis_july_5/kallisto_refs_USA/*/predictions_m0.tsv;
done
