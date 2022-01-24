#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-5:00
#SBATCH -p short
#SBATCH --mem=20G
#SBATCH -o logs/run_kallisto_%j.out
#SBATCH -e logs/run_kallisto_%j.err

dataset=$1
num_bootstraps=$2
ref_dir=$3
min_ab=$4

HDF5_USE_FILE_LOCKING=FALSE # prevent HDF5 problems (https://github.com/pachterlab/kallisto/issues/197)

outdir=kallisto/benchmarks/${dataset}_${ref_dir}
mkdir -p ${outdir}

for VOC in P.1_EPI_ISL_1239974 B.1.1.7_EPI_ISL_1064784 B.1.351_EPI_ISL_1038809 B.1.427_EPI_ISL_755182 B.1.429_EPI_ISL_1063907; do \
  for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do \
    /usr/bin/time -v kallisto quant -t 20 -b ${num_bootstraps} -i ${ref_dir}/sequences.kallisto_idx -o ${outdir}/${VOC}_ab${ab} benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_1.fastq benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_2.fastq > ${outdir}/${VOC}_ab${ab}.log 2>&1;
    /usr/bin/time -v python pipeline/output_abundances.py -m $min_ab -o ${outdir}/${VOC}_ab${ab}/predictions_m${min_ab}.tsv --metadata GISAID/downloads/${ref_dir}/metadata.tsv --voc B.1.1.7,B.1.351,B.1.427,B.1.429,P.1 ${outdir}/${VOC}_ab${ab}/abundance.tsv >> ${outdir}/${VOC}_ab${ab}.log 2>&1;
  done;
done;
