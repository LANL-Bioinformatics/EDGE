#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-11:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH -o logs/run_salmon_%j.out
#SBATCH -e logs/run_salmon_%j.err

dataset=$1

mkdir -p salmon/benchmarks/${dataset}

for VOC in P.1_EPI_ISL_1239974 B.1.1.7_EPI_ISL_1064784 B.1.351_EPI_ISL_1038809 B.1.427_EPI_ISL_755182 B.1.429_EPI_ISL_1063907; do \
  for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do \
    salmon quant -p 20 -i GISAID/downloads/ref_selection/sequences.salmon_idx -o salmon/benchmarks/${dataset}/${VOC}_ab${ab} -l A -1 benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_1.fastq -2 benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_2.fastq;
  done;
done;
