#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-11:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH -o logs/run_kallisto_%j.out
#SBATCH -e logs/run_kallisto_%j.err

sample_ids=$1
ref_dir=$2
batch_dir=$3
kallisto_dir=$4
num_bootstraps=$5

min_ab=0

mkdir -p ${batch_dir}/${kallisto_dir};

while read sample; do
    echo $sample;
    outdir=${batch_dir}/${kallisto_dir}/${sample};
    /usr/bin/time -v kallisto quant -t 20 -b $num_bootstraps -i ${ref_dir}/sequences.kallisto_idx -o ${outdir} -t 20 ${batch_dir}/ivar/${sample}.forward.fastq ${batch_dir}/ivar/${sample}.reverse.fastq > ${batch_dir}/${kallisto_dir}/${sample}.log 2>&1;
    python pipeline/output_abundances.py -m ${min_ab} -o ${outdir}/predictions_m${min_ab}.tsv --metadata ${ref_dir}/metadata.tsv --voc B.1.1.7,B.1.351,B.1.427,B.1.429,B.1.526,P.1 ${outdir}/abundance.tsv;
done < ${sample_ids}
