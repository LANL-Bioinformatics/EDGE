#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-2:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH -o call_variants_%j.out
#SBATCH -e call_variants_%j.err

ref_dir=$1 # directory with sequences for variant calling
single_ref_path=$2 # path to single reference sequence to compare against

cd ${ref_dir}

while read lineage; do \
    cd $lineage;
    for fasta in *.fa; do \
        # align and sort
        minimap2 -c -x asm20 --end-bonus 100 -t 20 --cs $single_ref_path $fasta 2>${fasta%.fa}.paftools.log | sort -k6,6 -k8,8n > ${fasta%.fa}.paf && paftools.js call -s ${fasta%.fa} -L 100 -f $single_ref_path ${fasta%.fa}.paf > ${fasta%.fa}.vcf 2>>${fasta%.fa}.paftools.log;
        bgzip -f ${fasta%.fa}.vcf;
        bcftools index -f ${fasta%.fa}.vcf.gz;
    done;
    cd ..;
    sample_count=$(ls ${lineage}/*.vcf.gz | wc -l);
    if [[ ${sample_count} -eq 1 ]]; then \
        cp ${lineage}/*.vcf.gz ${lineage}_merged.vcf.gz;
    else \
        bcftools merge -o ${lineage}_merged.vcf.gz -O z -0 ${lineage}/*.vcf.gz;
    fi;
    vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --site-pi;
    vcftools --gzvcf ${lineage}_merged.vcf.gz --out ${lineage}_merged --freq;
done < lineages.txt

cd ..
