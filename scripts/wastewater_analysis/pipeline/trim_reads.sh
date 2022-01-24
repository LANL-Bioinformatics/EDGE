#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH -o logs/trim_reads_%j.out
#SBATCH -e logs/trim_reads_%j.err

fastqdir=$1
outdir=$2
sample_ids=$3
adapters=$4
primers=$5
ref=$6

mkdir -p $outdir;
mkdir -p $outdir/trimmomatic;
mkdir -p $outdir/ivar;

while read sample; do
    echo $sample;
    # trim adapters with trimmomatic
    /usr/bin/time -v trimmomatic PE -threads 20 ${fastqdir}/*${sample}_*_001.fastq.gz ${outdir}/trimmomatic/${sample}_1.fastq ${outdir}/trimmomatic/${sample}_s1.fastq ${outdir}/trimmomatic/${sample}_2.fastq ${outdir}/trimmomatic/${sample}_s2.fastq ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads SLIDINGWINDOW:4:15 MINLEN:36 LEADING:3 TRAILING:3 > $outdir/trimmomatic/${sample}_trimmomatic.log 2>&1;
    # align reads with bwa mem
    /usr/bin/time -v bwa mem -L 100 -t 20 $ref ${outdir}/trimmomatic/${sample}_1.fastq ${outdir}/trimmomatic/${sample}_2.fastq | samtools view -bh | samtools sort > ${outdir}/trimmomatic/${sample}.bam;
    samtools index ${outdir}/trimmomatic/${sample}.bam;
    # trim primers with ivar (soft clipping)
    /usr/bin/time -v ivar trim -i ${outdir}/trimmomatic/${sample}.bam -b $primers -p ${outdir}/ivar/${sample} -e > ${outdir}/ivar/${sample}_ivar.log 2>&1;
    # remove soft-clipped primers
    java -jar ~/jvarkit/dist/biostar84452.jar --samoutputformat BAM <(samtools sort ${outdir}/ivar/${sample}.bam) > ${outdir}/ivar/${sample}.trimmed.bam;
    # extract fastqs
    samtools fastq -1 ${outdir}/ivar/${sample}.forward.fastq -2 ${outdir}/ivar/${sample}.reverse.fastq -s ${outdir}/ivar/${sample}.singles.fastq <(samtools sort -n ${outdir}/ivar/${sample}.trimmed.bam);
    # # extract reads aligning to spike (21063-25884)
    # samtools index ${outdir}/ivar/${sample}.trimmed.bam;
    # samtools view -bh ${outdir}/ivar/${sample}.trimmed.bam "MN908947.3:21063-25884" > ${outdir}/ivar/${sample}.trimmed.spike.bam
    # samtools fastq -1 ${outdir}/ivar/${sample}_spike.forward.fastq -2 ${outdir}/ivar/${sample}_spike.reverse.fastq -s ${outdir}/ivar/${sample}_spike.singles.fastq <(samtools sort -n ${outdir}/ivar/${sample}.trimmed.spike.bam);
done < ${sample_ids}
