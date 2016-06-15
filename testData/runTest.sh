#!/usr/bin/env bash

../runPipeline -p Ecoli_10x.1.fastq Ecoli_10x.2.fastq -c $PWD/config.txt -o $PWD/output -ref $PWD/Reference/NC_000913.fna -cpu 10 -primer $PWD/primers.fa -noColorLog
