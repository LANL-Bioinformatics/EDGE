#!/usr/bin/env bash

../runPipeline -p Ecoli_10x.1.fastq Ecoli_10x.2.fastq -c $PWD/config.txt -o $PWD/output -ref $PWD/Reference/NC_000913.gbk -cpu 4 -primer $PWD/primers.fa -noColorLog
