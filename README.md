# EDGE Bioinformatics

EDGE COVID-19 is a bioinformatics platform that provides standardized workflows for genome ‘assembly’ and preliminary analysis of Illumina or Nanopore data for SARS-CoV-2 genomes. The basic workflow includes data quality control and filtering, alignment of reads to a reference genome ([NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/)), creation of a consensus genome sequence based on the read alignments, and a preliminary Single Nucleotide Polymorphism analysis. 

The platform contains read alignment tools that accommodate Illumina or Nanopore data, including those generated using the SARS-CoV-2 ARTIC network sequencing workflows (https://artic.network/ncov-2019). It consists of a user-friendly GUI that drives a series of best-practice, open source bioinformatics software to aid in the reconstruction of complete genomes of SARS-CoV-2 from whole genome sequencing projects.

This light-weight version is available as a [Docker container](https://hub.docker.com/r/bioedge/edge_ncov), able to run on any local hardware infrastructure or in the cloud, and is based on the original [EDGE Bioinformatics platform](https://edgebioinformatics.org).  At this moment, no reference sequence database is required, and it should be able to run on a laptop. 

For users who want to use read- or assembly-based taxonomy classification tools to understand what organisms may be present within samples, we refer you to the original [EDGE Bioinformatics platform](https://edgebioinformatics.org) which harbors this workflow but requires a number of large databases to enable such a search. In initial tests, we recover Bat Coronavirus XXX as the genome in the databases as the closest near-neighbor to the SARS-CoV-2 genome sequence, but are actively working on generating new reference indexes to also allow hits to recent SARS-CoV-2 genomes. 

Users can input/upload Illumina or Nanopore sequencing fastq files (and/or download from NCBI SRA), perform *de novo* assembly, map reads to a reference genome to get the consensus SARS-CoV-2 sequence from the sequencing data. Users can also compare the assembled contigs to a reference sequence (default to NCBI refseq sequence [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/) and see the differences in both tabular and graphic form. 

## Documentation
	
[EDGE COVID-19 Guide](https://docs.google.com/document/d/e/2PACX-1vSNJntmg9Cc7dm1Agh0gXyp9LZCGQfNyiBgAiIYCY8CN8Lk2Ma8tBwEBhIut67ow5ItDuXtUsi0v2Du/pub)
[EDGE Docker Hub](https://hub.docker.com/r/bioedge/edge_ncov)

## Requirements

* Docker Engine version 19.03.2 or greater
* The image size is around 11.8GB. This image is built on top of the official Ubuntu 18.04.3 LTS Base Image.
* Recommended minimum computational resource : 8G memory, 4 CPUs, 20GB for the image.   

## Contact Info
Support: The EDGE user's group is a Google group for users and developers: [edge-users@googlegroups.com](mailto:edge-users@googlegroups.com)

## Citation

Po-E Li, Chien-Chi Lo, Joseph J. Anderson, Karen W. Davenport, Kimberly A. Bishop-Lilly, Yan Xu, Sanaa Ahmed, Shihai Feng, Vishwesh P. Mokashi, Patrick S.G. Chain; Enabling the democratization of the genomics revolution with a fully integrated web-based bioinformatics platform, Nucleic Acids Research, Volume 45, Issue 1, 9 January 2017, Pages 67–80, https://doi.org/10.1093/nar/gkw1027

## Copyright

Copyright (2020).  Triad National Security, LLC. All rights reserved.
 
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National 
Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National 
Nuclear Security Administration.
 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National 
Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, 
paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to 
the public, perform publicly and display publicly, and to permit others to do so.

This is open source software; you can redistribute it and/or modify it under the terms of the GPLv3 License. If software 
is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with 
the version available from LANL. Full text of the [GPLv3 License](https://github.com/LANL-Bioinformatics/edge/blob/master/LICENSE) can be found in the License file in the main development 
branch of the repository.
