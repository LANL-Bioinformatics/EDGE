# EDGE COVID-19

[EDGE COVID-19](https://edge-covid19.edgebioinformatics.org/) is a tailored bioinformatics platform based on the more flexible and fully open-source <a href="https://edgebioinformatics.org" target="_new">EDGE Bioinformatics</a> software (<a href="https://doi.org/10.1093/nar/gkw1027" target="_new">Li et al. 2017</a>). This mini-version consists of a user-friendly GUI that drives standardized workflows for genome reference-based 'assembly' and preliminary analysis of Illumina or Nanopore data for SARS-CoV-2 genome sequencing projects. <b>The result is a final SARS-CoV-2 genome ready for submission to GISAID and/or GenBank.</b>

The default workflow in EDGE COVID-19 includes:
                <ol>
                        <li> data quality control (QC) and filtering,</li>
                        <li> alignment of reads to the original (first available) reference genome (<a href="https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/" target="_new">NC_045512.2</a>),</li>
                        <li> creation of a consensus genome sequence based on the read alignments, and </li>
                        <li> a preliminary Single Nucleotide Polymorphism and Variant analyses, with some detail such as location and resulting coding differences if any. </li>
                </ol>

The [EDGE COVID-19](https://edge-covid19.edgebioinformatics.org/) platform can accommodate Illumina or ONT data, including ONT data from the <a href="https://artic.network/ncov-2019" target="_new">SARS-CoV-2 ARTIC network sequencing protocols</a>. Users can input/upload Illumina or Nanopore sequencing FASTQ files (and/or download from NCBI SRA). For Illumina data, default analyses include only read QC, read mapping to the reference, and SNP/variant analysis. For ONT data, the data must be demultiplexed prior to uploading; the samples will be processed individually.  The SNP/variant calling is not on by default for ONT. However, other functions (e.g. de novo assembly for whole genome data) are also available for both sequencing platforms.  While command line execution is possible (see <a href="https://github.com/LANL-Bioinformatics/EDGE/tree/SARS-CoV2" target="_new">here</a> and <a href="https://gitlab.com/chienchi/reference-based_assembly" target="_new">here</a>), the GUI provides an easy data submission and results viewing platform, with the graphical and tabular views of variant/SNP data and a genome browser to view read coverage and location of SNPs or variants, as well as the reference annotations.

This light-weight version is a <a href="https://hub.docker.com/r/bioedge/edge-covid19" target="_new">Docker container</a>, able to run on any local hardware infrastructure or in the cloud. We have tested this Docker container on laptops and cloud, using several Illumina (e.g. <a href="https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRR11177792" target="_new">SRR11177792</a>) and ONT (e.g. <a href="https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRR11300652" target="_new">SRR11300652</a>) datasets.
 
###### Note: For EDGE Bioinformatics users who would also like to use the phylogeny or read- and assembly-based taxonomy classification tools to identify all organisms that may be present within complex samples, we recommend using the original <a href="https://edgebioinformatics.org" target="_new">EDGE Bioinformatics platform</a> which harbors several tools and associated (large) databases that enable such a search. <em>In initial tests of taxonomy classification of SARS-CoV-2 samples (with no SARS-CoV-2 genomes in any of the databases), we recover SARS coronavirus and Bat Coronavirus as the nearest neighbor.</em> 

## Documentation
	
* [EDGE COVID-19 Documentation](https://github.com/LANL-Bioinformatics/EDGE/blob/SARS-CoV2/edge_ui/docs/EDGE_COVID-19_guide.pdf)

* [EDGE Docker Hub](https://hub.docker.com/r/bioedge/edge-covid19)

## Requirements

* Docker Engine version 19.03.2 or greater
* The image size is around 11.8GB. This image is built on top of the official Ubuntu 18.04.3 LTS Base Image.
* Recommended minimum computational resource : 8G memory, 4 CPUs, 20GB for the image.   

## Contact Info
Support: The EDGE user's group is a Google group for users and developers: [edge-users@googlegroups.com](mailto:edge-users@googlegroups.com) or email us at edge-covid19@lanl.gov

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
