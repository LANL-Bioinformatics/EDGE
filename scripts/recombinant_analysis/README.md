# ramifi

<ins>R</ins>ecombinant <ins>A</ins>nd <ins>M</ins>ix-<ins>I</ins>nfection <ins>FI</ins>nder

## Dependencies

### Programming/Scripting languages
- [Python >=v3.8](https://www.python.org/)
    - The pipeline has been tested in v3.8.10
    
### Python packages
- [pandas >=1.2.4](https://pandas.pydata.org/) 
- [pysam >= 0.16.0.1](https://github.com/pysam-developers/pysam)

#### Optional packages
- [plotly >=4.7.1](https://plotly.com/python/)
- [kaleido >= 0.2.1](https://github.com/plotly/Kaleido)
- [biopython >= 1.78](https://biopython.org/)


## Installation

### Install from source
Clone the `ramifi` repository.

```
git clone https://github.com/chienchi/ramifi
```

Then change directory to `ramifi` and install.

```
cd ramifi
pip install .
```

If the installation was succesful, you should be able to type `ramifi -h` and get a help message on how to use the tool.

```
ramifi -h
```


## Usage
```
usage: ramifi [-h] [--refacc [STR]] [--minMixAF [FLOAT]] [--maxMixAF [FLOAT]] [--minMixed_n [INT]] [--minReadCount [INT]]
              [--lineageMutation [FILE]] [--variantMutation [FILE]] [--verbose] --bam [FILE] [--vcf [File]] [--tsv [FILE]]
              [--outbam [File]] [-eo [PATH]] [--igv [PATH]]

Script to do recombinant read analysis

optional arguments:
  -h, --help            show this help message and exit
  --refacc [STR]        reference accession used in bam [default: NC_045512.2]
  --minMixAF [FLOAT]    minimum alleic frequency for checking mixed mutations on vcf [default:0.2]
  --maxMixAF [FLOAT]    maximum alleic frequency for checking mixed mutations on vcf [default:0.8]
  --minMixed_n [INT]    threadhold of mixed mutations count for vcf.
  --minReadCount [INT]  threadhold of read with variant count when no vcf provided.
  --lineageMutation [FILE]
                        lineage mutation json file [default: variant_mutation.json]
  --variantMutation [FILE]
                        variant mutation json file [default: lineage_mutation.json]
  --verbose             Show more infomration in log

Input:
  --bam [FILE]          <Required> bam file
  --vcf [File]          <Optional> vcf file which will infer the two parents of recombinant_variants

Output:
  --tsv [FILE]          output file name [default: recombinant_reads.tsv]
  --outbam [File]       output recombinant reads in bam file [default: recombinant_reads.bam]

EDGE COVID-19 Options:
  options specific used for EDGE COVID-19

  -eo [PATH], --ec19_projdir [PATH]
                        ec-19 project directory
  --igv [PATH]          igv.html relative path
```

## Test

```
cd tests
./runTest.sh
```

## Outputs 

-- recombinant_reads.stats:  counts

| total  | mapped | unmapped | mutation_reads | parents     | recomb_reads | parent1_reads | parent2_reads | recomb_perc |
|--------|--------|----------|----------------|-------------|--------------|---------------|---------------|-------------|
| 64355  | 64355  |   0      |  12075         |Delta,Omicron|   249        |  617          |     376       | 20.048309178|


-- recombinant_reads.tsv
|    read_name                | start | end | mutaions_json                                                                                                                                                                                                                                 |  note            |
|-----------------------------|-------|-----|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
|HMVN7DRXY:2:2153:21802:16078 |  21566|21859| {21575: ['ref of  Iota'], 21600: ['ref of  Epsilon'], 21614: ['ref of  Gamma'], 21618: ['Delta'], 21621: ['ref of  Gamma'], 21786: ['ref of  Lambda'], 21789: ['ref of  Lambda'], 21801: ['ref of  Beta'], 21846: ['Iota', 'Mu', 'Omicron']}  |  recombinant     |
|HMVN7DRXY:2:2166:28574:36229 |  21732|21883| {21762: ['Eta', 'Omicron'], 21786: ['ref of  Lambda'], 21789: ['ref of  Lambda'], 21801: ['ref of  Beta'], 21846: ['Iota', 'Mu', 'Omicron']}                                                                                                  |  parent Omicron  |
|HMVN7DRXY:2:2108:19795:22373 |  22605|22717| {22679: ['ref of  Omicron'], 22686: ['ref of  Omicron']}                                                                                                                                                                                      |  parent Delta    |
|
|etc ...

-- recombinant_reads.bam

-- recombinant_reads.bam.bai

## Data visualization

The `recombinant_reads.bam`, `ramifi/data/variants_mutation.gff` and `ramifi/data/NC_045512.fasta` can be loaded into IGV.

## Remove package:

```
pip uninstall ramifi
```
