# RAMIFI

<ins>R</ins>ecombinant <ins>A</ins>nd <ins>M</ins>ix-<ins>I</ins>nfection <ins>Fi</ins>nder for SARS-CoV-2 sample. It takes input from aligned bam file  (aligned to [NC_045512](https://github.com/chienchi/ramifi/blob/main/ramifi/data/NC_045512.fasta)) and output recombinant and parents reads in .bam and .tsv file with associated stats file. 

## Dependencies

### Programming/Scripting languages
- [Python >=v3.8](https://www.python.org/)
    - The pipeline has been tested in v3.8.10
    
### Python packages
- [pandas >=1.2.4](https://pandas.pydata.org/) 
- [pysam >= 0.16.0.1](https://github.com/pysam-developers/pysam)
- [importlib-resources>=5.7.1](https://pypi.org/project/importlib-resources/)

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
usage: ramifi.py [-h] [--refacc [STR]] [--minMixAF [FLOAT]] [--maxMixAF [FLOAT]] [--minMixed_n [INT]] [--minReadCount [INT]]
                 [--lineageMutation [FILE]] [--variantMutation [FILE]] [--mutations_af_plot] [--verbose] [--version] --bam [FILE]
                 [--vcf [File]] [--tsv [FILE]] [--outbam [File]] [-eo [PATH]] [--igv [PATH]] [--igv_variants]

Script to do recombinant read analysis

optional arguments:
  -h, --help            show this help message and exit
  --refacc [STR]        reference accession used in bam [default: NC_045512.2]
  --minMixAF [FLOAT]    minimum alleic frequency for checking mixed mutations on vcf [default:0.2]
  --maxMixAF [FLOAT]    maximum alleic frequency for checking mixed mutations on vcf [default:0.8]
  --minMixed_n [INT]    threshold of mixed mutations count for vcf.
  --minReadCount [INT]  threshold of read with variant count when no vcf provided.
  --lineageMutation [FILE]
                        lineage mutation json file [default: variant_mutation.json]
  --variantMutation [FILE]
                        variant mutation json file [default: lineage_mutation.json]
  --mutations_af_plot   generate mutations_af_plot (when --vcf provided)
  --verbose             Show more infomration in log
  --version             show program's version number and exit

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
  --igv_variants        add variants igv track
```

## Test

```
cd tests
./runTest.sh
```

## Outputs 

-- recombinant_reads.stats:  counts

| total  | mapped | unmapped | mutation_reads | parents     | recomb1_reads | recomb2_reads | recombx_reads | parent1_reads | parent2_reads | recomb1_perc| recomb2_perc | recombx_perc |
|--------|--------|----------|----------------|-------------|---------------|---------------|---------------|---------------|---------------|-------------|--------------|--------------|
| 64355  | 64355  |   0      |  5203          |Omicron,Delta|   162         | 175           |     18        |  489          |     730       | 10.29       | 11.11        | 1.14         |


-- recombinant_reads.tsv
|    read_name                | start | end | mutaions_json                                                                                                                                                                                                                                 |  note            |
|-----------------------------|-------|-----|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------|
|HMVN7DRXY:2:2153:21802:16078 |  21566|21859| {21618: ['Delta'], 21846: ['Iota', 'Mu', 'Omicron']}                                                                                                                                                                                          |  recombinant 2   |
|HMVN7DRXY:2:2166:28574:36229 |  21732|21883| {21762: ['Eta', 'Omicron'], 21846: ['Iota', 'Mu', 'Omicron']}                                                                                                                                                                                 |  parent Omicron  |
|HMVN7DRXY:2:2215:29749:15217 |  22867|22994| {22917: ['Delta', 'Epsilon', 'Kappa'], 22992: ['rev of Omicron']}                                                                                                                                                                             |  parent Delta    |
|HMVN7DRXY:2:2105:30572:25160 |  22865|23023| {22917: ['rev of Delta Epsilon Kappa'], 22992: ['rev of Omicron'], 22995: ['Delta', 'Omicron'], 23013: ['rev of Omicron']}                                                                                                                    |  recombinant 1   | 
|HMVN7DRXY:2:2127:18304:18850 |  24058|24518| {24130: ['Omicron'], 24469: ['rev of Omicron'], 24503: ['Omicron']}                                                                                                                                                                           |  recombinant x   |
|etc ...                      |       |     |

-- recombinant_reads_by_cross_region.tsv

| Cross_region  | Reads                                                                                                       |
|---------------|-------------------------------------------------------------------------------------------------------------|
|11201-11283    |{'recomb1': ['HMVN7DRXY:2:2150:13015:23750', 'HMVN7DRXY:2:2124:23746:28776', 'HMVN7DRXY:2:2232:6216:33395']} |
|11283-11314    |{'recomb2': ['HMVN7DRXY:2:2122:27624:23062']}                                                                |
|11537-11544    |{'recomb2': ['HMVN7DRXY:2:2126:12825:30154', 'HMVN7DRXY:2:2126:15302:29121']}                                |
|etc ...        |

-- recombinant_reads.parent1.bam

-- recombinant_reads.parent1.bam.bai

-- recombinant_reads.parent2.bam

-- recombinant_reads.parent2.bam.bai

-- recombinant_reads.recomb1.bam

-- recombinant_reads.recomb1.bam.bai

-- recombinant_reads.recomb2.bam

-- recombinant_reads.recomb2.bam.bai

-- recombinant_reads.recombx.bam

-- recombinant_reads.recombx.bam.bai

-- [recombinant_reads.mutations_af_plot.html](https://chienchi.github.io/ramifi/recombinant_reads.mutations_af_plot.html)

## Data visualization

The `recombinant_reads.bam`, `ramifi/data/variants_mutation.gff` and `ramifi/data/NC_045512.fasta` can be loaded into [IGV](https://software.broadinstitute.org/software/igv/).

Example:
IGV Link: [https://chienchi.github.io/ramifi/igv-webapp](https://chienchi.github.io/ramifi/igv-webapp)

![Screen Shot 2022-06-13 at 9 51 08 PM](https://user-images.githubusercontent.com/737589/173489713-18150a0d-176b-4526-a751-5a03d2047096.png)

## Remove package:

```
pip uninstall ramifi
```

## Citing RAMIFI

This work is currently unpublished. If you are making use of this package, we would appreciate if you gave credit to our repository.

## License

See the [LICENSE](https://github.com/chienchi/ramifi/blob/main/LICENSE) file included in the RAMIFI distribution.
