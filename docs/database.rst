Databases
#########

1. MvirDB: A Microbial database of protein toxins, virulence factors and antibiotic resistance genes for bio-defence applications
paper: http://www.ncbi.nlm.nih.gov/pubmed/?term=17090593
website: http://mvirdb.llnl.gov/

2. blast_nucl_db and bwa_index from NCBI
bacteria: ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz
          ftp://ftp.ncbi.nih.gov/genomes/IDS/Bacteria.ids
virus: http://www.ncbi.nlm.nih.gov/genomes/GenomesHome.cgi?
       ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz

NCBI 2014 Oct 7: ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz
     2786 genomes
NCBI 2014 Oct 7: ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz 
     5591 genomes

see id_mapping.txt for all gi/accession to genome name lookup table.

3. Krona_taxonomy
Download these files from ftp://ftp.ncbi.nih.gov/pub/taxonomy:
        gi_taxid_nucl.dmp.gz
        gi_taxid_prot.dmp.gz
        taxdump.tar.gz
Transfer the files to the taxonomy folder in the standalone KronaTools installation.
Run "./thirdParty/KronaTools-2.4/updateTaxonomy.sh --local".

website: http://sourceforge.net/p/krona/home/krona/
paper: http://www.ncbi.nlm.nih.gov/pubmed/?term=21961884

4. metaphlan database
The BowTie2 database file of the MetaPhlAn database 
MetaPhlAn relies on unique clade-specific marker genes identified from 3,000 reference genomes.

website: http://huttenhower.sph.harvard.edu/metaphlan
paper: http://www.ncbi.nlm.nih.gov/pubmed/?term=22688413

5. Human Genome
The human hs_ref_GRCh38 sequences from NCBI ftp site.
ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.*.fa.gz

6. MiniKraken DB
Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies. MiniKraken is a pre-built 4 GB database constructed from complete bacterial, archaeal, and viral genomes in RefSeq (as of Mar. 30, 2014).
website: http://ccb.jhu.edu/software/kraken/

7. GOTTCHA DB
A novel, annotation-independent and signature-based metagenomic taxonomic profiling tool. (manuscript in submission)

8. SNPdb
SNP database based on whole genome comparision. Current available db are Ecoli, Yersinia, Francisella, Brucella, Bacillus.

9. Vector based bwa_index

10. Other optional database not in the EDGE package.
* NCBI nr/nt blastDB  
   ftp://ftp.ncbi.nih.gov/blast/db/

Build bwa index
===============
Here take human genome as example.

1. Download the human hs_ref_GRCh38 sequences from NCBI ftp site.
ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.*.fa.gz

Or use a provided perl script in scripts/
perl download_human_refseq_genome.pl output_dir

2. Gunzip the downloaded fasta file and concatenate them into one human genome multifasta file.

>  gunzip hs_ref_GRCh38.*.fa.gz
>  cat hs_ref_GRCh38.*.fa > human_ref_GRCh38.all.fasta

3. Use the installed bwa to build the index
>  bin/bwa index human_ref_GRCh38.all.fasta

Now, you can configure the config file with "host=/path/human_ref_GRCh38.all.fasta" for host removal step.
