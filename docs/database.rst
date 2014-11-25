Databases
#########

EDGE provided database
======================
MvirDB
------ 

A Microbial database of protein toxins, virulence factors and antibiotic resistance genes for bio-defence applications

* paper: `http://www.ncbi.nlm.nih.gov/pubmed/?term=17090593 <http://www.ncbi.nlm.nih.gov/pubmed/?term=17090593>`_
* website: `http://mvirdb.llnl.gov/ <http://mvirdb.llnl.gov/>`_

NCBI Refseq
-----------

EDGE prebuilt the blast db and bwa_index of NCBI RefSeq genomes.

* Bacteria: `ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz <ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.fna.tar.gz>`_

  * Version: NCBI 2014 Oct 7
  * 2786 genomes
  
* Virus:  `ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz <ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.fna.tar.gz>`_

  * Version: NCBI 2014 Oct 7
  * 5591 genomes

see $EDGE_HOME/database/bwa_index/id_mapping.txt for all gi/accession to genome name lookup table.

Krona taxonomy
--------------

* paper: `http://www.ncbi.nlm.nih.gov/pubmed/?term=21961884 <http://www.ncbi.nlm.nih.gov/pubmed/?term=21961884>`_
* website: `http://sourceforge.net/p/krona/home/krona/ <http://sourceforge.net/p/krona/home/krona/>`_

Update Krona taxonomy db
^^^^^^^^^^^^^^^^^^^^^^^^

Download these files from `ftp://ftp.ncbi.nih.gov/pub/taxonomy <ftp://ftp.ncbi.nih.gov/pub/taxonomy>`_::

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    
Transfer the files to the taxonomy folder in the standalone KronaTools installation and run::

    $EDGE_HOME/thirdParty/KronaTools-2.4/updateTaxonomy.sh --local.



Metaphlan database
------------------

MetaPhlAn relies on unique clade-specific marker genes identified from 3,000 reference genomes.

* paper: `http://www.ncbi.nlm.nih.gov/pubmed/?term=22688413 <http://www.ncbi.nlm.nih.gov/pubmed/?term=22688413>`_
* website: `http://huttenhower.sph.harvard.edu/metaphlan <http://huttenhower.sph.harvard.edu/metaphlan>`_

Human Genome
------------
The bwa index is prebuilt in the EDGE.
The human hs_ref_GRCh38 sequences from NCBI ftp site.

* website `ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/ <ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/>`_

MiniKraken DB
-------------

Kraken is a system for assigning taxonomic labels to short DNA sequences, usually obtained through metagenomic studies. MiniKraken is a pre-built 4 GB database constructed from complete bacterial, archaeal, and viral genomes in RefSeq (as of Mar. 30, 2014).

* paper: `http://www.ncbi.nlm.nih.gov/pubmed/?term=24580807 <http://www.ncbi.nlm.nih.gov/pubmed/?term=24580807>`_
* website: `http://ccb.jhu.edu/software/kraken/ <http://ccb.jhu.edu/software/kraken/>`_

GOTTCHA DB
----------

A novel, annotation-independent and signature-based metagenomic taxonomic profiling tool. (manuscript in submission)

* website: `https://github.com/LANL-Bioinformatics/GOTTCHA <https://github.com/LANL-Bioinformatics/GOTTCHA>`_

SNPdb
-----

SNP database based on whole genome comparision. Current available db are Ecoli, Yersinia, Francisella, Brucella, Bacillus.

Invertebrate Vectors of Human Pathogens
---------------------------------------

The bwa index is prebuilt in the EDGE.

* paper: `http://www.ncbi.nlm.nih.gov/pubmed/?term=22135296 <http://www.ncbi.nlm.nih.gov/pubmed/?term=22135296>`_
* website: `https://www.vectorbase.org <https://www.vectorbase.org>`_

Version: 2014 July 24

Other optional database
-----------------------

Not in the EDGE but you can download.

* NCBI nr/nt blastDB: `ftp://ftp.ncbi.nih.gov/blast/db/ <ftp://ftp.ncbi.nih.gov/blast/db/>`_

.. _build-host-index:

Build bwa index
===============
Here take human genome as example.

1. Download the human hs_ref_GRCh38 sequences from NCBI ftp site.

  Go to `ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/ <ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/>`_
  Or use a provided perl script in $EDGE_HOME/scripts/ ::

    perl $EDGE_HOME/scripts/download_human_refseq_genome.pl output_dir

2. Gunzip the downloaded fasta file and concatenate them into one human genome multifasta file::

    gunzip hs_ref_GRCh38.*.fa.gz
    cat hs_ref_GRCh38.*.fa > human_ref_GRCh38.all.fasta

3. Use the installed bwa to build the index::

    $EDGE_HOME/bin/bwa index human_ref_GRCh38.all.fasta

  Now, you can configure the config file with "host=/path/human_ref_GRCh38.all.fasta" for host removal step.
