Installation
############

1. Please assure that your system has the :doc:`essential software building packages <system_requirement>`. installed properly before proceeding following installation.

2. Downloading main codes (72M), databases (7.9G , 13G, 40G, 101M for databases and GOTTCHA_db and bwa_index and SNPdb) and third party tools (~2G) from our SFTP server. The password per request. ::
    
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/edge_main_v1.0.tgz  ./  
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/edge_v1.0_thirdParty_softwares.tgz  ./  
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/edge_pipeline_v1.0.databases.tgz  ./  
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/GOTTCHA_db.tgz   ./  
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/bwa_index.tgz   ./  
    sftp -o "Port 33001" edge@img-gp.lanl.gov:/data/edge/SNPdb.tgz   ./   
    
.. warning:: Be patient, the database files are huge.

3. Unpack main archive::

    tar -xvzf edge_main_v1.0.tgz

.. note:: The main directory, edge_v1.0, will be created. 

4. Move the database and third party archives into main directory (edge_v1.0)::

    mv edge_v1.0_thirdParty_softwares.tgz edge_v1.0/
    mv edge_pipeline_v1.0.databases.tgz edge_v1.0/
    mv GOTTCHA_db.tgz edge_v1.0/
    mv bwa_index.tgz edge_v1.0/
    mv SNPdb.tgz edge_v1.0/
        
5. Change directory to main direcotry and Unpack databases and third party tools archive::
    
    cd edge_v1.0
    tar -xvzf edge_v1.0_thirdParty_softwares.tgz
    tar -xvzf edge_pipeline.v1.0.databases.tgz 
    tar -xvzf GOTTCHA_db.tgz
    tar -xzvf bwa_index.tgz
    tar -xvzf SNPdb.tgz
        
.. note:: To this point, you should see a database directory and a thirdParty directory in the main directory

6. Installing pipeline::

    ./INSTALL.sh

  It will install the following depended :doc:`tools <third_party>`.  
    
  * Assembly
  
    * idba
    * RATT

  * Annotation
  
    * prokka
    * tRNAscan
    * barrnap
    * BLAST+
    * blastall
    * phageFinder
    * glimmer
    * aragorn
    * prodigal
    * tbl2asn

  * Alignment
    
    * hmmer
    * infernal
    * bowtie2
    * bwa
    * mummer

  * Taxonomy
  
    * kraken
    * metaphlan
    * Metaphyler
    * kronatools
    * gottcha
    * metascope_plus

  * Phylogeny
  
    * FastTree
    * RAxML

  * Utility
  
    * bedtools
    * R
    * GNU_parallel
    * tabix
    * JBrowse
    * primer3
    * samtools

  * Perl_Modules
  
    * perl_parallel_forkmanager
    * perl_excel_writer
    * perl_archive_zip
    * perl_string_approx 
    * perl_pdf_api2
    * perl_tk
    * perl_html_template
    * perl_html_parser
    * perl_JSON
    * perl_bio_phylo
    * perl_xml_twig

7. Restart the Terminal Session to allow $EDGE_HOME to be exported.  

.. note:: After running INSTALL.sh successfully, the binaries and related scripts will be stored in the ./bin and ./scripts directory. It also writes EDGE_HOME environment vairable into .bashrc or .bash_profile. 
    
     
