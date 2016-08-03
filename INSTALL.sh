#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )
exec >  >(tee install.log)
exec 2>&1

cd $rootdir
cd thirdParty

mkdir -p $rootdir/bin

export PATH=$PATH:$rootdir/bin/

assembly_tools=( idba spades )
annotation_tools=( prokka RATT tRNAscan barrnap BLAST+ blastall phageFinder glimmer aragorn prodigal tbl2asn )
utility_tools=( bedtools R GNU_parallel tabix JBrowse primer3 samtools sratoolkit )
alignments_tools=( hmmer infernal bowtie2 bwa mummer )
taxonomy_tools=( kraken metaphlan kronatools gottcha )
phylogeny_tools=( FastTree RAxML )
perl_modules=( perl_parallel_forkmanager perl_excel_writer perl_archive_zip perl_string_approx perl_pdf_api2 perl_html_template perl_html_parser perl_JSON perl_bio_phylo perl_xml_twig perl_cgi_session )
all_tools=("${assembly_tools[@]}" "${annotation_tools[@]}" "${utility_tools[@]}" "${alignments_tools[@]}" "${taxonomy_tools[@]}" "${phylogeny_tools[@]}" "${perl_modules[@]}")

### Install functions ###
install_idba()
{
echo "------------------------------------------------------------------------------
                           Compiling IDBA 1.1.1
------------------------------------------------------------------------------
"
tar xvzf idba-1.1.1.tar.gz
cd idba-1.1.1
sed -i.bak 's/kMaxShortSequence = 128/kMaxShortSequence = 351/' src/sequence/short_sequence.h
sed -i.bak 's/kNumUint64 = 4/kNumUint64 = 6/' src/basic/kmer.h
#src/sequence/short_sequence.h:    static const uint32_t kMaxShortSequence = 128
./configure --prefix=$rootdir
make 
make install
cp bin/idba_ud $rootdir/bin/.
cp bin/fq2fa $rootdir/bin/.
cd $rootdir/thirdParty
if [[ "$OSTYPE" == "darwin"* ]]
then
{
      cp -f idba_ud_mac $rootdir/bin/idba_ud
}
fi
echo "
------------------------------------------------------------------------------
                           IDBA compiled
------------------------------------------------------------------------------
"
}


install_spades(){
echo "------------------------------------------------------------------------------
                           Installing SPAdes 3.5.0
------------------------------------------------------------------------------
"
tar xvzf SPAdes-3.5.0-Linux.tar.gz 
ln -sf $rootdir/thirdParty/SPAdes-3.5.0-Linux/bin/spades.py $rootdir/bin/spades.py
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           SPAdes installed
------------------------------------------------------------------------------
"
}

install_tRNAscan()
{
echo "------------------------------------------------------------------------------
                           Installing tRNAscan-SE 1.3.1
------------------------------------------------------------------------------
"
tar xvzf tRNAscan-SE-1.3.1.tar.gz
cd tRNAscan-SE-1.3.1
sed -i.bak 's,home,'"$rootdir"',' Makefile
make
make install
make clean
cd $rootdir/thirdParty
chmod -R +x $rootdir/bin/tRNAscanSE
chmod -R +r $rootdir/bin/tRNAscanSE
echo "
------------------------------------------------------------------------------
                           tRNAscan-SE 1.3.1 installed
------------------------------------------------------------------------------
"
}

install_prokka()
{
echo "------------------------------------------------------------------------------
                           Installing prokka-1.11
------------------------------------------------------------------------------
"
tar xvzf prokka-1.11.tar.gz
cd prokka-1.11
cd $rootdir/thirdParty
ln -sf $rootdir/thirdParty/prokka-1.11/bin/prokka $rootdir/bin/prokka
$rootdir/thirdParty/prokka-1.11/bin/prokka --setupdb
echo "
------------------------------------------------------------------------------
                           prokka-1.11 installed
------------------------------------------------------------------------------
"
}

install_barrnap()
{
echo "------------------------------------------------------------------------------
                           Installing barrnap-0.4.2
------------------------------------------------------------------------------
"
tar xvzf barrnap-0.4.2.tar.gz
cd barrnap-0.4.2
cd $rootdir/thirdParty
ln -sf $rootdir/thirdParty/barrnap-0.4.2/bin/barrnap $rootdir/bin/barrnap
echo "
------------------------------------------------------------------------------
                           barrnap-0.4.2 installed
------------------------------------------------------------------------------
"
}


install_bedtools()
{
echo "------------------------------------------------------------------------------
                           Installing bedtools-2.19.1
------------------------------------------------------------------------------
"
tar xvzf bedtools-2.19.1.tar.gz
cd bedtools2-2.19.1
make 
cp -fR bin/* $rootdir/bin/. 
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           bedtools-2.19.1 installed
------------------------------------------------------------------------------
"
}

install_sratoolkit()
{
echo "------------------------------------------------------------------------------
                           Installing sratoolkit.2.4.4-linux64
------------------------------------------------------------------------------
"
tar xvzf sratoolkit.2.4.4-linux64.tgz
cd sratoolkit.2.4.4-linux64
ln -sf $rootdir/thirdParty/sratoolkit.2.4.4-linux64/bin/fastq-dump $rootdir/bin/fastq-dump
ln -sf $rootdir/thirdParty/sratoolkit.2.4.4-linux64/bin/vdb-dump $rootdir/bin/vdb-dump
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           sratoolkit.2.4.3-linux64 installed
------------------------------------------------------------------------------
"
}

install_R()
{
echo "------------------------------------------------------------------------------
                           Compiling R 2.15.3
------------------------------------------------------------------------------
"
tar xvzf R-2.15.3.tar.gz
cd R-2.15.3
./configure --prefix=$rootdir --with-readline=no 
make
make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           R compiled
------------------------------------------------------------------------------
"
}

install_GNU_parallel()
{
echo "------------------------------------------------------------------------------
                           Compiling GNU parallel
------------------------------------------------------------------------------
"
tar xvzf parallel-20140622.tar.gz
cd parallel-20140622
./configure --prefix=$rootdir 
make
make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           GNU parallel compiled
------------------------------------------------------------------------------
"
}

install_BLAST+()
{
echo "------------------------------------------------------------------------------
                           Install ncbi-blast-2.2.28+-x64
------------------------------------------------------------------------------
"
BLAST_ZIP=ncbi-blast-2.2.29+-x64-linux.tar.gz
if [[ "$OSTYPE" == "darwin"* ]]
then
{
    BLAST_ZIP=ncbi-blast-2.2.29+-universal-macosx.tar.gz
}
fi

tar xvzf $BLAST_ZIP
cd ncbi-blast-2.2.29+
cp -fR bin/* $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           ncbi-blast-2.2.28+-x64 installed
------------------------------------------------------------------------------
"
}

install_blastall()
{
echo "------------------------------------------------------------------------------
                           Install blast-2.2.26-x64-linux
------------------------------------------------------------------------------
"
BLAST_ZIP=blast-2.2.26-x64-linux.tar.gz
if [[ "$OSTYPE" == "darwin"* ]]
then
{
    BLAST_ZIP=blast-2.2.26-universal-macosx.tar.gz
}
fi

tar xvzf $BLAST_ZIP
cd blast-2.2.26
cp -fR bin/* $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           blast-2.2.26-x64 installed
------------------------------------------------------------------------------
"
}

install_kraken()
{
echo "------------------------------------------------------------------------------
                           Install kraken-0.10.4-beta
------------------------------------------------------------------------------
"
tar xvzf kraken-0.10.4-beta.tgz
cd kraken-0.10.4-beta
./install_kraken.sh $rootdir/bin/

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           kraken-0.10.4-beta installed
------------------------------------------------------------------------------
"
}

install_JBrowse()
{
echo "------------------------------------------------------------------------------
                           Installing JBrowse-1.11.6
------------------------------------------------------------------------------
"
tar xvzf JBrowse-1.11.6.tar.gz
cd JBrowse-1.11.6
./setup.sh
mkdir -p -m 775 data
cd $rootdir/thirdParty
if [ -e $rootdir/edge_ui/JBrowse ]
then
  rm $rootdir/edge_ui/JBrowse
fi
ln -sf $rootdir/thirdParty/JBrowse-1.11.6 $rootdir/edge_ui/JBrowse
echo "
------------------------------------------------------------------------------
                           JBrowse-1.11.6 installed
------------------------------------------------------------------------------
"
}

install_tabix()
{
echo "------------------------------------------------------------------------------
                           Compiling tabix bgzip
------------------------------------------------------------------------------
"
tar xvzf tabix.tgz
cd tabix
make
cp tabix $rootdir/bin/.
cp bgzip $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           tabix bgzip  compiled
------------------------------------------------------------------------------
"
}

install_hmmer()
{
echo "------------------------------------------------------------------------------
                           Compiling hmmer-3.1b1
------------------------------------------------------------------------------
"
tar xvzf hmmer-3.1b1.tar.gz
cd hmmer-3.1b1/
./configure --prefix=$rootdir && make && make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           hmmer-3.1b1 compiled
------------------------------------------------------------------------------
"
}

install_infernal()
{
echo "------------------------------------------------------------------------------
                           Installing infernal-1.1rc4
------------------------------------------------------------------------------
"
tar xzvf infernal-1.1rc4.tar.gz
cd infernal-1.1rc4/
./configure --prefix=$rootdir && make && make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           infernal-1.1rc4 installed
------------------------------------------------------------------------------
"
}

install_phageFinder()
{
echo "------------------------------------------------------------------------------
                           Installing phage_finder_v2.1
------------------------------------------------------------------------------
"
tar xvzf phage_finder_v2.1.tar.gz
cd phage_finder_v2.1
cd $rootdir/thirdParty
chmod -R +x phage_finder_v2.1
chmod -R +r phage_finder_v2.1
echo "
------------------------------------------------------------------------------
                           phage_finder_v2.1 installed
------------------------------------------------------------------------------
"
}

install_bowtie2()
{
echo "------------------------------------------------------------------------------
                           Compiling bowtie2 2.1.0
------------------------------------------------------------------------------
"
tar xvzf bowtie2-2.1.0.tar.gz
cd bowtie2-2.1.0
make
cp bowtie2 $rootdir/bin/.
cp bowtie2-build $rootdir/bin/.
cp bowtie2-align $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           bowtie2 compiled
------------------------------------------------------------------------------
"
}

install_gottcha()
{
echo "------------------------------------------------------------------------------
                           Compiling gottcha-1.0b
------------------------------------------------------------------------------
"
tar xvzf gottcha.tar.gz
cd gottcha
./INSTALL.sh
ln -sf $PWD/bin/gottcha.pl $rootdir/bin/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           gottcha-1.0b compiled
------------------------------------------------------------------------------
"
}

install_metaphlan()
{
echo "------------------------------------------------------------------------------
                           Compiling metaphlan-1.7.7
------------------------------------------------------------------------------
"
tar xvzf metaphlan-1.7.7.tar.gz
cd metaphlan-1.7.7
cp -fR metaphlan.py $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           metaphlan-1.7.7 compiled
------------------------------------------------------------------------------
"
}

install_RATT()
{
echo "------------------------------------------------------------------------------
                           Installing RATT 
------------------------------------------------------------------------------
"
tar xvzf RATT.tar.gz
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           RATT installed 
------------------------------------------------------------------------------
"
}

install_glimmer()
{
echo "------------------------------------------------------------------------------
                           Compiling glimmer 3.02
------------------------------------------------------------------------------
"
tar xvzf glimmer302b.tar.gz
cd glimmer3.02/SimpleMake
make
cp ../bin/* $rootdir/bin/.
cp ../scripts/* $rootdir/scripts/.
for i in $rootdir/scripts/*.csh
do 
 sed -i.bak 's!/fs/szgenefinding/Glimmer3!'$rootdir'!' $i
 sed -i.bak 's!/fs/szgenefinding/Glimmer3!'$rootdir'!' $i
done
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           glimmer3.02 compiled
------------------------------------------------------------------------------
"
}

install_aragorn()
{
echo "------------------------------------------------------------------------------
                           Compiling aragorn1.2.36
------------------------------------------------------------------------------
"
tar xvzf aragorn1.2.36.tgz
cd aragorn1.2.36
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.36.c
cp -fR aragorn $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           aragorn1.2.36 compiled
------------------------------------------------------------------------------
"
}

install_prodigal()
{
echo "------------------------------------------------------------------------------
                           Compiling prodigal.v2_60
------------------------------------------------------------------------------
"
tar xvzf prodigal.v2_60.tar.gz
cd prodigal.v2_60/
make
cp -fR prodigal $rootdir/bin
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           prodigal.v2_60 compiled
------------------------------------------------------------------------------
"
}

install_tbl2asn()
{
echo "------------------------------------------------------------------------------
                           Installing NCBI tbl2asn
------------------------------------------------------------------------------
"
if [[ "$OSTYPE" == "darwin"* ]]
then
{
    tar xvzf mac.tbl2asn.tgz
    chmod +x mac.tbl2asn
    ln -sf $rootdir/thirdParty/mac.tbl2asn $rootdir/bin/tbl2asn
}
else
{
    tar xvzf linux64.tbl2asn.tgz
    chmod +x linux64.tbl2asn
    ln -sf $rootdir/thirdParty/linux64.tbl2asn $rootdir/bin/tbl2asn.orig
}
fi

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           NCBI tbl2asn installed
------------------------------------------------------------------------------
"
}

install_bwa()
{
echo "------------------------------------------------------------------------------
                           Compiling bwa 0.7.12
------------------------------------------------------------------------------
"
tar xvzf bwa-0.7.12.tar.gz
cd bwa-0.7.12
make clean && make
cp bwa $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           bwa compiled
------------------------------------------------------------------------------
"
}

install_mummer()
{
echo "------------------------------------------------------------------------------
                           Compiling MUMmer3.23 64bit
------------------------------------------------------------------------------
"
tar xvzf MUMmer3.23.tar.gz
cd MUMmer3.23
#for 64bit MUMmer complie
make CPPFLAGS="-O3 -DSIXTYFOURBITS"
cp nucmer $rootdir/bin/.
cp show-coords $rootdir/bin/.
cp show-snps $rootdir/bin/.
cp mgaps $rootdir/bin/.
cp delta-filter $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           MUMmer3.23 compiled
------------------------------------------------------------------------------
"
}

install_primer3()
{
echo "------------------------------------------------------------------------------
                           Compiling primer3 2.3.5
------------------------------------------------------------------------------
"
tar xvzf primer3-2.3.5.tar.gz
cd primer3-2.3.5/src
make
cp primer3_core $rootdir/bin/.
cp oligotm $rootdir/bin/.
cp ntthal $rootdir/bin/.
cp ntdpal $rootdir/bin/.
cp long_seq_tm_test $rootdir/bin/.
cp -R primer3_config $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                          primer3-2.3.5 compiled
------------------------------------------------------------------------------
"
}

install_kronatools()
{
echo "------------------------------------------------------------------------------
               Installing KronaTools-2.4
------------------------------------------------------------------------------
"
tar xvzf KronaTools-2.4.tar.gz
cd KronaTools-2.4
perl install.pl --prefix $rootdir --taxonomy $rootdir/database/Krona_taxonomy
#./updateTaxonomy.sh --local
cp $rootdir/scripts/microbial_profiling/script/ImportBWA.pl scripts/
ln -sf $rootdir/thirdParty/KronaTools-2.4/scripts/ImportBWA.pl $rootdir/bin/ktImportBWA 
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        KronaTools-2.4 Installed
------------------------------------------------------------------------------
"
}

install_samtools()
{
echo "------------------------------------------------------------------------------
                           Compiling samtools 0.1.19
------------------------------------------------------------------------------
"
tar xvzf samtools-0.1.19.tar.gz
cd samtools-0.1.19
make
cp samtools $rootdir/bin/.
cp bcftools/bcftools $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           samtools compiled
------------------------------------------------------------------------------
"
}

install_FastTree()
{
echo "------------------------------------------------------------------------------
                           Compiling FastTree
------------------------------------------------------------------------------
"
gcc -DOPENMP -DUSE_DOUBLE -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm
cp -f FastTreeMP $rootdir/bin/.
echo "
------------------------------------------------------------------------------
                           FastTree compiled
------------------------------------------------------------------------------
"
}

install_RAxML()
{
echo "------------------------------------------------------------------------------
                           Compiling RAxML-8.0.26
------------------------------------------------------------------------------
"
tar xvzf RAxML-8.0.26.tar.gz
cd RAxML-8.0.26
make -f Makefile.PTHREADS.gcc
cp -f raxmlHPC-PTHREADS $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           RAxML-8.0.26 compiled
------------------------------------------------------------------------------
"
}


install_perl_parallel_forkmanager()
{
echo "------------------------------------------------------------------------------
               Installing Perl Module Parallel-ForkManager-1.03
------------------------------------------------------------------------------
"
tar xvzf Parallel-ForkManager-1.03.tar.gz
cd Parallel-ForkManager-1.03
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        Parallel-ForkManager-1.03 Installed
------------------------------------------------------------------------------
"
}

install_perl_excel_writer()
{
echo "------------------------------------------------------------------------------
               Installing Perl Module Excel-Writer-XLSX-0.71
------------------------------------------------------------------------------
"
tar xvzf Excel-Writer-XLSX-0.71.tar.gz
cd Excel-Writer-XLSX-0.71
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cp blib/script/extract_vba $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        Excel-Writer-XLSX-0.71 Installed
------------------------------------------------------------------------------
" 
}

install_perl_archive_zip()
{
echo "------------------------------------------------------------------------------
               Installing Perl Module Archive-Zip-1.37
------------------------------------------------------------------------------
"
tar xvzf Archive-Zip-1.37.tar.gz
cd Archive-Zip-1.37
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        Archive-Zip-1.37 Installed
------------------------------------------------------------------------------
"
}


install_perl_string_approx()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module String-Approx-3.27
------------------------------------------------------------------------------
"
tar xvzf String-Approx-3.27.tar.gz
cd String-Approx-3.27
perl Makefile.PL 
make
cp -fR blib/lib/* $rootdir/lib/
mkdir -p $rootdir/lib/auto
mkdir -p $rootdir/lib/auto/String
mkdir -p $rootdir/lib/auto/String/Approx
cp -fR blib/arch/auto/String/Approx/Approx.* $rootdir/lib/auto/String/Approx/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        String-Approx-3.27 Installed
------------------------------------------------------------------------------
"
}

install_perl_pdf_api2()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module PDF-API2-2.020
------------------------------------------------------------------------------
"
tar xvzf PDF-API2-2.020.tar.gz
cd PDF-API2-2.020
perl Makefile.PL 
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        PDF-API2-2.020 Installed
------------------------------------------------------------------------------
"
}

install_perl_JSON()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module JSON-2.90
------------------------------------------------------------------------------
"
tar xvzf JSON-2.90.tar.gz
cd JSON-2.90
perl Makefile.PL 
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        JSON-2.90 Installed
------------------------------------------------------------------------------
"
}

install_perl_html_parser()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module HTML-Parser-3.71
------------------------------------------------------------------------------
"
tar xvzf HTML-Parser-3.71.tar.gz
cd HTML-Parser-3.71
perl Makefile.PL 
make
cp -fR blib/lib/* $rootdir/lib/.
cp -fR blib/arch/auto/* $rootdir/lib/auto/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        HTML-Parser-3.71 Installed
------------------------------------------------------------------------------
"
}

install_perl_html_template()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module HTML-Template-2.6
------------------------------------------------------------------------------
"
tar xvzf HTML-Template-2.6.tar.gz
cd HTML-Template-2.6
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                         HTML-Template-2.6 Installed
------------------------------------------------------------------------------
"
}

install_perl_bio_phylo()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module Bio-Phylo-0.58
------------------------------------------------------------------------------
"
tar xvzf Bio-Phylo-0.58.tar.gz
cd Bio-Phylo-0.58
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                         Bio-Phylo-0.58 Installed
------------------------------------------------------------------------------
"
}

install_perl_xml_twig()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module XML-Twig-3.48
------------------------------------------------------------------------------
"
tar xvzf XML-Twig-3.48.tar.gz
cd XML-Twig-3.48
perl Makefile.PL -y
make
cp -fR blib/lib/* $rootdir/lib/.
cp -fR blib/script/* $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                         XML-Twig-3.48 Installed
------------------------------------------------------------------------------
"
}

install_perl_cgi_session()
{
echo "------------------------------------------------------------------------------
                 Installing Perl Module CGI-Session-4.48
------------------------------------------------------------------------------
"
tar xvzf CGI-Session-4.48.tar.gz
cd CGI-Session-4.48
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                         CGI-Session-4.48 Installed
------------------------------------------------------------------------------
"
}

checkSystemInstallation()
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkLocalInstallation()
{
    IFS=:
    for d in $rootdir/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkPerlModule()
{
   perl -e "use lib \"$rootdir/lib\"; use $1;"
   return $?
}


containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}

print_usage()
{
cat << EOF
usage: $0 options
    If no options, it will check existing installation and run tools installation for those uninstalled.
    options:
    help            show this help
    list            show available tools for updates
    tools_name      install/update individual tool
    force           force to install all list tools locally
    
    ex: To update bowtie2 only
        $0 bowtie2
    ex: To update bowtie2 and bwa
        $0 bowtie2 bwa
    ex: RE-install Phylogeny tools
        $0 Phylogeny
        
EOF

}

print_tools_list()
{

   
   echo "Available tools for updates/re-install"
   echo -e "\nAssembly"
   for i in "${assembly_tools[@]}"
   do
	   echo "* $i"
   done
   echo -e "\nAnnotation"
   for i in "${annotation_tools[@]}"
   do
	   echo "* $i"
   done
   echo -e "\nAlignment"
   for i in "${alignments_tools[@]}"
   do
	   echo "* $i"
   done
   echo -e "\nTaxonomy"
   for i in "${taxonomy_tools[@]}"
   do
	  echo "* $i"
   done
   echo -e "\nPhylogeny"
   for i in "${phylogeny_tools[@]}"
   do
	   echo "* $i"
   done
   echo -e "\nUtility"
   for i in "${utility_tools[@]}"
   do
	   echo "* $i"
   done
   echo -e "\nPerl_Modules"
   for i in "${perl_modules[@]}"
   do
	   echo "* $i"
   done
}


### Main ####
if ( checkSystemInstallation csh )
then
  #echo "csh is found"
  echo ""
else
  echo "csh is not found"
  echo "Please Install csh first, then INSTALL the package"
  exit 1
fi

if perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"' >/dev/null 2>&1 
then 
  #perl -MBio::Root::Version -e 'print "BioPerl Version ", $Bio::Root::Version::VERSION," is found\n"'
  echo ""
else 
  echo "Cannot find a perl Bioperl Module installed" 1>&2
  echo "Please install Bioperl (http://www.bioperl.org/)"
  exit 1
fi

if [ "$#" -ge 1 ]
then
  for f in $@
  do
    case $f in
      help)
        print_usage
        exit 0;;
      list)
        print_tools_list
        exit 0;;
      Assembly)
        for tool in "${assembly_tools[@]}"
        do
            install_$tool
        done
        echo -e "Assembly tools installed.\n"
        exit 0;;  
      Annotation)
        for tool in "${annotation_tools[@]}"
        do
            install_$tool
        done
        echo -e "Annotation tools installed.\n"
        exit 0;;  
      Alignment)
        for tool in "${alignments_tools[@]}"
        do
            install_$tool
        done
        echo -e "Alignment tools installed.\n"
        exit 0;; 
      Taxonomy)
        for tool in "${taxonomy_tools[@]}"
        do
            install_$tool
        done
        echo -e "Taxonomy tools installed.\n"
        exit 0;; 
      Phylogeny)
        for tool in "${phylogeny_tools[@]}"
        do
            install_$tool
        done
        echo -e "Phylogeny tools installed.\n"
        exit 0;; 
      Utility)
        for tool in "${utility_tools[@]}"
        do
            install_$tool
        done
        echo -e "Utility tools installed.\n"
        exit 0;; 
      Perl_Modules)
        for tool in "${perl_modules[@]}"
        do
            install_$tool
        done
        echo -e "Perl_Modules installed.\n"
        exit 0;;
      force)
        for tool in "${all_tools[@]}" 
        do
            install_$tool
        done
        ;;
      *)
        if ( containsElement "$f" "${assembly_tools[@]}" || containsElement "$f" "${annotation_tools[@]}" || containsElement "$f" "${alignments_tools[@]}" || containsElement "$f" "${taxonomy_tools[@]}" || containsElement "$f" "${phylogeny_tools[@]}" || containsElement "$f" "${utility_tools[@]}" || containsElement "$f" "${perl_modules[@]}" )
        then
            install_$f
        else
            echo "$f: no this tool in the list"
            print_tools_list
        fi
        exit 0;;
    esac
  done
fi

if ( checkSystemInstallation inkscape )
then
  echo "inkscape is found"
else
  echo "inkscape is not found"
 # echo "Please Install inkscape, then INSTALL the package"
 # exit 1
fi

if [[ "$OSTYPE" == "darwin"* ]]
then
{
    if ( checkSystemInstallation R )
    then
    {
        echo "R is found"
    }
    else
    {
        echo "R is not found"
        echo "Please install R from http://cran.r-project.org/bin/macosx/";
        exit 1
    }
    fi
}
else
{
    if ( checkLocalInstallation R )
    then
    {
        echo "R is found"
    }
    else
    {
        install_R
    }
    fi
}
fi

echo "if(\"gridExtra\" %in% rownames(installed.packages()) == FALSE)  {install.packages(\"gridExtra_0.9.1.tar.gz\", repos = NULL, type=\"source\")}" | Rscript -  

if ( checkSystemInstallation bedtools )
then
  echo "bedtools is found"
else
  echo "bedtools is not found"
  install_bedtools 
fi

if ( checkSystemInstallation fastq-dump )
then
  echo "sratoolkit is found"
else
  echo "sratoolkit is not found"
  install_sratoolkit 
fi

if ( checkSystemInstallation parallel )
then
  echo "GNU parallel is found"
else
  echo "GNU parallel is not found"
  install_GNU_parallel
fi

if ( checkSystemInstallation blastn )
then
  echo "BLAST+ is found"
else
  echo "BLAST+ is not found"
  install_BLAST+
fi

if ( checkSystemInstallation blastall )
then
  echo "blastall is found"
else
  echo "blastall is not found"
  install_blastall
fi

if ( checkSystemInstallation tRNAscan-SE )
then
  echo "tRNAscan-SE is found"
else
  echo "tRNAscan-SE is not found"
  install_tRNAscan
fi

if ( checkLocalInstallation ktImportBLAST )
then
  echo "KronaTools is found"
else
  echo "KronaTools is not found"
  install_kronatools
fi


if ( checkLocalInstallation hmmpress )
then
  echo "hmmer3 is found"
else
  echo "hmmer3 is not found"
  install_hmmer
fi

if ( checkLocalInstallation cmbuild )
then
  echo "infernal-1.1 is found"
else
  echo "infernal-1.1 is not found"
  install_infernal
fi

if ( checkLocalInstallation prokka )
then
  echo "prokka is found"
else
  echo "prokka is not found"
  install_prokka
fi

if [ -x $rootdir/thirdParty/RATT/start.ratt.sh   ]
then
  echo "RATT is found"
else
  echo "RATT is not found"
  install_RATT
fi

if ( checkSystemInstallation barrnap )
then
  echo "barrnap is found"
else
  echo "barrnap is not found"
  install_barrnap
fi

if ( checkSystemInstallation glimmer3 )
then
  echo "glimmer is found"
else
  echo "glimmer is not found"
  install_glimmer
fi

if ( checkSystemInstallation prodigal )
then
  echo "prodigal is found"
else
  echo "prodigal is not found"
  install_prodigal
fi

if ( checkSystemInstallation aragorn )
then
  echo "aragorn is found"
else
  echo "aragorn is not found"
  install_aragorn
fi

if ( checkSystemInstallation tbl2asn.orig )
then
  echo "tbl2asn is found"
else
  echo "tbl2asn is not found"
  install_tbl2asn
fi


if ( checkLocalInstallation kraken )
then
  echo "kraken is found"
else
  echo "kraken is not found"
  install_kraken
fi

if ( checkSystemInstallation tabix )
then
  echo "tabix is found"
else
  echo "tabix is not found"
  install_tabix
fi

if ( checkSystemInstallation bgzip )
then
  echo "bgzip is found"
else
  echo "bgzip is not found"
  install_tabix
fi

if ( checkSystemInstallation bowtie2 )
then
  echo "bowtie2 is found"
else
  echo "bowtie2 is not found"
  install_bowtie2
fi

if ( checkLocalInstallation bwa )
then
  echo "bwa is found"
else
  echo "bwa is not found"
  install_bwa
fi

if ( checkLocalInstallation samtools )
then
  echo "samtools is found"
else
  echo "samtools is not found"
  install_samtools
fi

if ( checkLocalInstallation nucmer )
then
  echo "nucmer is found"
else
  echo "nucmer is not found"
  install_mummer
fi

if ( checkSystemInstallation wigToBigWig )
then
  echo "wigToBigWig is found"
else
  echo "wigToBigWig is not found, intall wigToBigWig"
  if [[ "$OSTYPE" == "darwin"* ]]
  then
  {
      ln -sf $rootdir/thirdParty/wigToBigWig_mac $rootdir/bin/wigToBigWig
  }
  else
  {
      ln -sf $rootdir/thirdParty/wigToBigWig $rootdir/bin/wigToBigWig
  }
  fi
fi

if ( checkLocalInstallation idba_ud )
then
  echo "idba is found"
else
  echo "idba is not found"
  install_idba
fi

if ( checkLocalInstallation spades.py )
then
  echo "SPAdes is found"
else
  echo "SPAdes is not found"
  install_spades
fi

if [ -x $rootdir/thirdParty/phage_finder_v2.1/bin/phage_finder_v2.1.sh  ]
then
  echo "phage_finder_v2.1 is found"
else
  echo "phage_finder_v2 is not found"
  install_phageFinder
fi

if ( checkLocalInstallation gottcha.pl  )
then
  echo "gottcha.pl  is found"
else
  echo "gottcha.pl  is not found"
  install_gottcha
fi

if ( checkLocalInstallation metaphlan.py  )
then
  echo "metaphlan  is found"
else
  echo "metaphlan  is not found"
  install_metaphlan
fi

if ( checkLocalInstallation primer3_core  )
then
   echo "primer3  is found"
else
   echo "primer3  is not found"
   install_primer3
fi

if ( checkSystemInstallation FastTreeMP )
then
  FastTree_VER=`FastTreeMP  2>&1 | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`;
  if  ( echo $FastTree_VER | awk '{if($1>="2.1.8") exit 0; else exit 1}' )
  then 
    echo "FastTreeMP is found"
  else
    install_FastTree
  fi
else
  echo "FastTreeMP is not found"
  install_FastTree
fi

if ( checkSystemInstallation raxmlHPC-PTHREADS )
then
  echo "RAxML is found"
else
  echo "RAxML is not found"
  install_RAxML
fi

#if [ -f $rootdir/lib/Parallel/ForkManager.pm ]
if ( checkPerlModule Parallel::ForkManager )
then
  echo "Perl Parallel::ForkManager is found"
else
  echo "Perl Parallel::ForkManager is not found"
  install_perl_parallel_forkmanager
fi

#if [ -f $rootdir/lib/Excel/Writer/XLSX.pm ]
if ( checkPerlModule Excel::Writer::XLSX )
then
  echo "Perl Excel::Writer::XLSX is found"
else
  echo "Perl Excel::Writer::XLSX is not found"
  install_perl_excel_writer
fi

#if [ -f $rootdir/lib/Archive/Zip.pm ]
if ( checkPerlModule Archive::Zip )
then
  echo "Perl Archive::Zip is found"
else
  echo "Perl Archive::Zip is not found"
  install_perl_archive_zip
fi

#if [ -f $rootdir/lib/JSON.pm ]
if ( checkPerlModule JSON )
then
  echo "Perl JSON is found"
else
  echo "Perl JSON is not found"
  install_perl_JSON
fi

#if [ -f $rootdir/lib/HTML/Parser.pm ]
if ( checkPerlModule HTML::Parser )
then
  echo "Perl HTML::Parser is found"
else
  echo "Perl HTML::Parser is not found"
  install_perl_html_parser
fi

#if [ -f $rootdir/lib/String/Approx.pm ]
if ( checkPerlModule String::Approx )
then
  echo "Perl String::Approx is found"
else
  echo "Perl String::Approx is not found"
  install_perl_string_approx
fi

#if [ -f $rootdir/lib/PDF/API2.pm ]
if ( checkPerlModule PDF::API2 )
then
  echo "Perl PDF:API2 is found"
else
  echo "Perl PDF:API2 is not found"
  install_perl_pdf_api2
fi

#if [ -f $rootdir/lib/HTML/Template.pm ]
if ( checkPerlModule HTML::Template )
then
  echo "Perl HTML::Template is found"
else
  echo "Perl HTML::Template is not found"
  install_perl_html_template
fi

if ( checkPerlModule Bio::Phylo )
then
  echo "Perl Bio::Phylo is found"
else
  echo "Perl Bio::Phylo is not found"
  install_perl_bio_phylo
fi

if ( checkPerlModule XML::Twig )
then
  echo "Perl XML::Twig is found"
else
  echo "Perl XML::Twig is not found"
  install_perl_xml_twig
fi

if ( checkPerlModule CGI::Session )
then
  echo "Perl CGI::Session is found"
else
  echo "Perl CGI::Session is not found"
  install_perl_cgi_session
fi

if [ -x $rootdir/thirdParty/JBrowse-1.11.6/bin/prepare-refseqs.pl ]
then
  echo "JBrowse is found"
else
  echo "JBrowse is not found"
  install_JBrowse
fi

if [[ "$OSTYPE" == "darwin"* ]]
then
	  ln -sf $rootdir/thirdParty/gottcha/bin/splitrim $rootdir/scripts/microbial_profiling/script/splitrim
      cp -fR $rootdir/start_edge_ui.sh $HOME/Desktop/EDGE_Python_Server_startup
else
	  ln -sf $rootdir/thirdParty/gottcha/bin/splitrim $rootdir/scripts/microbial_profiling/script/splitrim
	  mkdir -p $HOME/Desktop
      sed -e 's,<edge_home>,'"$rootdir"',g'  $rootdir/scripts/EDGE.desktop > $HOME/Desktop/EDGE.desktop
fi

cd $rootdir

mkdir -p $rootdir/edge_ui/data
perl $rootdir/edge_ui/cgi-bin/edge_build_list.pl $rootdir/edge_ui/data/Host/* > $rootdir/edge_ui/data/host_list.json
perl $rootdir/edge_ui/cgi-bin/edge_build_list.pl -sort_by_size -basename $rootdir/database/NCBI_genomes/  > $rootdir/edge_ui/data/Ref_list.json

if [ -f $HOME/.bashrc ]
then
{
  echo "#Added by EDGE pipeline installation" >> $HOME/.bashrc
  echo "export EDGE_HOME=$rootdir" >> $HOME/.bashrc
  echo "export PATH=$rootdir/bin/:$PATH:$rootdir/scripts" >> $HOME/.bashrc
}
else
{
  echo "#Added by EDGE pipeline installation" >> $HOME/.bash_profile
  echo "export EDGE_HOME=$rootdir" >> $HOME/.bash_profile
  echo "export PATH=$rootdir/bin/:$PATH:$rootdir/scripts" >> $HOME/.bash_profile
}
fi

sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/cgi-bin/edge_config.tmpl
sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/apache_conf/edge_apache.conf

TOLCPU=`cat /proc/cpuinfo | grep processor | wc -l`;
if [ $TOLCPU -gt 0 ]
then
{
	sed -i.bak 's,%TOTAL_NUM_CPU%,'"$TOLCPU"',g' $rootdir/edge_ui/cgi-bin/edge_config.tmpl
	DEFAULT_CPU=`echo -n $((TOLCPU/3))`;
	if [ $DEFAULT_CPU -lt 1 ]
	then
	{
		sed -i.bak 's,%DEFAULT_CPU%,'"1"',g' $rootdir/edge_ui/index.html
	}
	else
	{
		sed -i.bak 's,%DEFAULT_CPU%,'"$DEFAULT_CPU"',g' $rootdir/edge_ui/index.html
	}
	fi
}
fi

echo "

All done! Please Restart the Terminal Session.

Run
./runPipeline
for usage.

Read the README
for more information!

Thanks!
"

