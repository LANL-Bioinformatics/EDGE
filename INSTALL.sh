#!/usr/bin/env bash
set -e
rootdir=$( cd $(dirname $0) ; pwd -P )
exec >  >(tee install.log)
exec 2>&1

cd $rootdir
cd thirdParty

mkdir -p $rootdir/bin

export PATH=$PATH:$rootdir/bin/:$rootdir/thirdParty/Anaconda2/bin
export CPLUS_INCLUDE_PATH=$rootdir/thirdParty/Anaconda2/include/:$CPLUS_INCLUDE_PATH

if [ ! -d $HOME ]; then export HOME=$rootdir; fi	

anaconda3bin=$rootdir/thirdParty/Anaconda3/bin
anaconda2bin=$rootdir/thirdParty/Anaconda2/bin


assembly_tools=( idba spades megahit lrasm racon )
annotation_tools=( prokka RATT tRNAscan barrnap BLAST+ blastall phageFinder glimmer aragorn prodigal tbl2asn ShortBRED )
utility_tools=( FaQCs bedtools R GNU_parallel tabix JBrowse bokeh primer3 samtools bcftools sratoolkit ea-utils omics-pathway-viewer NanoPlot Porechop Rpackages )
alignments_tools=( hmmer infernal bowtie2 bwa mummer diamond minimap2 )
taxonomy_tools=( kraken2 metaphlan2 kronatools gottcha gottcha2 centrifuge miccr )
phylogeny_tools=( FastTree RAxML )
perl_modules=( perl_parallel_forkmanager perl_excel_writer perl_archive_zip perl_string_approx perl_pdf_api2 perl_html_template perl_html_parser perl_JSON perl_bio_phylo perl_xml_twig perl_cgi_session perl_email_valid perl_mailtools )
python_packages=( Anaconda2 Anaconda3 )
metagenome_tools=( MaxBin )
pipeline_tools=( DETEQT reference-based_assembly PyPiReT )
all_tools=( "${pipeline_tools[@]}" "${python_packages[@]}" "${assembly_tools[@]}" "${annotation_tools[@]}" "${utility_tools[@]}" "${alignments_tools[@]}" "${taxonomy_tools[@]}" "${phylogeny_tools[@]}" "${metagenome_tools[@]}" "${perl_modules[@]}")

### Install functions ###
install_MaxBin(){
local VER=2.2.5
echo "------------------------------------------------------------------------------
                           Installing MaxBin
------------------------------------------------------------------------------
"
tar xzf FragGeneScan1.31.tar.gz 
cd FragGeneScan1.31
make clean
make fgs
cd $rootdir/thirdParty

tar xzf MaxBin-$VER.tar.gz 
cd MaxBin-$VER
cd src
make
cd ..
cat << _EOM > setting
[FragGeneScan] $rootdir/thirdParty/FragGeneScan1.31
[Bowtie2] $rootdir/bin
[HMMER3] $rootdir/bin
[IDBA_UD] $rootdir/bin
_EOM
ln -sf $rootdir/thirdParty/MaxBin-$VER $rootdir/bin/MaxBin
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           MaxBin installed
------------------------------------------------------------------------------
"

}


install_diamond(){
local VER=0.9.10
echo "------------------------------------------------------------------------------
                           Installing diamond aligner
------------------------------------------------------------------------------
"
tar xvzf diamond-linux64.tar.gz
rm -f diamond_manual.pdf
mv -f diamond $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           diamond aligner installed
------------------------------------------------------------------------------
"
}

install_reference-based_assembly(){
echo "------------------------------------------------------------------------------
                           Installing reference-based_assembly package
------------------------------------------------------------------------------
"
tar xvzf reference-based_assembly.tgz
cd reference-based_assembly
./INSTALL.sh
ln -sf $rootdir/thirdParty/reference-based_assembly $rootdir/bin/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           reference-based_assembly package installed
------------------------------------------------------------------------------
"
}

install_FaQCs(){
local VER=2.09
echo "------------------------------------------------------------------------------
                           Installing FaQCs $VER
------------------------------------------------------------------------------
"
tar xvzf FaQCs-$VER.tar.gz
cd FaQCs
make
cp -f FaQCs $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           FaQCs $VER installed
------------------------------------------------------------------------------
"
}

install_omics-pathway-viewer(){
local VER=0.3
echo "------------------------------------------------------------------------------
                           Installing omics-pathway-viewer $VER
------------------------------------------------------------------------------
"
tar xvzf omics-pathway-viewer.tgz
cp $rootdir/thirdParty/omics-pathway-viewer/scripts/opaver_anno.pl $rootdir/bin/opaver_anno.pl
cp -fR $rootdir/thirdParty/omics-pathway-viewer/opaver_web $rootdir/edge_ui/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           omics-pathway-viewer $VER installed
------------------------------------------------------------------------------
"
}

install_DETEQT(){
local VER=0.3.0
echo "------------------------------------------------------------------------------
                           Installing DETEQT $VER
------------------------------------------------------------------------------
"
tar xvzf DETEQT-$VER.tgz
cd DETEQT
./INSTALL.sh
ln -sf $rootdir/thirdParty/DETEQT $rootdir/bin/DETEQT
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           DETEQT $VER installed
------------------------------------------------------------------------------
"
}

install_PyPiReT(){
local VER=0.3.2
echo "------------------------------------------------------------------------------
                           Installing PyPiReT $VER
------------------------------------------------------------------------------
"
tar xvzf PyPiReT-$VER.tgz
cd PyPiReT
Org_PATH=$PATH;
export PATH=$rootdir/thirdParty/Anaconda3/bin:$rootdir/bin:$PATH;
if [ -e $rootdir/thirdParty/PyPiReT/thirdParty/miniconda/envs/piret ]
then
  rm -rf $rootdir/thirdParty/PyPiReT/thirdParty/miniconda/envs/piret
fi 
./INSTALL.sh
ln -sf $rootdir/thirdParty/PyPiReT $rootdir/bin/PyPiReT
export PATH=$Org_PATH
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           PyPiReT $VER installed
------------------------------------------------------------------------------
"
}

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
./configure --prefix=$rootdir CXXFLAGS='-g -O2 -std=c++03'
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
local VER=3.11.1
echo "------------------------------------------------------------------------------
                           Installing SPAdes $VER
------------------------------------------------------------------------------
"
tar xvzf SPAdes-$VER-Linux.tar.gz 
ln -sf $rootdir/thirdParty/SPAdes-$VER-Linux/bin/spades.py $rootdir/bin/spades.py
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           SPAdes $VER installed
------------------------------------------------------------------------------
"
}

install_megahit(){
local VER=1.1.3
## --version MEGAHIT v1.1.3
echo "------------------------------------------------------------------------------
                           Installing megahit $VER
------------------------------------------------------------------------------
"
tar xvzf megahit-$VER.tar.gz 
cd megahit-$VER
make
cp -f megahit* $rootdir/bin/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           megahit $VER installed
------------------------------------------------------------------------------
"
}

install_racon(){
local VER=1.3.1
echo "------------------------------------------------------------------------------
                           Installing racon $VER
------------------------------------------------------------------------------"
tar xvzf racon-v$VER.tar.gz
cd racon-v$VER
mkdir -p build
cd build
$anaconda2bin/cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp -f bin/racon* $rootdir/bin/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           racon $VER installed
------------------------------------------------------------------------------
"
}

install_lrasm(){
local VER=0.1.0
echo "------------------------------------------------------------------------------
                           Installing long_read_assembly $VER
------------------------------------------------------------------------------
"
tar xvzf long_read_assembly-$VER.tgz
cd long_read_assembly
./INSTALL.sh
ln -sf $rootdir/thirdParty/long_read_assembly/lrasm $rootdir/bin/lrasm
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           long_read_assembly $VER installed
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
local VER=1.13
echo "------------------------------------------------------------------------------
                           Installing prokka-$VER
------------------------------------------------------------------------------
"
tar xvzf prokka-$VER.tar.gz
cd prokka-$VER
# remove old perl modules
# https://github.com/tseemann/prokka/issues/293
rm -rf perl5
cd $rootdir/thirdParty
ln -sf $rootdir/thirdParty/prokka-$VER/bin/prokka $rootdir/bin/prokka
$rootdir/thirdParty/prokka-$VER/bin/prokka --setupdb
echo "
------------------------------------------------------------------------------
                           prokka-$VER installed
------------------------------------------------------------------------------
"
}

install_barrnap()
{
local VER=0.9
echo "------------------------------------------------------------------------------
                           Installing barrnap-$VER
------------------------------------------------------------------------------
"
tar xvzf barrnap-$VER.tar.gz
cd barrnap-$VER
cd $rootdir/thirdParty
ln -sf $rootdir/thirdParty/barrnap-$VER/bin/barrnap $rootdir/bin/barrnap
echo "
------------------------------------------------------------------------------
                           barrnap-$VER installed
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
local VER=2.9.2
echo "------------------------------------------------------------------------------
                           Installing sratoolkit.$VER-linux64
------------------------------------------------------------------------------
"
tar xvzf sratoolkit.$VER-linux64.tgz
cd sratoolkit.$VER-linux64
ln -sf $rootdir/thirdParty/sratoolkit.$VER-linux64/bin/fastq-dump $rootdir/bin/fastq-dump
ln -sf $rootdir/thirdParty/sratoolkit.$VER-linux64/bin/vdb-dump $rootdir/bin/vdb-dump
./bin/vdb-config --restore-defaults
./bin/vdb-config -s /repository/user/default-path=$rootdir/edge_ui/ncbi
./bin/vdb-config -s /repository/user/main/public/root=$rootdir/edge_ui/ncbi/public
if [[ -n ${HTTP_PROXY} ]]; then
	proxy_without_protocol=${HTTP_PROXY#http://}
        ./bin/vdb-config --proxy $proxy_without_protocol
fi
if [[ -n ${http_proxy} ]]; then
	proxy_without_protocol=${http_proxy#http://}
        ./bin/vdb-config --proxy $proxy_without_protocol
fi

ln -sf $HOME/.ncbi $rootdir/.ncbi

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           sratoolkit.$VER-linux64 installed
------------------------------------------------------------------------------
"
}

install_ea-utils(){
echo "------------------------------------------------------------------------------
                           Installing ea-utils.1.1.2-537
------------------------------------------------------------------------------
"
tar xvzf ea-utils.1.1.2-537.tar.gz
cd ea-utils.1.1.2-537
PREFIX=$rootdir make install

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           ea-utils.1.1.2-537 installed
------------------------------------------------------------------------------
"
}

install_R()
{
local VER=3.5.1
echo "------------------------------------------------------------------------------
                           Compiling R $VER
------------------------------------------------------------------------------
"
tar xvzf R-$VER.tar.gz
cd R-$VER
./configure --prefix=$rootdir
make
make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           R compiled
------------------------------------------------------------------------------
"
}
install_Rpackages()
{
echo "------------------------------------------------------------------------------
                           installing R packages
------------------------------------------------------------------------------
"
local VER=3.5.1
tar xzf R-${VER}-Packages.tgz  
echo "if(\"gridExtra\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"gridExtra\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file - 
echo "if(\"devtools\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"devtools\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"phyloseq\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"phyloseq\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"dplyr\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"dplyr\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"Cairo\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"Cairo\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"plotly\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"plotly\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"MetaComp\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"MetaComp\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
echo "if(\"gplots\" %in% rownames(installed.packages()) == FALSE)  {install.packages(c(\"gplots\"), repos = NULL, type=\"source\", contriburl=\"file:Rpackages/\")}" | $rootdir/bin/Rscript --no-init-file  - 
rm -r Rpackages/
echo "
------------------------------------------------------------------------------
                           R packages installed
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
local VER=2.5.0
echo "------------------------------------------------------------------------------
                           Install ncbi-blast-$VER+-x64
------------------------------------------------------------------------------
"
BLAST_ZIP=ncbi-blast-$VER+-x64-linux.tar.gz
if [[ "$OSTYPE" == "darwin"* ]]
then
{
    VER=2.2.29
    BLAST_ZIP=ncbi-blast-$VER+-universal-macosx.tar.gz
}
fi

tar xvzf $BLAST_ZIP
cd ncbi-blast-$VER+
cp -fR bin/* $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           ncbi-blast-$VER+-x64 installed
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

install_kraken2()
{
local VER=2.0.7
echo "------------------------------------------------------------------------------
                           Install kraken-$VER-beta
------------------------------------------------------------------------------
"
tar xvzf kraken-$VER.tgz
cd kraken2-$VER-beta
./install_kraken2.sh $rootdir/bin/

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           kraken-$VER installed
------------------------------------------------------------------------------
"
}

install_centrifuge()
{
local VER=1.0.4
echo "------------------------------------------------------------------------------
                           Install centrifuge-$VER
------------------------------------------------------------------------------
"
tar xvzf centrifuge-$VER-beta.tgz
cd centrifuge-$VER-beta
make 
cp centrifuge $rootdir/bin/
cp centrifuge-* $rootdir/bin/

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           centrifuge-$VER installed
------------------------------------------------------------------------------
"
}

install_JBrowse()
{
local VER=1.12.3
echo "------------------------------------------------------------------------------
                           Installing JBrowse-$VER
------------------------------------------------------------------------------
"
tar xvzf JBrowse-$VER.tar.gz
if [ -e $rootdir/edge_ui/JBrowse/data ]
then
  mv $rootdir/edge_ui/JBrowse/data $rootdir/edge_ui/JBrowse_olddata
fi
if [ -e $rootdir/edge_ui/JBrowse ]
then
  rm -rf $rootdir/edge_ui/JBrowse
fi

mv JBrowse-$VER $rootdir/edge_ui/JBrowse
cd $rootdir/edge_ui/JBrowse
./setup.sh
if [ -e $rootdir/edge_ui/JBrowse_olddata ]
then
  mv $rootdir/edge_ui/JBrowse_olddata $rootdir/edge_ui/JBrowse/data
else
  mkdir -p -m 775 data
fi

cd $rootdir/thirdParty
#ln -sf $rootdir/thirdParty/JBrowse-1.11.6 $rootdir/edge_ui/JBrowse
echo "
------------------------------------------------------------------------------
                           JBrowse-$VER installed
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
                           Compiling bowtie2 2.2.6
------------------------------------------------------------------------------
"
tar xvzf bowtie2-2.2.6.tar.gz
cd bowtie2-2.2.6
make
cp bowtie2* $rootdir/bin/.
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
tar xzf gottcha.tar.gz
rm -rf gottcha/ext/opt/dmd2/
tar xzf dmd.2.082.0.linux.tgz -C gottcha/ext/opt
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

install_gottcha2()
{
echo "------------------------------------------------------------------------------
                           Compiling gottcha-2.1 BETA
------------------------------------------------------------------------------
"
tar xvzf gottcha2.tar.gz
cd gottcha2
#ln -sf $PWD/gottcha.py $rootdir/bin/
ln -sf $rootdir/database/GOTTCHA2 ./database
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           gottcha-2.1 BETA installed
------------------------------------------------------------------------------
"
}

install_miccr()
{
local VER=0.0.2
echo "------------------------------------------------------------------------------
                           Installing miccr-$VER
------------------------------------------------------------------------------
"
tar xvzf miccr-$VER.tgz
cd miccr
ln -sf $rootdir/database/miccrDB ./database
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           miccr-$VER installed
------------------------------------------------------------------------------
"
}


install_pangia()
{
local VER=2.4.11
echo "------------------------------------------------------------------------------
                           Installing PANGIA $VER BETA
------------------------------------------------------------------------------
"
if [ -e $rootdir/thirdParty/pangia/pangia-vis/data ]
then
  mv $rootdir/thirdParty/pangia/pangia-vis/data $rootdir/thirdParty/pangia-vis-data
fi
if [ -e $rootdir/thirdParty/pangia ]
then
  rm -rf $rootdir/thirdParty/pangia
fi

tar xvzf pangia-$VER.tar.gz
cd pangia
if [ -e $rootdir/thirdParty/pangia-vis-data ]
then
  cp -f $rootdir/thirdParty/pangia-vis-data/* $rootdir/thirdParty/pangia/pangia-vis/data/
rm -rf $rootdir/thirdParty/pangia-vis-data
fi 

cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           pangia-$VER BETA
------------------------------------------------------------------------------
"
}


install_metaphlan2()
{
local VER=2.7.7
echo "------------------------------------------------------------------------------
                           Installing metaphlan-$VER
------------------------------------------------------------------------------
"
tar xvzf metaphlan-$VER.tgz
cd metaphlan-$VER
cp -fR metaphlan2.py $rootdir/bin/.
cp -fR utils/read_fastx.py $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           metaphlan-$VER installed
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

install_ShortBRED()
{
echo "------------------------------------------------------------------------------
                           Installing ShortBRED
------------------------------------------------------------------------------
"
tar xvzf ShortBRED-0.9.4M.tgz
ln -sf $rootdir/thirdParty/ShortBRED-0.9.4M $rootdir/bin/ShortBRED
echo "
------------------------------------------------------------------------------
                           ShortBRED installed
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

install_minimap2()
{
local VER=2.10
echo "------------------------------------------------------------------------------
                           Compiling minimap2 $VER
------------------------------------------------------------------------------
"
tar xvzf minimap2-${VER}_x64-linux.tgz
cd minimap2-${VER}_x64-linux
cp minimap2 $rootdir/bin/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           minimap2 compiled
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
cp mummerplot $rootdir/bin/.
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
make CC_OPTS='-fpermissive -g -Wall -D__USE_FIXED_PROTOTYPES__'
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
local VER=2.7
echo "------------------------------------------------------------------------------
               Installing KronaTools-$VER
------------------------------------------------------------------------------
"
tar xvzf KronaTools-$VER.tar.gz
cd KronaTools-$VER
perl install.pl --prefix $rootdir --taxonomy $rootdir/database/Krona_taxonomy
#./updateTaxonomy.sh --local
cp $rootdir/scripts/microbial_profiling/script/ImportBWA.pl scripts/
ln -sf $rootdir/thirdParty/KronaTools-$VER/scripts/ImportBWA.pl $rootdir/bin/ktImportBWA 
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        KronaTools-$VER Installed
------------------------------------------------------------------------------
"
}

install_samtools()
{
local VER=1.9
echo "------------------------------------------------------------------------------
                           Compiling samtools-$VER
------------------------------------------------------------------------------
"
tar xvzf samtools-$VER.tar.gz
cd samtools-$VER
#make CFLAGS='-g -fPIC -Wall -O2'
./configure --prefix=$rootdir
make
make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           samtools $VER compiled
------------------------------------------------------------------------------
"
}

install_bcftools()
{
local VER=1.9
echo "------------------------------------------------------------------------------
                           Compiling bcftools-$VER
------------------------------------------------------------------------------
"
tar xvzf bcftools-$VER.tar.gz
cd bcftools-$VER
./configure --prefix=$rootdir
make
make install
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                           bcftools $VER compiled
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
local VER=1.19
echo "------------------------------------------------------------------------------
               Installing Perl Module Parallel-ForkManager-$VER
------------------------------------------------------------------------------
"
tar xvzf Parallel-ForkManager-$VER.tar.gz
cd Parallel-ForkManager-$VER
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
                        Parallel-ForkManager-$VER Installed
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

install_perl_email_valid(){
local VER=1.202
echo "-----------------------------------------------------------------------------
		Installing Perl Module Email-Valid-$VER
 ------------------------------------------------------------------------------
"
tar xvzf Email-Valid-$VER.tar.gz
cd Email-Valid-$VER
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
		Email-Valid-$VER Installed
------------------------------------------------------------------------------
"
}

install_perl_mailtools(){
local VER=2.20
echo "-----------------------------------------------------------------------------
		Installing Perl Module MailTools-$VER
 ------------------------------------------------------------------------------
"
tar xvzf MailTools-$VER.tar.gz
cd MailTools-$VER
perl Makefile.PL
make
cp -fR blib/lib/* $rootdir/lib/.
cd $rootdir/thirdParty
echo "
------------------------------------------------------------------------------
		MailTools-$VER Installed
------------------------------------------------------------------------------
"
}

install_Anaconda2()
{
echo "------------------------------------------------------------------------------
                 Installing Python Anaconda2 4.1.1
------------------------------------------------------------------------------
"
if [ ! -f $rootdir/thirdParty/Anaconda2/bin/python ]; then
    bash Anaconda2-4.1.1-Linux-x86_64.sh -b -p $rootdir/thirdParty/Anaconda2/
fi
ln -fs $anaconda2bin/python $rootdir/bin
ln -fs $anaconda2bin/pip $rootdir/bin
ln -fs $anaconda2bin/conda $rootdir/bin
tar -xvzf Anaconda2Packages.tgz
$anaconda2bin/conda install Anaconda2Packages/biopython-1.68-np111py27_0.tar.bz2 
$anaconda2bin/conda install Anaconda2Packages/blast-2.5.0-boost1.60_1.tar.bz2 
$anaconda2bin/conda install Anaconda2Packages/icu-58.1-0.tar.bz2 
$anaconda2bin/conda install Anaconda2Packages/libgcc-5.2.0-0.tar.bz2 
$anaconda2bin/conda install Anaconda2Packages/mysql-connector-python-2.0.4-py27_0.tar.bz2 
$anaconda2bin/conda install Anaconda2Packages/prodigal-2.60-1.tar.bz2 
#$anaconda2bin/conda install Anaconda2Packages/rgi-3.1.1-py27_1.tar.bz2
$anaconda2bin/conda install Anaconda2Packages/subprocess32-3.2.7-py27_0.tar.bz2
$anaconda2bin/conda install Anaconda2Packages/cmake-3.6.3-0.tar.bz2
$anaconda2bin/conda install Anaconda2Packages/matplotlib-2.0.0-np111py27_0.tar.bz2
$anaconda2bin/conda install Anaconda2Packages/rapsearch-2.24-1.tar.bz2
ln -s $anaconda2bin/rapsearch $anaconda2bin/rapsearch2
$anaconda2bin/pip install --no-index --find-links=./Anaconda2Packages qiime
$anaconda2bin/pip install --no-index --find-links=./Anaconda2Packages xlsx2csv
$anaconda2bin/pip install --no-index --find-links=./Anaconda2Packages h5py
matplotlibrc=`$anaconda2bin/python -c 'import matplotlib as m; print m.matplotlib_fname()' 2>&1`
perl -i.orig -nle 's/(backend\s+:\s+\w+)/\#${1}\nbackend : Agg/; print;' $matplotlibrc
rm -r Anaconda2Packages/
echo "
------------------------------------------------------------------------------
                         Python Anaconda2 4.1.1 Installed
------------------------------------------------------------------------------
"
}

install_Anaconda3()
{
local VER=5.1.0
echo "------------------------------------------------------------------------------
                 Installing Python Anaconda3 $VER
------------------------------------------------------------------------------
"
if [ ! -f $rootdir/thirdParty/Anaconda3/bin/python3 ]; then
    bash Anaconda3-$VER-Linux-x86_64.sh -b -p $rootdir/thirdParty/Anaconda3/
fi
ln -fs $anaconda3bin/python3 $rootdir/bin

tar -xvzf Anaconda3Packages.tgz
$anaconda3bin/pip install --no-index --find-links=./Anaconda3Packages CairoSVG 
$anaconda3bin/pip install --no-index --find-links=./Anaconda3Packages pymc3
$anaconda3bin/pip install --no-index --find-links=./Anaconda3Packages lzstring
$anaconda3bin/conda config --add channels defaults
$anaconda3bin/conda config --add channels bioconda
$anaconda3bin/conda config --add channels conda-forge
$anaconda3bin/conda install -c bioconda rgi=4.2.2
$anaconda3bin/conda install -c conda-forge pandas
ln -fs $anaconda3bin/cairosvg $rootdir/bin

echo "
------------------------------------------------------------------------------
                         Python Anaconda3 $VER Installed
------------------------------------------------------------------------------
"
}

install_bokeh()
{
local VER=0.12.10
echo "------------------------------------------------------------------------------
                        Installing bokeh $VER
------------------------------------------------------------------------------
"
$anaconda3bin/pip install  --no-index --find-links=./Anaconda3Packages bokeh==$VER
echo "
------------------------------------------------------------------------------
                         bokeh $VER Installed
------------------------------------------------------------------------------
"
}

install_NanoPlot()
{
echo "------------------------------------------------------------------------------
                 	Installing NanoPlot
------------------------------------------------------------------------------
"
$anaconda3bin/pip install --no-index --find-links=./Anaconda3Packages NanoPlot
ln -fs $anaconda3bin/NanoPlot $rootdir/bin
echo "
------------------------------------------------------------------------------
                         NanoPlot Installed
------------------------------------------------------------------------------
"
}
install_Porechop()
{
echo "------------------------------------------------------------------------------
                 	Installing Porechop
------------------------------------------------------------------------------
"
$anaconda3bin/conda install Anaconda3Packages/porechop-0.2.3_seqan2.1.1-py36_2.tar.bz2
ln -fs $anaconda3bin/porechop $rootdir/bin
echo "
------------------------------------------------------------------------------
                         Porechop Installed
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
   echo -e "\nMetagenome"
   for i in "${metagenome_tools[@]}"
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
   echo -e "\nPython_Packages"
   for i in "${python_packages[@]}"
   do
           echo "* $i"
   done
   echo -e "\nPipeline_Tools"
   for i in "${pipeline_tools[@]}"
   do
           echo "* $i"
   done
}


### Main ####
if ( checkSystemInstallation csh )
then
  #echo "csh is found"
  echo -n ""
else
  echo "csh is not found"
  echo "Please Install csh first, then INSTALL the package"
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
      Metagenome)
        for tool in "${metagenome_tools[@]}"
        do
            install_$tool
        done
        echo -e "Metagenome tools installed.\n"
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
        exit 0 ;;
      Python_Packages)
        for tool in "${python_packages[@]}"
        do
            install_$tool
        done
        echo -e "Python_Packages installed.\n"
        exit 0 ;;
      Pipeline_Tools)
        for tool in "${pipeline_tools[@]}"
        do
            install_$tool
        done
        echo -e "Pipeline_Tools installed.\n"
        exit 0 ;;
      force)
        for tool in "${all_tools[@]}"
        do
            install_$tool
        done
        ;;
      *)
        if ( containsElement "$f" "${assembly_tools[@]}" || containsElement "$f" "${annotation_tools[@]}" || containsElement "$f" "${alignments_tools[@]}" || containsElement "$f" "${taxonomy_tools[@]}" || containsElement "$f" "${phylogeny_tools[@]}" || containsElement "$f" "${metagenome_tools[@]}" || containsElement "$f" "${utility_tools[@]}" || containsElement "$f" "${perl_modules[@]}" || containsElement "$f" "${python_packages[@]}" || containsElement "$f" "${pipeline_tools[@]}" )
        then
            install_$f
        else
            echo "$f: no this tool in the list"
            print_tools_list
        fi
        exit;;
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

if perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"' >/dev/null 2>&1 
then 
  perl -MBio::Root::Version -e 'print "BioPerl Version ", $Bio::Root::Version::VERSION," is found\n"'
else 
  echo "Cannot find a perl Bioperl Module installed" 1>&2
  echo "Please install Bioperl (http://www.bioperl.org/)"
  exit 1
fi

if $rootdir/bin/python -c 'import Bio; print Bio.__version__' >/dev/null 2>&1
then
  $rootdir/bin/python -c 'import Bio; print "BioPython Version", Bio.__version__, "is found"'
else
  install_Anaconda2
fi

if $rootdir/bin/python3 -c 'import sys; sys.exit("Python > 3.0 required.") if sys.version_info < ( 3, 0) else ""' >/dev/null 2>&1
then
  $rootdir/bin/python3 -c 'import sys; print( "Python3 version %s.%s found." % (sys.version_info[0],sys.version_info[1]))'
else
  install_Anaconda3
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
	R_VER=`$rootdir/bin/R --version | perl -nle 'print $& if m{version \d+\.\d+}'`;
	if  ( echo $R_VER | awk '{if($2>="3.5") exit 0; else exit 1}' )
	then
	{
        	echo "R $R_VER found"
	}
	else
	{
		install_R
	}
	fi
    }
    else
    {
        install_R
    }
    fi
}
fi

install_Rpackages

if ( checkSystemInstallation bedtools )
then
  echo "bedtools is found"
else
  echo "bedtools is not found"
  install_bedtools 
fi

if ( checkSystemInstallation FaQCs )
then
  FaQCs_VER=`FaQCs --version 2>&1| awk '{print $2}'`;
  if  ( echo $FaQCs_VER | awk '{if($1>="2.08") exit 0; else exit 1}' )
  then
    echo "FaQCs $FaQCs_VER found"
  else
    install_FaQCs
  fi
else
  echo "FaQCs is not found"
  install_FaQCs
fi

if ( checkSystemInstallation fastq-dump )
then
  sratoolkit_VER=`fastq-dump --version | perl -nle 'print $& if m{\d\.\d\.\d}'`;
  if  ( echo $sratoolkit_VER | awk '{if($1>="2.8.1") exit 0; else exit 1}' )
  then
    echo "sratoolkit $sratoolkit_VER found"
  else
    install_sratoolkit
  fi
else
  echo "sratoolkit is not found"
  install_sratoolkit
fi

if ( checkSystemInstallation fastq-join )
then
  echo "fastq-join is found"
else
  install_ea-utils
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
   BLAST_VER=`blastn -version | grep blastn | perl -nle 'print $& if m{\d\.\d\.\d}'`;
   if ( echo $BLAST_VER | awk '{if($1>="2.4.0") exit 0; else exit 1}' )
   then
     echo "BLAST+ $BLAST_VER found"
   else
     install_BLAST+
   fi
else
  echo "BLAST+ is not found"
  install_BLAST+
fi

if ( checkLocalInstallation blastall )
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
  Krona_VER=`$rootdir/bin/ktGetLibPath | perl -nle 'print $& if m{KronaTools-\d\.\d}' | perl -nle 'print $& if m{\d\.\d}'`;
  if  ( echo $Krona_VER | awk '{if($1>="2.6") exit 0; else exit 1}' )
  then
    echo "KronaTools $Krona_VER found"
  else
    install_kronatools
  fi
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

if ( checkLocalInstallation ShortBRED/shortbred_quantify.py )
then
  echo "ShortBRED is found"
else
  echo "ShortBRED is not found"
  install_ShortBRED
fi


if ( checkLocalInstallation kraken2 )
then
  kraken_VER=`kraken2 --version | grep version | perl -nle 'print $1 if m{(\d+\.\d+)}'`;
  if  ( echo $kraken2_VER | awk '{if($1>=2.0) exit 0; else exit 1}' )
  then
    echo "kraken2 $kraken2_VER found"
  else
    install_kraken2
  fi
else
  echo "kraken2 is not found"
  install_kraken2
fi

if ( checkLocalInstallation centrifuge )
then
  centrifuge_VER=`centrifuge --version | grep "centrifuge-class version" | perl -nle 'print $1 if m{(\d+\.\d+\.\d+)}'`;
  if  ( echo $centrifuge_VER | awk '{if($1>=1.0.4) exit 0; else exit 1}' )
  then
    echo "centrifuge $centrifuge_VER found"
  else
    install_centrifuge
  fi
else
  echo "centrifuge is not found"
  install_centrifuge
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
  bowtie_VER=`bowtie2 --version | grep bowtie | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`;
  if  ( echo $bowtie_VER | awk '{if($2>="2.2.4") exit 0; else exit 1}' )
  then 
    echo "bowtie2 $bowtie_VER found"
  else
    install_bowtie2
  fi
else
  echo "bowtie2 is not found"
  install_bowtie2
fi

if ( checkSystemInstallation diamond )
then
  diamond_VER=`diamond --version 2>&1| perl -nle 'print $1 if m{(\d+\.\d+.\d+)}'`;
  if  ( echo $diamond_VER | awk '{if($1>="0.9.10") exit 0; else exit 1}' )
  then
    echo "diamond $diamond_VER found"
  else
    install_diamond
  fi
else
  echo "diamond is not found"
  install_diamond
fi

if ( checkSystemInstallation minimap2 )
then
  minimap2_VER=`minimap2 --version 2>&1| perl -nle 'print $1 if m{(\d+\.\d+)}'`;
  if  ( echo $minimap2_VER | awk '{if($1>="2.10") exit 0; else exit 1}' )
  then
    echo "minimap2 $minimap2_VER found"
  else
    install_minimap2
  fi
else
  echo "minimap2 is not found"
  install_minimap2
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
  samtools_installed_VER=`samtools 2>&1| grep 'Version'|perl -nle 'print $1 if m{Version: (\d+\.\d+.\d+)}'`;
  if [ -z "$samtools_installed_VER" ]
  then
      samtools_installed_VER=`samtools 2>&1| grep 'Version'|perl -nle 'print $1 if m{Version: (\d+\.\d+)}'`; 
  fi
  if  ( echo $samtools_installed_VER | awk '{if($1>=1.7) exit 0; else exit 1}' )
  then
      echo "samtools is found"
  else
      echo "samtools is not found"
      install_samtools
  fi
else
  install_samtools
fi

if ( checkLocalInstallation bcftools )
then
  bcftools_installed_VER=`bcftools 2>&1| grep 'Version'|perl -nle 'print $1 if m{Version: (\d+\.\d+.\d+)}'`;
  if [ -z "$bcftools_installed_VER" ]
  then
      bcftools_installed_VER=`bcftools 2>&1| grep 'Version'|perl -nle 'print $1 if m{Version: (\d+\.\d+)}'`;
  fi
  if  ( echo $bcftools_installed_VER | awk '{if($1>=1.7) exit 0; else exit 1}' )
  then
      echo "bcftools is found"
  else
      echo "bcftools is not found"
      install_bcftools
  fi
else
  install_bcftools
fi

if ( checkLocalInstallation nucmer )
then
  echo "nucmer is found"
else
  echo "nucmer is not found"
  install_mummer
fi

if ( checkLocalInstallation wigToBigWig )
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

if ( checkSystemInstallation spades.py )
then
  spades_VER=`spades.py 2>&1 | perl -nle 'print $& if m{\d\.\d\.\d}'`;
  if ( echo $spades_VER | awk '{if($1>="3.9.0") exit 0; else exit 1}' )
  then
    echo "SPAdes $spades_VER found"
  else
    install_spades
  fi
else
  echo "SPAdes is not found"
  install_spades
fi

if ( checkSystemInstallation megahit  )
then
  ## --version MEGAHIT v1.1.3
  megahit_VER=`megahit --version | perl -nle 'print $& if m{\d\.\d.\d}'`;
  if  ( echo $megahit_VER | awk '{if($1>="1.1.3") exit 0; else exit 1}' )
  then
    echo "megahit $megahit_VER found"
  else
    install_megahit
  fi
else
  echo "megahit is not found"
  install_megahit
fi

if ( checkSystemInstallation lrasm  )
then
  racon_installed_VER=`racon --version | perl -nle 'print $1 if m{v(\d+\.\d+\.*\d*)}'`;
  if ( echo $racon_installed_VER | awk '{if($1>=1.3.1) exit 0; else exit 1}' )
  then
    echo "racon $racon_installed_VER found"
  else
    install_racon
  fi
else
  echo "racon is not found"
  install_racon
fi

if ( checkSystemInstallation lrasm  )
then
  lrasm_VER=`lrasm --version | perl -nle 'print $& if m{\d\.\d.\d}'`;
  if  ( echo $lrasm_VER | awk '{if($1>="0.1.0") exit 0; else exit 1}' )
  then
    echo "lrasm $lrasm_VER found"
  else
    install_lrasm
  fi
else
  echo "lrasm is not found"
  install_lrasm
fi

if [ -x $rootdir/thirdParty/phage_finder_v2.1/bin/phage_finder_v2.1.sh  ]
then
  echo "phage_finder_v2.1 is found"
else
  echo "phage_finder_v2 is not found"
  install_phageFinder
fi

if [ -x $rootdir/bin/MaxBin/src/MaxBin ]
then
  MaxBin_VER=`$rootdir/bin/MaxBin/run_MaxBin.pl -v | head -1 |perl -nle 'print $& if m{\d\.\d}'`;
  if ( echo $MaxBin_VER | awk '{if($1>="2.2") exit 0; else exit 1}' )
  then
    echo "MaxBin2 is found"
  else
    install_MaxBin
  fi
else
  echo "MaxBin2 is not found"
  install_MaxBin
fi

if ( checkLocalInstallation gottcha.pl  )
then
  echo "gottcha.pl  is found"
else
  echo "gottcha.pl  is not found"
  install_gottcha
fi

if [ -x $rootdir/thirParty/gottcha2/gottcha.py ]
then
  gottcha2_VER=`$rootdir/thirParty/gottcha2/gottcha.py -h | grep VERSION |perl -nle 'print $& if m{\d\.\d}'`;
  if ( echo $gottcha2_VER | awk '{if($1>="2.1") exit 0; else exit 1}' )
  then
    echo "GOTTCHA2 $gottcha2_VER is found"
  else
    install_gottcha2
  fi
else
  echo "GOTTCHA2 is not found"
  install_gottcha2
fi

if [ -x $rootdir/thirParty/miccr/miccr.py ]
then
  miccr_VER=`$rootdir/thirParty/miccr/miccr.py -h | grep MICCR |perl -nle 'print $& if m{\d\.\d.\d}'`;
  if ( echo $miccr_VER | awk '{if($1>="0.0.2") exit 0; else exit 1}' )
  then
    echo "miccr $miccr_VER is found"
  else
    install_miccr
  fi
else
  echo "miccr is not found"
  install_miccr
fi

#if [ -x $rootdir/thirParty/pangia/pangia.py ]
#then
#  pangia_VER=`$rootdir/thirParty/pangia/pangia.py -h | grep 'PanGIA Bioinformatics' |perl -nle 'print $& if m{\d\.\d\.\d}'`;
#  if ( echo $pangia_VER | awk '{if($1>="2.4.5") exit 0; else exit 1}' )
#  then
#    echo "PANGIA $pangia_VER is found"
#  else
#    install_pangia
#  fi
#else
#  echo "PANGIA is not found"
#  install_pangia
#fi

if ( checkLocalInstallation metaphlan2.py  )
then
  metaphlan_VER=`metaphlan2.py -v 2>&1 |perl -nle 'print $& if m{\d\.\d\.\d}'`;
  if ( echo $metaphlan_VER | awk '{if($1>="2.7.7") exit 0; else exit 1}' )
  then
    echo "metaphlan2  is found"
  else
    install_metaphlan2
  fi
else
  echo "metaphlan2 is not found"
  install_metaphlan2
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
  if  ( echo $FastTree_VER | awk '{if($2>="2.1.8") exit 0; else exit 1}' )
  then
    echo "FastTreeMP is found"
  else
   install_FastTree
  fi
else
  echo "FastTreeMP is not found"
  install_FastTree
fi

if ( checkLocalInstallation reference-based_assembly )
then
    echo "reference-based_assembly is found"
else
    echo "reference-based_assembly is not found"
    install_reference-based_assembly
fi

if ( checkLocalInstallation porechop )
then
    echo "porechop is found"
else
    echo "porechop is not found"
    install_Porechop
fi

if ( checkLocalInstallation NanoPlot )
then
    echo "NanoPlot is found"
else
    echo "NanoPlot is not found"
    install_NanoPlot
fi

if [ -x $anaconda3bin/bokeh ]
then
  bokeh_VER=`$anaconda3bin/bokeh --version |perl -nle 'print $& if m{\d\.\d+\.\d+}'`;
  if ( echo $bokeh_VER | awk '{if($1=="0.12.10") exit 0; else exit 1}' )
  then
    echo "bokeh $bokeh_VER is found"
  else
    install_bokeh
  fi
else
  echo "bokeh is not found"
  install_bokeh
fi

if ( checkLocalInstallation DETEQT )
then
  DETEQT_VER=`$rootdir/bin/DETEQT/DETEQT -V | perl -nle 'print $& if m{Version \d+\.\d+\.\d+}'`;
  if  ( echo $DETEQT_VER | awk '{if($2>="0.3.0") exit 0; else exit 1}' )
  then
    echo "DETEQT is found"
  else
   install_DETEQT
  fi
else
  echo "DETEQT is not found"
  install_DETEQT
fi


if ( checkLocalInstallation PyPiReT)
then
  PiReT_VER=`$rootdir/bin/PyPiReT/bin/runPiReT -V | perl -nle 'print $& if m{Version \d+\.\d+\.*\d*}'`;
  if  ( echo $PiReT_VER | awk '{if($2>="0.3.1") exit 0; else exit 1}' )
  then
    echo "PyPiReT is found"
  else
   install_PyPiReT
  fi
else
  echo "PyPiReT is not found"
  install_PyPiReT
fi

if ( checkLocalInstallation opaver_anno.pl )
then
  opaver_VER=`opaver_anno.pl -h | perl -nle 'print $1 if m{Version: v(\d+\.\d+)}'`;
  if  ( echo $opaver_VER | awk '{if($1>="0.3") exit 0; else exit 1}' )
  then
    echo "omics-pathwar-viewer is found"
  else
   install_omics-pathway-viewer
  fi
else
  echo "omics-pathway-viewer is not found"
  install_omics-pathway-viewer
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
  Parallel_ForManager_installed_VER=`perl -e "use lib '$rootdir/lib'; use Parallel::ForkManager; print \\\$Parallel::ForkManager::VERSION;"`
  if  ( echo $Parallel_ForManager_installed_VER | awk '{if($1>="1.03") exit 0; else exit 1}' )
  then
    echo "Perl Parallel::ForkManager $Parallel_ForManager_installed_VER is found"
  else 
    install_perl_parallel_forkmanager
  fi
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

if ( checkPerlModule Email::Valid )
then
  echo "Perl Email::Valid is found"
else
  echo "Perl Email::Valid is not found"
  install_perl_email_valid
fi

if ( checkPerlModule Mail::Address )
then
  echo "Perl Mail::Address is found"
else
  echo "Perl Mail::Address is not found"
  install_perl_mailtools
fi

if [ -x $rootdir/edge_ui/JBrowse/bin/prepare-refseqs.pl ]
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
#perl $rootdir/edge_ui/cgi-bin/edge_build_list.pl -sort_by_size -basename $rootdir/database/NCBI_genomes/  > $rootdir/edge_ui/data/Ref_list.json

echo "Setting up EDGE_input"
if [ -d $rootdir/edge_ui/EDGE_input/ ]
then
    rsync -a $rootdir/deployment/public $rootdir/edge_ui/EDGE_input/
    ln -sf $rootdir/testData $rootdir/edge_ui/EDGE_input/public/data/
else
	mkdir -p $HOME/EDGE_input
	rm -rf $rootdir/edge_ui/EDGE_input
	ln -sf $HOME/EDGE_input $rootdir/edge_ui/EDGE_input
	rsync -a $rootdir/deployment/public $rootdir/edge_ui/EDGE_input/
   	ln -sf $rootdir/testData $rootdir/edge_ui/EDGE_input/public/data/
fi
if [ ! -d $rootdir/edge_ui/EDGE_output/ ]
then
   	echo "Setting up EDGE_output"
   	mkdir -p $HOME/EDGE_output
	rm -rf $rootdir/edge_ui/EDGE_output
	ln -sf $HOME/EDGE_output $rootdir/edge_ui/EDGE_output
fi
if [ ! -d $rootdir/edge_ui/EDGE_report/ ]
then
	echo "Setting up EDGE_report"
	mkdir -p $HOME/EDGE_report
	rm -rf $rootdir/edge_ui/EDGE_report
	ln -sf $HOME/EDGE_report $rootdir/edge_ui/EDGE_report
fi

if [ -f $HOME/.bashrc ]
then
{
  echo "#Added by EDGE pipeline installation" >> $HOME/.bashrc
  echo "export EDGE_HOME=$rootdir" >> $HOME/.bashrc
  echo "export EDGE_PATH=$rootdir/bin/:$rootdir/scripts" >> $HOME/.bashrc
  echo "export PATH=\$EDGE_PATH:\$PATH:" >> $HOME/.bashrc
}
else
{
  echo "#Added by EDGE pipeline installation" >> $HOME/.bash_profile
  echo "export EDGE_HOME=$rootdir" >> $HOME/.bash_profile
  echo "export EDGE_PATH=$rootdir/bin/:$rootdir/scripts" >> $HOME/.bashrc
  echo "export PATH=\$EDGE_PATH:\$PATH:" >> $HOME/.bashrc
}
fi
# Check for existing sys.properties file before regenerating a new one from original
if [ -e $rootdir/edge_ui/sys.properties ] 
then
{
	cp $rootdir/edge_ui/sys.properties $rootdir/edge_ui/sys.properties.bak
	cp $rootdir/edge_ui/sys.properties.original $rootdir/edge_ui/sys.properties
}
else
{
	cp $rootdir/edge_ui/sys.properties.original $rootdir/edge_ui/sys.properties
}
fi

sed -i 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/sys.properties
sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/apache_conf/edge_apache.conf
sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/apache_conf/edge_httpd.conf
#sed -i.bak 's,%EDGE_HOME%,'"$rootdir"',g' $rootdir/edge_ui/apache_conf/pangia-vis.conf

TOLCPU=`cat /proc/cpuinfo | grep processor | wc -l`;
if [ $TOLCPU -gt 0 ]
then
{
	sed -i 's,%TOTAL_NUM_CPU%,'"$TOLCPU"',g' $rootdir/edge_ui/sys.properties
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



## Cleanup
rm -r $rootdir/thirdParty/Anaconda3Packages/
$anaconda2bin/conda clean -y -a
$anaconda3bin/conda clean -y -a

# set up a cronjob for project old files clena up
echo "01 00 * * * perl $rootdir/edge_ui/cgi-bin/edge_data_cleanup.pl" | crontab -
(crontab -l ; echo "* * * * * perl $rootdir/edge_ui/cgi-bin/edge_auto_run.pl > /dev/null 2>&1") | crontab -
#(crontab -l ; echo "*/1 * * * * bash $rootdir/scripts/pangia-vis-checker.sh > /dev/null 2>&1") | crontab -

echo "

All done! Please Restart the Terminal Session.

Run
./runPipeline
for usage.

Read the README
for more information!

Thanks!
"

