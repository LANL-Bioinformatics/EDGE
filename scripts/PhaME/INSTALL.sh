#!/bin/bash

rootdir=$( cd $(dirname $0) ; pwd -P )

export PATH=$rootdir/ext/bin/:$PATH;
export PERL5LIB=$rootdir/ext/lib/:$PERL5LIB;
mkdir -p bin;

done_message () {
   if [ $? == 0 ]; then
      if [ "x$1" != "x" ]; then
          echo $1;
      else
          echo "done.";
      fi
   else
      echo "Installation failed." $2
      exit 1;
   fi  
}

download_ext () {
   if [ -e $2 ]; then
      echo "$2 existed. Skiping downloading $1."
      return;
   fi; 
  
   if hash curl 2>/dev/null; then
		if [[ -n ${HTTP_PROXY} ]]; then
   	  		curl --proxy $HTTP_PROXY -L $1  -o $2;
   	  	else
			curl -L $1 -o $2; 
		fi;
   else
      wget -O $2 $1; 
   fi; 

   if [ ! -r $2 ]; then
      echo "ERROR: $1 download failed."
   fi; 
}
 
if ! hash unzip 2>/dev/null; then
   echo "WARNNING: Unzip not found. You might experience error in further installation.";
fi;

echo "Checking NUCMER ..."

NUCMER_VER=`nucmer --version 2>&1 | perl -nle 'print $& if m{version \d+\.\d+}' | perl -nle 'print $& if m{\d+\.\d+}'`;

if ( hash nucmer 2>/dev/null ) && ( echo $NUCMER_VER | awk '{if($1>="3.07") exit 0; else exit 1}' )

then
   echo "NUCmer >=3.07 found.";
else
   echo "NUCmer >=3.07 or above not found/. Trying to download from https://github.com/chienchi/MUMmer/archive/master.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/chienchi/MUMmer/archive/master.zip ext/opt/nucmer.zip;
   unzip ext/opt/nucmer.zip -d ext/opt/;
   cd ext/opt/MUMmer-master;
   make CPPFLAGS="-O3 -DSIXTYFOURBITS";
   cp nucmer $rootdir/ext/bin/.
   cp show-coords $rootdir/ext/bin/.
   cp show-snps $rootdir/ext/bin/.
   cp mgaps $rootdir/ext/bin/.
   cp delta-filter $rootdir/ext/bin/.	
   cd $rootdir;
fi;
done_message " Done." "";

echo "Checking BWA ..."

BWA_VER=`bwa 2>&1 | perl -nle 'print $& if m{Version: \d+\.\d+}' | perl -nle 'print $& if m{\d+\.\d+}'`;

if ( hash bwa 2>/dev/null ) && ( echo $BWA_VER | awk '{if($1>="0.7") exit 0; else exit 1}' )
then
   echo "BWA >=0.7 found.";
else
   echo "BWA >=0.7 or above not found. Trying to download from https://github.com/lh3/bwa/archive/master.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/lh3/bwa/archive/master.zip ext/opt/bwa.zip;
   unzip ext/opt/bwa.zip -d ext/opt/;
   cd ext/opt/bwa-master;
   make;
   cd $rootdir;
   cp ext/opt/bwa-master/bwa ext/bin/;
fi;
done_message " Done." "";

echo "Checking bowtie ..."

BOWTIE_VER=`bowtie2 2>&1 | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'| perl -nle 'print $& if m{\d+\.\d+\.\d+}' `;

if ( hash bowtie2 2>/dev/null ) && ( echo $BOWTIE_VER | awk '{if($1>="2.1.0") exit 0; else exit 1}' )
then
   echo "bowtie2 >=2.1.0 found.";
else
   echo "bowtie >=2.1.0 or above not found. Trying to download from https://github.com/BenLangmead/bowtie2/archive/master.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/BenLangmead/bowtie2/archive/master.zip ext/opt/bowtie2.zip;
   unzip ext/opt/bowtie2.zip -d ext/opt/;
   cd ext/opt/bowtie2-master;
   make;
   cp bowtie2-build* $rootdir/ext/bin/.
   cp bowtie2-align* $rootdir/ext/bin/.
   cp bowtie2 $rootdir/ext/bin/.
   cd $rootdir;
fi;
done_message " Done." "";

echo "Checking SAMtools ..."

SAMTOOLS_VER=`samtools 2>&1 | perl -nle 'print $& if m{\d+\.\d+\.\d+}' | perl -nle 'print $& if m{\d+\.\d+\.\d+}' `;

if ( hash samtools 2>/dev/null ) && ( echo $SAMTOOLS_VER | awk '{if($1 =="0.1.20") exit 0; else exit 1}' )
then
   echo "samtools == 0.1.20 found.";
else
   echo "samtools == 0.1.20 not found. Trying to download from https://github.com/samtools/samtools/archive/0.1.20.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/samtools/samtools/archive/0.1.20.zip ext/opt/samtools-0.1.20.zip;
   unzip ext/opt/samtools-0.1.20.zip -d ext/opt/;
   cd ext/opt/samtools-0.1.20;
   make;
   cd $rootdir;
   cp ext/opt/samtools-0.1.20/samtools ext/bin/;
   cp ext/opt/samtools-0.1.20/bcftools/bcftools ext/bin/;
   cp ext/opt/samtools-0.1.20/bcftools/vcfutils.pl ext/bin/;
fi;
done_message " Done." "";

echo "Checking FastTree ..."

FASTTREE_VER=`FastTreeMP 2>&1 | perl -nle 'print $& if m{\d+\.\d+\.\d+}' | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`;

if ( hash FastTreeMP 2>/dev/null ) && ( echo $FASTTREE_VER | awk '{if($1>="2.1.8") exit 0; else exit 1}' )
then
   echo "FastTree >=2.1.8 found.";
else
   echo "FastTree >=2.1.8 not found. Trying to download from http://www.microbesonline.org/fasttree/FastTree.c ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext http://www.microbesonline.org/fasttree/FastTree.c ext/opt/FastTree.c;
   cd ext/opt/;
   if [[ "$OSTYPE" == "darwin"* ]]
	then
   		gcc -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm
   	else
   		gcc -DOPENMP -DUSE_DOUBLE -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm ;
   fi;
   cd $rootdir;
   cp ext/opt/FastTreeMP ext/bin/;
fi;
done_message " Done." "";

echo "Checking RAxML ..."

RAXML_VER=`raxmlHPC-PTHREADS -version 2>&1 | perl -nle 'print $& if m{\d+\.\d+\.\d+}' | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`;

if ( hash raxmlHPC-PTHREADS 2>/dev/null ) && ( echo $RAXML_VER | awk '{if($1>="8.0") exit 0; else exit 1}' )
then
   echo "RAxML >=8.0 found.";
else
   echo "RAxML >=8.0 not found. Trying to download from https://github.com/stamatak/standard-RAxML/archive/master.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/stamatak/standard-RAxML/archive/master.zip ext/opt/RAxML8.zip;
   unzip ext/opt/RAxML8.zip -d ext/opt/;
   cd ext/opt/standard-RAxML-master;
   make -f Makefile.SSE3.PTHREADS.gcc
   cd $rootdir;
   cp ext/opt/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 ext/bin/raxmlHPC-PTHREADS;
fi;
done_message " Done." "";

echo "Checking MUSCLE ..."

MUSCLE_VER=`muscle -version 2>&1 | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`;

if ( hash muscle 2>/dev/null ) && ( echo $MUSCLE_VER | awk '{if($1>="3.8.0") exit 0; else exit 1}' )
then
   echo "MUSCLE >=3.8.0 found.";
else
   echo "MUSCLE >=3.8.0 not found. Trying to download from http://www.drive5.com/muscle/downloads3.8.31/ ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   #determine platform
   MUSCULE_URL="http://www.drive5.com/muscle/downloads3.8.31"
   FILENAME="muscle3.8.31_i86linux64"
   UNAME=`uname`
   UNAMEM=`uname -m`
   if echo $UNAME | grep -i "linux" | grep -i "x86_64" > /dev/null; then
     MUSCULE_URL="$MUSCULE_URL/${FILENAME}.tar.gz"
   elif echo $UNAME | grep -i "^linux"  > /dev/null; then
     FILENAME="muscle3.8.31_i86linux32"
     MUSCULE_URL="$MUSCULE_URL/${FILENAME}.tar.gz"
   elif echo $UNAME | grep -i "darvin" |  grep -i "x86_64" > /dev/null; then
     FILENAME="muscle3.8.31_i86darwin64"
     MUSCULE_URL="$MUSCULE_URL/${FILENAME}.tar.gz"
   elif echo $UNAME | grep -i "darvin" > /dev/null; then
   FILENAME="muscle3.8.31_i86darwin32"
     MUSCULE_URL="$MUSCULE_URL/${FILENAME}.tar.gz"
   fi; 

   download_ext $MUSCULE_URL ext/opt/muscle.tar.gz;
   cd ext/opt/
   tar xvzf muscle.tar.gz
   mv $FILENAME $rootdir/ext/bin/
   cd $rootdir;
fi;
done_message " Done." "";

echo "Checking MAFFT ..."

MAFFT_VER=`mafft -version 2>&1 | grep "MAFFT v" | perl -nle 'print $& if m{\d+\.\d+}' `;

if ( hash mafft 2>/dev/null ) && ( echo $MAFFT_VER | awk '{if($1 >= 7.221) exit 0; else exit 1}' )
then
   echo "MAFFT >=7.221 found.";
else
   echo "MAFFT >=7.221 not found. Trying to download from http://mafft.cbrc.jp/alignment/software/mafft-7.221-without-extensions-src.tgz ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext http://mafft.cbrc.jp/alignment/software/mafft-7.221-without-extensions-src.tgz ext/opt/mafft-7.221-without-extensions-src.tgz;
   cd ext/opt/
   tar xvzf mafft-7.221-without-extensions-src.tgz
   cd mafft-7.221-without-extensions/core;
   sed -i.bak 's,PREFIX = \/usr\/local,'"PREFIX = $rootdir/ext"',' Makefile 
   make clean && make && make install
   cd $rootdir;
fi;
done_message " Done." "";


echo "Checking PAL2NAL ..."

PAL2NAL_VER=`pal2nal.pl 2>&1 | grep "pal2nal.pl" |  perl -nle 'print $& if m{v\d+}' | perl -nle 'print $& if m{\d+}' `;

if ( hash pal2nal.pl 2>/dev/null ) && ( echo $PAL2NAL_VER | awk '{if($1>="14") exit 0; else exit 1}' )
then
   echo "PAL2NAL >=14 found.";
else
   echo "PAL2NAL >=14 not found. Trying to download from http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz ext/opt/pal2nal.v14.tar.gz;
   cd ext/opt/
   tar xvzf pal2nal.v14.tar.gz
   cd $rootdir;
   cp ext/opt/pal2nal.v14/pal2nal.pl ext/bin/;
fi;
done_message " Done." "";

echo "Checking PAML ..."

PAML_VER=`evolver \0 2>&1 | grep "version" |  perl -nle 'print $& if m{\d+\.\d+}'`;

if ( hash evolver 2>/dev/null ) && ( echo $PAML_VER | awk '{if($1>=4.8) exit 0; else exit 1}' )
then
   echo "PAML >=4.8 found.";
else
   echo "PAML >=4.8 not found. Trying to download from http://abacus.gene.ucl.ac.uk/software/pamlX1.3.1+paml4.8a-win32.tgz ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext http://abacus.gene.ucl.ac.uk/software/pamlX1.3.1+paml4.8a-win32.tgz ext/opt/pamlX1.3.1+paml4.8a-win32.tgz;
   cd ext/opt/
   tar xvzf pamlX1.3.1+paml4.8a-win32.tgz
   cd paml4.8/;
   rm bin/*.exe
   cd src
   make -f Makefile
   cp baseml basemlg codeml evolver pamp yn00 mcmctree chi2 $rootdir/ext/bin/;
   cd $rootdir
fi;
done_message " Done." "";

echo "Checking CMake the cross-platform, open-source build system. ..."

CMake_VER=`cmake -version 2>&1 | perl -nle 'print $& if m{\d+\.\d+\.\d+}' `;

if ( hash cmake 2>/dev/null ) && ( echo $CMake_VER | awk '{if($_>="3.0.0") exit 0; else exit 1}' )
then
   echo "CMake >=3.0.0 found."
else
   echo "CMake >=3.0.0 not found. Trying to download from https://github.com/Kitware/CMake/archive/master.zip ..."
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/Kitware/CMake/archive/master.zip ext/opt/CMake.zip;
   unzip ext/opt/CMake.zip -d ext/opt/;
   cd ext/opt/CMake-master;
   ./bootstrap --prefix=$rootdir/ext && make && make install
   cd $rootdir
fi;
done_message " Done." "";

echo "Checking HyPhy ..."

HyPhy_VER=`echo -e "1\n2\n3" | HYPHYMP  2>&1 | grep HYPHY | perl -nle 'print $& if m{\d+\.\d+}'`;

if ( hash HYPHYMP 2>/dev/null ) && ( echo $HyPhy_VER | awk '{if($_>="2.2") exit 0; else exit 1}' )
then
   echo "HyPhy >=2.2 found.";
else
   echo "HyPhy >=2.2 not found. Trying to download from https://github.com/veg/hyphy/archive/master.zip ...";
   mkdir -p ext/opt;
   mkdir -p ext/bin;
   download_ext https://github.com/veg/hyphy/archive/master.zip ext/opt/HyPhy.zip;
   unzip ext/opt/HyPhy.zip -d ext/opt/;
   cd ext/opt/hyphy-master;
   cmake -DINSTALL_PREFIX=$rootdir/ext
   make MP2 && make install
   make GTEST
    ./HYPHYGTEST
   cd $rootdir
fi;
done_message " Done." "";

echo "Installing Perl dependencies..."

if hash cpanm 2>/dev/null; then
  echo "cpanm found. Start installing Perl dependencies..."
else
  echo "cpanm not found. Downloading from http://cpanmin.us..."
  download_ext http://cpanmin.us ext/bin/cpanm;
  chmod a+x ext/bin/cpanm;
fi
( set -xe;
  perl -MGetopt::Long -e 1 > /dev/null 2>&1       || cpanm -v --notest -l ext Getopt::Long;
  perl -MTime::HiRes -e 1 > /dev/null 2>&1        || cpanm -v --notest -l ext Time::HiRes;
  perl -MFile::Path -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext File::Path;
  perl -MFile::Basename -e 1 > /dev/null 2>&1     || cpanm -v --notest -l ext File::Basename;
  perl -MFile::Copy -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext File::Copy;
  perl -MIO::Handle -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext IO::Handle;
  perl -MParallel::ForkManager -e 1 > /dev/null 2>&1           || cpanm -v --notest -l ext Parallel::ForkManager;
  perl -MStatistics::Distributions -e 1 > /dev/null 2>&1         || cpanm -v --notest -l ext Statistics::Distributions;
)
done_message "Done installing Perl dependencies." "Failed installing Perl dependencies.";

echo -n "Moving scripts...";
( set -e;
  cp src/*.pl bin/
  chmod a+x bin/*
)
done_message " Done." "Failed installing PhaME scripts.";

rm -f ext/opt/*zip ext/opt/*gz
rm -f evolver.out messages.log

echo "
================================================================================
                 PhaME installed successfully.
================================================================================
Check phame.ctl for the control file

Quick start:
    bin/runPhaME.pl phame.ctl
Check our github site for update:
    https://github.com/LANL-Bioinformatics/PhaME
";
