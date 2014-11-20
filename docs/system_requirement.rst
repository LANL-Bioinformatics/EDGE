.. _sys_requirement:

System requirements
###################

The current version of EDGE pipeline has been extensively tested on a Linux Server with Ubuntu 14.04 and Centos 6.5 and 7.0 operating system and will work on 64bit Linux environments. Perl v5.8 or above is required. Python 2.7 is required. Due to the involvement of several memory/time consuming steps, it requires at least 16Gb memory and at least 8 computing CPUs. A higher computer spec is recommended: 128Gb memory and 16 computing CPUs.

Please assure that your system has the essential software building packages installed properly before running the installing script.

The following are required installed by system administrator.

Ubuntu 14.04 
============

.. image:: https://design.ubuntu.com/wp-content/uploads/ubuntu-logo14.png
    :width: 200px

1. Install build essential, libraries and dependancies::
    
    sudo apt-get install build-essential
    sudo apt-get install libreadline5-dev
    sudo apt-get install libx11-dev
    sudo apt-get install libxt-dev
    sudo apt-get install libncurses5-dev 
    sudo apt-get install gfortran
    sudo apt-get install inkscape
    sudo apt-get install libwww-perl
    sudo apt-get install zlib1g-dev zip unzip
    sudo apt-get install libpng-dev
    sudo apt-get install cpanminus
    sudo apt-get install firefox

2. Install python packages for Metaphlan (Taxonomy assignment software)::
   
    sudo apt-get install python-numpy python-matplotlib python-scipy libpython2.7-stdlib 
    sudo apt-get install ipython ipython-notebook python-pandas python-sympy python-nose
  
3. Install BioPerl::
   
    sudo apt-get install bioperl  
        or
    sudo cpan -i -f CJFIELDS/BioPerl-1.6.923.tar.gz

CentOS 6
========

.. image:: https://scottlinux.com/wp-content/uploads/2011/07/centos6.png
    :width: 200px
    
1. Install dependancies using yum::

    sudo yum -y install perl-CPAN expat-devel gcc gcc-c++ kernel-devel inkscape perl-JSON 
    sudo yum -y install libXt-devel.x86_64 gcc-gfortran.x86_64 readline-devel.x86_64
    sudo yum -y install git zip zlib-devel perl-ExtUtils-Embed perl-CGI firefox curl wget 
    sudo yum -y install perl-App-cpanminus libX11-devel.x86_64 unzip 
    sudo yum -y install perl-libwww-perl fuse-exfat exfat–utils

2. Install Bioperl using CPAN::

    sudo cpan -i -f CJFIELDS/BioPerl-1.6.923.tar.gz

3. Install Anaconda python distribution (`https://store.continuum.io/cshop/anaconda/ <https://store.continuum.io/cshop/anaconda/>`_)
   The installation is interactive. Type in /opt/apps/anaconda when the script asks for the location to install python.::
             
    bash Anaconda-2.0.1-Linux-x86.sh
    ln -s /opt/apps/anaconda/bin/python /opt/apps/edge/bin/
    
  You have to symlink anaconda python to edge/bin. So system will use your python over the system’s. 


CentOS 7
========

.. image:: https://scottlinux.com/wp-content/uploads/2010/07/centos.png
    :width: 200px

EPEL Repo

1. Install libraries and dependancies by yum::

    sudo yum -y install libX11-devel readline-devel libXt-devel ncurses-devel inkscape 
    sudo yum -y install freetype freetype-devel zlib zlib-devel perl-App-cpanminus git 
    sudo yum -y install blas-devel atlas-devel lapack-devel  libpng12 libpng12-devel
    sudo yum -y install perl-XML-Simple perl-JSON expat expat-devel perl-Test-Most
    sudo yum -y install python-pip numpy numpy-f2py scipy  
 
2. Upgrade python packages by pip::
   
    pip install --upgrade six
    pip install --upgrade scipy
    pip install --upgrade matplotlib

3. Detect outdated CPAN modules::

    cpanm App::cpanoutdated
    cpan-outdated -p | cpanm

4. Install perl modules by cpanm::
    
    cpanm Algorithm::Munkres Archive::Tar Array::Compare Clone CGI Convert::Binary::C GD 
    cpanm GraphViz HTML::Template HTML::TableExtract List::MoreUtils PostScript::TextBlock 
    cpanm SOAP::Lite SVG SVG::Graph Set::Scalar Sort::Naturally Spreadsheet::ParseExcel 
    cpanm XML::Parser::PerlSAX XML::SAX XML::SAX::Writer XML::Simple XML::Twig XML::Writer 
    cpanm Graph Time::Piece Tk BioPerl

