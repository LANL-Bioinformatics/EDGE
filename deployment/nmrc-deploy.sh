#!/bin/bash
####################################
## This is the EDGE Bioinformatics Deployment Configurator Script.
## It assumes a bare install of CentOS 7 with Development Tools and GNOME desk-
## top, rebooted through first boot.  This is the deployment script used for DoD
## installations.
##
## A network connection with port 80 open is required for installation.  You
## must be able to route to centos.org, cpan.org, and fedoraproject.org as a
## minimum.  On a low-speed connection, this script may take several hours.
##
## WARNING: If you run this script on a server that's already been configured to
## do something else, it will probably break.  If you're running this script
## and you haven't met Joe from NMRC and Patrick from LANL, it's probably not
## time to run this script yet.
##
## This is the work of an employee of the United States Government, performed as
## part of their official duties.  No copyright is claimed, this work is entered
## into the public domain without restriction.
##
####################################

## Update changePassword to the appropriate password for the system
## Update changePassword in createDatabase.sql to the appropriate password

## Install pre-reqs
sudo yum install -y epel-release
sudo yum install -y libX11-devel readline-devel libXt-devel ncurses-devel inkscape scipy expat expat-devel freetype freetype-devel zlib zlib-devel perl-App-cpanminus perl-Test-Most python-pip blas-devel atlas-devel lapack-devel numpy numpy-f2py libpng12 libpng12-devel perl-XML-Simple perl-JSON csh gcc gcc-c++ make binutils gd gsl-devel git graphviz java-1.7.0-openjdk perl-Archive-Zip perl-CGI perl-CGI-Session perl-CPAN-Meta-YAML perl-DBI perl-Data-Dumper perl-GD perl-IO-Compress perl-Module-Build perl-XML-LibXML perl-XML-Parser perl-XML-SAX perl-XML-SAX-Writer perl-XML-Twig perl-XML-Writer perl-YAML perl-PerlIO-gzip python-matplotlib python-six libstdc++-static 



## Update existing python and perl tools
sudo pip install --upgrade six scipy matplotlib
sudo cpanm App::cpanoutdated
sudo su -
cpan-outdated -p | cpanm
exit

## Install more Perl modules
## Some of these may fail, that's okay.  BioPerl needs to go smoothly.
cpanm Graph Time::Piece BioPerl
cpanm Algorithm::Munkres Archive::Tar Array::Compare Clone Convert::Binary::C
cpanm HTML::Template HTML::TableExtract List::MoreUtils PostScript::TextBlock
cpanm SOAP::Lite SVG SVG::Graph Set::Scalar Sort::Naturally Spreadsheet::ParseExcel
cpanm CGI CGI::Simple GD Graph GraphViz XML::Parser::PerlSAX XML::SAX XML::SAX::Writer XML::Simple XML::Twig XML::Writer 

## Install Apache for the Web GUI
sudo yum install -y httpd httpd-tools

## Configure firewall for ssh, http, https, and smtp:
sudo firewall-cmd --permanent --add-service=ssh
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --permanent --add-service=smtp

## Try to clone from the EDGE repo if SSH works
cd ~
if [ -d /home/edge/edge ]; then 
  echo "EDGE directory exists at ~/edge, skipping."
elif [ $1 == "SSH" ]; then
  echo "SSH deployment from BitBucket, requires registered SSH key."
  git clone git@bitbucket.org:nmrcjoe/edge.git
elif [ $1 == "HTTPS" ]; then
  echo "HTTPS deployment from BitBucket, expect timeouts."
  git clone https://nmrcjoe@bitbucket.org/nmrcjoe/edge.git
elif [ $1 == "NMRC" ]; then
  echo "Will use rsync to copy from bigsilver."
  rsync -avzhr --progress joe@192.168.20.3:~/edge ~/
else 
  echo "Program directory not present, and no valid source specified.  Exiting."
  exit 0
fi

## Do the same for the database files (if not done already)
if ![ -d /home/edge/database ]; then
  echo "Retrieving database files from bigsilver."
  rsync -avzhr --progress joe@192.168.20.3:~/database ~/
fi


## Disable SELINUX
# sudo sed -i 's/SELINUX=enforcing/SELINUX=disabled/g' /etc/sysconfig/selinux
# sudo sed -i 's/SELINUX=enforcing/SELINUX=disabled/g' /etc/selinux/config
sudo setenforce 0

## Set default NMRC EDGE params
sed -i 's/opt\/apps/home\/edge/g' /home/edge/edge/edge_ui/cgi-bin/edge_config.tmpl
sed -i 's/user_management=1/user_management=0/g' /home/edge/edge/edge_ui/cgi-bin/edge_config.tmpl


## Fix the database directories (Assumes dbs are in ~/database)
## Make sure the database is linked prior to running INSTALL.sh
rm -rf ~/edge/database
ln -s ~/database ~/edge/database
sudo ln -s ~/database/ /database
## This should be done by the INSTALL.sh
#ln -s ~/database/Krona_taxonomy ~/edge/thirdParty/KronaTools-2.4/taxonomy

## Install LANL EDGE
~/edge/INSTALL.sh

## Copy the EDGE httpd conf files to the appropriate directories
## This should be done after INSTALL.sh, the script inserts the appropriate paths
sudo cp ~/edge/edge_ui/apache_conf/edge_httpd.conf /etc/httpd/conf.d/
sudo cp ~/edge/deployment/httpd.conf /etc/httpd/conf/

## Setup userManagement
## Install database
sudo yum install mariadb-server mariadb
sudo systemctl start mariadb.service && sudo systemctl enable mariadb.service
## Setup root password on database
sudo mysql_secure_installation

## Create userManagement database and Load schema/constrains
## Update changePassword within createDatabase.sql to the appropriate password
echo createDatabase.sql | mysql -u root -p

## Install php
sudo yum install php php-pear
sudo yum install php-mysql
sudo httpd -k restart

## Install and Configure tomcat
sudo yum install tomcat
sudo cp ~/userManagement/mariadb-java-client-1.2.0.jar /usr/share/tomcat/lib
sudo sed -i 's@<!-- <role rolename="admin"/> -->@<!-- <role rolename="admin"/> -->\n<role rolename="admin"/>\n<user username="edge" password="changePassword" roles="admin"/>@g' /usr/share/tomcat/conf/tomcat-users.xml
sudo sed -i 's@<session-timeout>.*</session-timeout>@<session-timeout>4320</session-timeout>@g' /usr/share/tomcat/conf/web.xml
sudo sed -i 's@#JAVA_OPTS@JAVA_OPTS="-Xms256m -Xmx1024m -XX:PermSize=256m -XX:MaxPermSize=512m"\n#JAVA_OPTS@g' /usr/share/tomcat/conf/tomcat.conf

## Deploy userManagement to tomcat server
sudo cp ~/userManagement/userManagement*.war /usr/share/tomcat/webapps/

## Edit ~/edge/userManagement/userManagementWS.xml then deploy it to /usr/share/tomcat/conf/Catalina/localhost
sed -i 's@username=.*$@username="edge"@' ~/edge/userManagement/userManagementWS.xml
sed -i 's@password=.*$@password="changePasword"@' ~/edge/userManagement/userManagementWS.xml
sed -i 's@driverClassName=.*$@driverClassName="org.mariadb.jdbc.Driver"@' ~/edge/userManagement/userManagementWS.xml
sudo cp ~/userManagement/userManagementWS.xml /usr/share/tomcat/conf/Catalina/localhost/

## Edit /usr/share/tomcat/webapps/userManagement/WEB-INF/classes/sys.properties to match the appropriate settings for the server
sudo sed -i 's@host_url=.*$@host_url=http://localhost:8080/userManagement@g'  /usr/share/tomcat/webapps/userManagement/WEB-INF/classes/sys.properties
sudo sed -i 's@wsURL=.*$@wsURL=http://localhost:8080/userManagementWS@g'  /usr/share/tomcat/webapps/userManagement/WEB-INF/classes/sys.properties
sudo sed -i 's@email_notification=.*$@email_notification=off@g'  /usr/share/tomcat/webapps/userManagement/WEB-INF/classes/sys.properties

## Create admin account for EDGE
## Should this script fail, something is not set up correctly
perl ~/edge/userManagement/createAdminAccount.pl -e admin@edge.com -p changePassword -fn admin -ln edge

## Enable userManagement in ~/edge/edge_ui/sys.properties
sed -i 's@user_management=.*$@user_management=1@g' ~/edge/edge_ui/sys.properties
sed -i 's@edge_user_management_url=.*$@edge_user_management_url=http://localhost:8080/userManagement@g' ~/edge/edge_ui/sys.properties

