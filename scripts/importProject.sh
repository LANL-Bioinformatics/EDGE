#!/usr/bin/sh
#projectNameList=$(ls */config.txt | sed 's$\(.*\)/config.txt$\1$g')
projectNameList=$(ls *.tar.gz | sed 's$\(.*\).tar.gz$\1$g')
echo "This program will import projects into EDGE with a user mangement system.  The projects need to be within your current directory and be in '*.tar.gz' format."
echo 
read -p "Please enter your EDGE login email address: " email

read -s -p "Please enter your EDGE Password: " password
echo


username=$(perl -swe 'use Digest::MD5 qw(md5_hex);print md5_hex($un);' -- -un=$email)

for projectName in $projectNameList
do
	echo "Unzipping $projectName"
	tar -xf $projectName".tar.gz"
	echo "Importing $projectName"
	idOut=$(perl edgeWSClientProjectAdd.pl $email $password $projectName)
	id=$(echo $idOut | sed 's${"id":\(.*\)}$\1$g')
	echo "Project ID: "$id", Project Name: "$projectName
	cp -r $projectName $EDGE_HOME"/edge_ui/EDGE_output/"$id
	ln -s $EDGE_HOME"/edge_ui/EDGE_output/"$id $EDGE_HOME"/edge_ui/EDGE_input/"$username"/MyProjects/"$projectName
	ln -s $EDGE_HOME"/edge_ui/EDGE_output/"$id $EDGE_HOME"/edge_ui/JBrowse/data/"$id
	rm -r $projectName/
	echo "Done importing "$porjectName". Please check EDGE to confirm it worked."
done

echo "Done importing projects"
