#!/bin/bash

CAIROSVG_VER=`cairosvg --version 2>&1 | perl -nle 'print $& if m{\d+\.\d+}'`;

if ( hash inkscape 2>/dev/null ) 
then
	for i in $@; do
  		inkscape --without-gui --export-pdf="$(dirname $i)/$(basename $i .svg).pdf" $i
	done
elif ( hash cairosvg 2>/dev/null ) && ( echo $CAIROSVG_VER | awk '{if($1>="2.0") exit 0; else exit 1}' )
then
	for i in $@; do
                cairosvg --output="$(dirname $i)/$(basename $i .svg).pdf" $i
        done
else
	echo "No svg to pdf software" 
	exit 1;
fi

