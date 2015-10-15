#!/bin/bash

for i in $@; do
  inkscape --without-gui --export-pdf="$(dirname $i)/$(basename $i .svg).pdf" $i
done
