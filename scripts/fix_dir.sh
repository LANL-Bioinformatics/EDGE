#!/usr/bin/env bash

outdir=$(grep edgeui_output  $EDGE_HOME/edge_ui/cgi-bin/edge_config.tmpl | awk -F"=" '{print $2}')

cd $outdir

for i in */SNP_Phylogeny/*/SNPtree.finished; do name=${i%/*}; mv $name/* $name/../ ;done

for i in *
do
if [ -f "$i/ReadsBasedAnalysis/runReadsToGenome.finished" ]
then
    mkdir -p $outdir/$i/ReadsBasedAnalysis/readsMappingToRef
    mv $i/ReadsBasedAnalysis/runReadsToGenome.finished $i/ReadsBasedAnalysis/readsMappingToRef/
    mv $i/ReadsBasedAnalysis/readsToRef* $i/ReadsBasedAnalysis/readsMappingToRef/
fi

if [[ -d "$i/Reference" ]] && [[ ! -d "$i/ReferenceBasedAnalysis" ]]
then
    echo "Fix project $i dir ..."
    mkdir -p $i/ReferenceBasedAnalysis
    if [ -d "$i/AssemblyBasedAnalysis/Blast" ]
    then 
        ln -sf $outdir/$i/AssemblyBasedAnalysis/Blast $i/ReferenceBasedAnalysis/UnmapppedContigs
    fi
    ln -sf $outdir/$i/AssemblyBasedAnalysis/contigMappingToRef $i/ReferenceBasedAnalysis/contigMappingToRef
  #  mkdir -p $outdir/$i/ReadsBasedAnalysis/readsMappingToRef
  #  if [ -f "$i/ReadsBasedAnalysis/readsToRef.sort.bam" ]
  #  then
  #      mv $i/ReadsBasedAnalysis/* $i/ReadsBasedAnalysis/readsMappingToRef/
  #      if [ -d "$i/ReadsBasedAnalysis/readsMappingToRef/Taxonomy" ]
  #      then
  #          mv $i/ReadsBasedAnalysis/readsMappingToRef/Taxonomy $i/ReadsBasedAnalysis/
  #      fi
  #      if [ -d "$i/ReadsBasedAnalysis/readsMappingToRef/UnmappedReads" ]
  #      then
  #          mv $i/ReadsBasedAnalysis/readsMappingToRef/UnmappedReads $i/ReadsBasedAnalysis/UnmappedReads
  #          ln -sf $outdir/$i/ReadsBasedAnalysis/UnmappedReads $i/ReferenceBasedAnalysis/
  #      fi
  #  fi
    ln -sf $outdir/$i/ReadsBasedAnalysis/readsMappingToRef $i/ReferenceBasedAnalysis/readsMappingToRef
fi
done
