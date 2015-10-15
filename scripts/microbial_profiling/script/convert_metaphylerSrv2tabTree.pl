#!/usr/bin/perl -w
use strict;

my %taxa;

while (<STDIN>) {
	chomp;
	next if /@/;
	#"prefix".classify.tab
	#column 0: query sequence id
	#column 1: phylogenetic marker gene name
    #    column 2: best reference gene hit
    #    column 3: % similarity with best hit
    #    column 4: classification rule
	#column 5-9: taxonomic label at genus,family,order,class,phylum level

    #"prefix".genus|family|order|class|phylum.tab
    #    column 1: taxonomic clade name
    #    column 2: % relative abundances
    #    column 3: depth of coverage of genomes
    #    column 4: number of sequences binned to this clade
    #    column 5: similarity with reference genes (only available at the genus level)	

	my @temp = split /\t/, $_;
	my @taxas = grep {$_ ne "cutoff"} reverse @temp[5..9];

    my $lineage = join "\t", @taxas;
    $lineage =~ s/{\w+}//g;
    $lineage =~ s/([\(\),:])//g;

    $taxa{$lineage} ||= 0;
    $taxa{$lineage} ++;
}

foreach my $lineage ( sort {$taxa{$b}<=>$taxa{$a}} keys %taxa ){
    print "$taxa{$lineage}\t$lineage\n";
}
