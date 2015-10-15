#!/usr/bin/perl
use strict;

while(<STDIN>){
	chomp;
	next if /LEVEL/;
	next unless /\|s__/;
	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;

    $phylo =~ s/^\w__/\t/g;
    $phylo =~ s/\|\w__/\t/g;
    $phylo =~ s/_/ /g;

    print "$count\t$phylo\n"
}
