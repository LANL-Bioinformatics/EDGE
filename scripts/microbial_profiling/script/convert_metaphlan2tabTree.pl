#!/usr/bin/env perl
use strict;

my $tax;


while(<STDIN>){
	chomp;
	next if /LEVEL/;
	next unless /\|s__/;
	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;

    $phylo =~ s/^\w__/\t/g;
    $phylo =~ s/\|\w__/\t/g;
    $phylo =~ s/_/ /g;

    $tax->{$phylo} = $count;

}

foreach my $phylo ( sort {$a cmp $b} keys %$tax ){
    print "$tax->{$phylo}\t$phylo\n";
}
