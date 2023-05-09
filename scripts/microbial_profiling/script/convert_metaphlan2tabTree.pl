#!/usr/bin/env perl
use strict;

my $tax;


while(<STDIN>){
	chomp;
	next if /LEVEL/;
	next unless /\|s__/;
	next if (/^#/);
        my @result = split /\t/, $_;
        my ($phylo,$count);
        if ( scalar(@result) > 2 ){
                $phylo=$result[0];
                $count=$result[2];

        }else{
                $phylo=$result[0];
                $count=$result[1];
        }

#	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;

    $phylo =~ s/^\w__/\t/g;
    $phylo =~ s/\|\w__/\t/g;
    $phylo =~ s/_/ /g;

    $tax->{$phylo} = $count;

}

foreach my $phylo ( sort {$a cmp $b} keys %$tax ){
    print "$tax->{$phylo}\t$phylo\n";
}
