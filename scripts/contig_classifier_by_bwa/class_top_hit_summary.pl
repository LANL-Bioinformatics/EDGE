#!/usr/bin/perl
use strict;

my $proj = $ARGV[0];
my $rank = $ARGV[1];
my $contig;

while(<STDIN>){
    my @temp = split /\t/, $_;
    next if defined $rank && $temp[1] ne $rank;
    next if defined $contig->{$temp[0]}->{$temp[1]};
    $contig->{$temp[0]}->{$temp[1]} = 1;
    my $out = $_;
    $out =~ s/^[^#]\S+\t/$proj\t/ if defined $proj && $proj;
    print $out;
}
