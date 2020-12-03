#!/usr/bin/env perl

my $col = 5; #genus
$col = 6 if $ARGV[0] eq "species";
$col = 7 if $ARGV[0] eq "strain";

my $r;

while(<STDIN>){
	chomp;
	my @temp = split /\t/;
	my $lineage = join "\t", @temp[1..$col];
	$r->{$lineage} ||= 0;
	$r->{$lineage} += $temp[0];
}

foreach my $l ( sort {$a cmp $b} keys %$r){
	print "$r->{$l}\t$l\n";
}
