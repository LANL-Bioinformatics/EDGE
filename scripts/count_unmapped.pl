#!/usr/bin/env perl

use strict;
use warnings;
my @bam=@ARGV;
my $ref_count=scalar(@bam);
my $unmapped=0;

my %filter;
foreach my $bam (@bam){
	open (my $fh, "samtools view -f 4 $bam | ") or die "Cannot read $bam";

	while(<$fh>){
		chomp;
        	my @samFields=split /\t/,$_;
        	my $R1_R2 = 1;
        	$R1_R2 = 2 if ($samFields[1] & 128);
        	my $unique_id=$samFields[0]."_$R1_R2";
		#print $unique_id,"\n";
        	$filter{$unique_id}++;
        	next if ($filter{$unique_id} != $ref_count);
		$unmapped++;
	}
	close $fh;
}
print $unmapped,"\n";
