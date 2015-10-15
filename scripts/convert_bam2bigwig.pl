#!/usr/bin/env perl
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
# 2014/02/21

use strict;

my @files = @ARGV;
my @bw;

foreach my $file ( @files )
{
	my $tempfile = "$file.genomeSize";
	open GSZ, ">$tempfile" or die $!;
	open REF, "samtools view -H $file |" or die $!;
	while(<REF>){
	    if( /SN:(\S+)\tLN:(\d+)/ ){
	        print GSZ "$1\t$2\n";
	    }
	}
	close REF;
	close GSZ;

	`genomeCoverageBed -split -bg -ibam $file -g $tempfile > $file.bedgraph`;
	`wigToBigWig $file.bedgraph $tempfile $file.bw`;
}
