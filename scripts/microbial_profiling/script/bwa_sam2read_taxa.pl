#! /usr/bin/env perl
#
# This script is used to output read mapping results in 2 columns: READ_NAME and TAXONOMY.
# 
# USAGE:
#   
#   bwa_sam2read_taxa.pl ([RANK]) ("preload") < [sam file] > [output]
#  
# The default RANK is species. That can be changed to any rank. The "preload" option will load the
# entire GI-taxonomy mapping table into memory. That can increase the speed if your input is a
# very large SAM file (>10G). Note that "preload" option can only work properly with perl 5.16 or 
# above.
#
# EXAMPLE:
#   
#   samtools view sample1.bam | bwa_sam2classification.pl > output.classification
#   cat sample2.sam | bwa_sam2classification.pl genus preload > output.genus_classification
# 
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
# 2014/02/21
#

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage;

my $rank=$ARGV[0];
$rank ||= "species";

my $preload=$ARGV[1];
loadTaxonomy($preload);

while(<STDIN>)
{
	next if /^\@SQ/;
	chomp;
	my @fields = split /\t/, $_;

	#unmapped
	if( $fields[1] & 4 ){
		print "$fields[0]\t\n";
	}
	else{#mapped
		#my ($gi) = $fields[2] =~ /gi\|(\d+)/;
		my $acc = getAccFromSeqID($fields[2]);
		my $name = acc2rank($acc,$rank);
		if($name){
			print "$fields[0]\t$name\n";
		}
		else{
			print "$fields[0]\tno $rank ($fields[2])\n";
		}
	}
}
