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
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage;

my $list;
my $score_cutoff = $ARGV[0];

while(<STDIN>)
{
	next if /^\@/;
	my ($flag, $seqid) = $_ =~ /^\S+\t(\d+)\t(\S+)/;
	my ($as)           = $_ =~ /AS:i:(\d+)/;
	my ($xs)           = $_ =~ /XS:i:(\d+)/;
	next unless defined $flag;
	unless( $flag & 4 )
	{
		next if $as < $score_cutoff && $score_cutoff > 0;
		#my ($gi) = $seqid =~ /gi\|(\d+)/;
		my $acc = getAccFromSeqID($seqid);
		my $gi = $acc;
		$list->{$gi}->{MAPPED} ||= 0;
		$list->{$gi}->{MAPPED}++;
		if( !$xs || $as>$xs ){
			$list->{$gi}->{UNIQUE} ||= 0;
			$list->{$gi}->{UNIQUE}++;
		}
	}
}

foreach my $gi ( keys %$list ){
	$list->{$gi}->{UNIQUE}||=0;
	print "$gi\t$list->{$gi}->{MAPPED}\t$list->{$gi}->{UNIQUE}\n";
}
