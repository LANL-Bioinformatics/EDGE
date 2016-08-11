#!/usr/bin/env perl
#
use strict;
use File::Basename;
use Getopt::Long;

my @input;
my $usage = qq{ 
  Usage: perl $0 -i <gff gff ...> 
};

GetOptions('i=s{,}'=> \@input,
          )||die $usage;

die $usage if (!@input);

my $seq_end=0;
my $new_start=0;
foreach my $input_index (0..$#input){
	$new_start = $seq_end + $new_start;
	my ($tmp,$seq_id,$seq_start);
        open (my $fh, $input[$input_index]) or die "$!\n";
	while(<$fh>){
		my @gff_features = split(/\t/,$_);
		if(/sequence-region/){
			($tmp,$seq_id,$seq_start,$seq_end)=split(/\s+/,$_);
		}
		if ($input_index > 0){
			if ( scalar(@gff_features) == 9 ){
				$gff_features[3] = $gff_features[3] + $new_start;
				$gff_features[4] = $gff_features[4] + $new_start;
				print join ("\t",@gff_features);
			}else{
				print $_;
			}
		}else{
			print $_;
		}
	}
	close $fh;
}
