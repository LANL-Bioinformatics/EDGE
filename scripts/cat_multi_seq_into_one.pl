#! /usr/bin/perl

use File::Basename;
use strict;
use Getopt::Long;

my @input;
my $seq_name;
my $n_num=0;
my $usage = qq{
  Usage: perl $0 -i <fasta fasta ...> -name <seq_id>;
         -name  The output sequence name. (Default: using the prefix of file name)
         -n     Insert <INT> number of "N" between sequences. (default:0)   
};
GetOptions('i=s{,}'=> \@input,
	   'name=s'=> \$seq_name,
           'n=i'   => \$n_num,
          )||die $usage;

die $usage if (!@input);
my ($file_name, $file_path, $file_suffix)=fileparse("$input[0]", qr/\.[^.]*/);
$seq_name=$file_name if (!$seq_name);
print ">",$seq_name,"\n";

my $seq="";
foreach my $input (@input){
	open (IN, $input) or die "$!\n";
	while(<IN>)
	{
   		chomp;
   		if ($_=~/>/)
   		{ 
     		# ($species) = $_ =~ /^>\d+\s\w+\s(.*)/;
     			if ($seq)
     			{
         			$seq = $seq.("N" x $n_num) ; 
     			}
   		}
   		else
   		{
     			if ($_=~/\w/){
       				$seq .=$_;
     			} 
   		}
	}
	if ($seq)
	{
		 $seq = $seq.("N" x $n_num) ; 
	}
	close IN;
}

$seq =~ s/(.{80})/$1\n/g;
$seq =~ s/\n$//g; 
print $seq,"\n";

#print "$id\t$species\n";



