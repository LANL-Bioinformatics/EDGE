#!/usr/bin/env perl

use strict;
use File::Basename;

if (scalar(@ARGV)<2)
{
   print "usage: $0 <Coverage_file> <fold>\n";
   exit
}
my ($prefix, $tmp1, $tmp2) = fileparse("$ARGV[0]", qr/\.[^.]*/);
my $cut_off = $ARGV[1];
$cut_off=0 if (!defined $cut_off);
my ($ref_id) = ($prefix =~ /readsToRef_(\S+)/);
print "Start\tEnd\tLength\tRef_ID\n";
my @gap_array;
my $genome_bases;
my $covered_base_num;
my $gap_length;
my @cov_array;
my ($total_length, $count)=(0,0);
open (IN,"$ARGV[0]");
while(<IN>)
{
  next if ($_ !~ /^\d/);
  my ($pos, $cov) = split;
  if ($cov > $cut_off){
     if (@gap_array){
      $gap_length = $gap_array[-1] - $gap_array[0]+1;
      print $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t$ref_id\n";
      $count++;
      $total_length += $gap_length;
      @gap_array=();
     }
     $covered_base_num++;
  }
  elsif ($cov<=$cut_off)
  {
     push @gap_array, $pos;
  }
  $genome_bases++;
  push @cov_array,$cov;

}

 if (@gap_array){
      $gap_length = $gap_array[-1] - $gap_array[0]+1;
      print $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t$ref_id\n";
      $count++;
      $total_length += $gap_length;
      @gap_array=();
 }
close IN;

#my ($std_cov,$avg_cov)= &standard_deviation(@cov_array);
#printf STDERR ("Genome_coverage:\t%.2f%%\n",$covered_base_num/$genome_bases*100);
#printf STDERR ("Mean_fold(x):\t%.2f\t%s %.2f\n",$avg_cov,"Std",$std_cov);
#print STDERR "$count ($total_length)\n";
#print STDERR "Gap:\t$count\n";
#print STDERR "Total_Gap_bases:\t$total_length\n";

sub standard_deviation {
my(@numbers) = @_;
#Prevent division by 0 error in case you get junk data
return undef unless(scalar(@numbers));

# Step 1, find the mean of the numbers
my $total1 = 0;
foreach my $num (@numbers) {
$total1 += $num;
}
my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
my $total2 = 0;
foreach my $num (@numbers) {
$total2 += ($mean1-$num)**2;
}
my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
my $std_dev = sqrt($mean2);
return ($std_dev,$mean1);
}

