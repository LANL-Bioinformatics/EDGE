#!/usr/bin/perl

######################################################
# Written by Sanaa Ahmed
# Jan. 18, 2013

# Given a fasta files
#  Removes all repeats from the sequence
#  Generates a multi-fasta files for each sequence.

######################################################

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use File::Basename;

my $identity=95;
my $len_cutoff=0;
my $outfasta= "non_repeat.fna";
my $repeat_stats="non_repeat_stats.txt";
my ($file,$coords);

#my %regions;

GetOptions(
   'f=s'      => \$file,
   'c=s'      => \$coords,
   'o=s'      => \$outfasta,
   's=s'      => \$repeat_stats,
   'help|?'   => sub{Usage()},
);

sub usage
{
print "
$0 [options] 
   
   -f STRING  reference sequence filename
   -c STRING  tab-delimited coords filename
   -o STRING  output multifasta filename (default: non_repeats.fna)
   -s STRING  output repeats stats filename (default: non_repeat_stats.txt)
";
exit;
}

&usage unless ($file);


my ($header,@seq,$sequence);
my $length=0;
my $fh= FileHandle->new($file)||die "$!";
if ($fh->open("< $file")){
   $/=">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless($_){next;};
      ($header,@seq)=split /\n/,$_;
      $sequence= join "",@seq;
      $length = length $sequence;
   }
   $/="\n";
   $fh->close;
}

if (-z $coords){
   open (OUT, ">$outfasta");
   print OUT ">${header}_1_$length\n$sequence\n";
   close OUT;
   exit;
}

my %repeats;
open (IN, $coords) || die "$!";
while (<IN>){
   chomp;
   my ($ref,$rep,$rstart,$rend,$length)= split /\s/,$_;
   $repeats{$rstart}=$rend;
}
&get_sequence_file($header);

sub get_sequence_file
{
my $reference=shift;
my $first=1;
my ($start,$end);
my $contig;
my ($lastend,$laststart)=(0,0);
my ($r_start,$r_end)=(0,0);
my $count=0;
$repeat_stats= 'non_repeats_'.$reference.'_stats.txt';

open (OUT, ">$outfasta");
open (STAT,">$repeat_stats");
print STAT "Reference\tLength\tNon-Repeat\tRepeat size\n$reference\t$length\t";
foreach my $repeat(sort {$a<=>$b} keys %repeats){
#   print "$repeat\n";
   $r_start=$repeat;
   $r_end= $repeats{$repeat};
   if ($r_start>$lastend && $r_end >$lastend){
#      print "$gap\t$gaps{$gap}\n";
      $lastend=$r_end;
      if ($first){
         if ($r_start==1){$start=$r_end;}
         else{
            $start=1;
            $end=$r_start+1;
            $contig=$reference.'_'.$start.'_'.$end;
#            print "$start\t$end\n";
            my $output= substr($sequence,$start-1,$end-$start+1);
            print OUT ">$contig\n$output\n";
            $count= $count+ length $output;
         }
         $first=0;
      }
      else{
         $end=$r_start+1;
         $contig=$reference.'_'.$start.'_'.$end;
#         print "$start\t$end\n";
         my $output= substr($sequence,$start-1,$end-$start+1);
         print OUT ">$contig\n$output\n";
         $count= $count+ length $output;
      }
      $start=$r_end;
   }
}
if ($length>$r_end){
   $start= $r_end;
   $end= $length;
   $contig=$reference.'_'.$start.'_'.$end;
#   print "$start\t$end\n";
   my $output= substr($sequence,$start-1,$end-$start+1);
   print OUT ">$contig\n$output\n";
   $count= $count+ length $output;
}
my $difference= $length-$count;
print STAT "$count\t$difference\n";
}
close OUT;
close IN;


