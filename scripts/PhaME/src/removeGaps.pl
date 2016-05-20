#!/usr/bin/perl 

######################################################
# Written by Sanaa Ahmed
# Jan. 03, 2013

# Given a fasta file and a gap file
#  pulls non-gap sequence regions. 
# Gap file must consist of a tab-delimited file
#  with the following template
#  Fasta header, start coord, end coord
######################################################

use strict;
use FileHandle;
use File::Basename;

my $ref_file= $ARGV[0];  #fasta file
open (IN,$ARGV[1]); #gap file

my %gaps;
my $lastend=0;
my $laststart=0;
my ($start,$end);
my ($contig,$reference);
my ($header,@seq,$ref_sequence);
my %sequence;
my ($gap_start,$gap_end);
my $len_seq=0;

my $fh= FileHandle->new($ref_file)|| die "$!";
if ($fh->open("< $ref_file")){
   $/=">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless($_){next;};
      ($header,@seq)=split /\n/,$_;
      $ref_sequence= join "",@seq;
      $sequence{$header}=$ref_sequence;
#      print "$header\n";
      $len_seq= length $ref_sequence;
   }
   $/="\n";
   $fh->close;   
}

my ($name,$path,$suffix)=fileparse("$ARGV[1]",qr/\.[^.]*/);
while (<IN>){
   chomp;
   my ($ref,$gstart,$gend,$temp)=split /\t/,$_,4;
#   print "Start:  $gstart\nEnd:  $gend\n";
   $reference= $ref;
   $gaps{$gstart}=$gend;
}

open (OUT, ">$path/temp/$reference.fna");
my $first=1;
foreach my $gap(sort {$a<=>$b} keys %gaps){
#   print "$gap\t$gaps{$gap}\n";
   my $length= $gaps{$gap}-$gap;
   $gap_start= $gap;
   $gap_end= $gaps{$gap};
   if ($length >50){
      if ($gap_start>$lastend && $gap_end>$lastend){
#         print "$gap\t$gaps{$gap}\n";
         $lastend=$gap_end;
         if ($first){
            if ($gap_start==1){
#               $start=$gaps{$gap};
               if ($gap_end>100){$start=$gap_end-99;}
               elsif($gap_end<100){$start=1;}
            }
            else{
               $start=1;
#               $end=$gap;
               $end=$gap_start+99;
               $contig=$reference.'_'.$start.'_'.$end;
#               print "$start\t$end\n";
               my $output= substr($sequence{$reference},$start-1,$end-$start+1);
               my $length= length $output;
               if ($length>=250){print OUT ">$contig\n$output\n";}
            }
            $first=0;
         }
         else{
#            $end=$gap;
            $end=$gap+99;
            $contig=$reference.'_'.$start.'_'.$end;
            my $output= substr($sequence{$reference},$start-1,$end-$start+1);
#            print ">$contig\n$output\n";
            my $length= length $output;
            if ($length>=250){print OUT ">$contig\n$output\n";}
         }
#         $start=$gaps{$gap};
         $start=$gaps{$gap}-99;
      }
   }
}
if ($len_seq>$gap_end){
   $start= $gap_end;
   $end= $len_seq;
   $contig=$reference.'_'.$start.'_'.$end;
   my $output= substr($sequence{$reference},$start-1,$end-$start+1);
   print OUT ">$contig\n$output\n";
}

close IN;
close OUt;
exit;
