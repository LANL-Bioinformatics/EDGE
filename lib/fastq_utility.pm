#!/usr/bin/env perl

use strict;
use warnings;

package fastq_utility;

sub quality_encoding_coversion 
{  
    # given quality acsii string, input offset and output offset
    my $q_string=shift;
    my $input_offset=shift;
    my $out_offset=shift;
    $q_string=~ s/(\S)/chr(ord($1)-$input_offset+$out_offset)/eg;
    return($q_string);
}



sub low_complexity_filter
{
    my ($seq,$cut_off_ratio) = @_;
    my $seq_len = length ($seq);
    my @low_complex_array=("A","T","C","G","AT","AC","AG","TA","TC","TG","CA","CT","CG","GA","GT","GC");
    my $filter=0;
    for my $item (@low_complex_array)
    {
      my $item_len = length ($item);
      my $num_low_complex = $seq =~ s/$item/$item/g;
      if (($num_low_complex*$item_len/$seq_len) >= $cut_off_ratio)
      {
        $filter=1;
        last; #end for loop
      }
    }
    return ($filter);

}

sub is_fastq {
    #given a file it will return 1 for fastq, 0 for others.
    my $file=shift;
    my $fastq=0;
    open (IN, "< $file") or die "$!\n";
    my $header=<IN>;
    close IN;
    $fastq=1 if ($header =~ /^@/);
    return $fastq;
}

sub checkQualityFormat {
    # $offset_value=&checkQualityFormat($fastq_file)
    # $offset_value = -1 means not proper fastq format.
    my $fastq_file=shift;
    # open the files
    if ($fastq_file =~ /gz$/)
    {
        open FQ, "zcat $fastq_file | " or die $!;
    }
    else
    {
        open FQ, "<", $fastq_file or die $!;
    }

    # initiate
    my @line;
    my $l;
    my $number;
    my $offset;
    my $line_num=0;
    # go thorugh the file
    my $first_line=<FQ>;
    if ($first_line !~ /^@/) {$offset=-1; return $offset;}
    OUTER:while(<FQ>){
      $line_num++;
      # if it is the line before the quality line
      if($_ =~ /^\+/){
    
       $l = <FQ>; # get the quality line
       chomp $l;
       @line = split(//,$l); # divide in chars
       for(my $i = 0; $i <= $#line; $i++){ # for each char
        $number = ord($line[$i]); # get the number represented by the ascii char
      
        # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
        if($number > 104){
          $offset=33;  # pacbio CCS reads
          last OUTER;
        }
        elsif($number > 74){ # if solexa/illumina
          $offset=64;
          #die "This file is solexa/illumina format\n"; # print result to terminal and die
          # read a few more quality line to make sure the offset 
          last OUTER if ($line_num > 100);
        }elsif($number < 59){ # if sanger
          $offset=33;
          #die "This file is sanger format\n"; # print result to terminal and die
          last OUTER;
        }
       }
      }
    }
    return $offset;
}

1;

