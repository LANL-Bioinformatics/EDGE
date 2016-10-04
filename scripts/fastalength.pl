#! /usr/bin/env perl 

my $seq_len=0;
my $id;
my @qual_line;
my $qual_len=0;
if (!@ARGV){print "perl $0 <fasta>\n";exit;}
while(<>)
{
    chomp;
     if(/^>(.*)/) 
     {  
        if ($seq_len){
             print "$id\t".$seq_len."\n";
        }
        if ($qual_len)
        { 
            print "$id\t".$qual_len."\n";
        }
        $id =$1;
        $seq_len = 0;
        $qual_len = 0;
     }
     else 
     {  
        if ($_ =~ /\d+/)
        {
          @qual_line=split;
          $qual_len += scalar(@qual_line);
        }
        else
        {
          $seq_len += length ($_);
        }
        
     }
}
       if ($seq_len){
             print "$id\t".$seq_len."\n";
       }
        if ($qual_len)
        { 
            print "$id\t".$qual_len."\n";
        }
