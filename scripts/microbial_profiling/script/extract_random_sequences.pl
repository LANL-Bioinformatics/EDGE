#! /usr/bin/perl
# chienchi at lanl.gov
# 20110412

use strict;
use Getopt::Long;

my $file;
my $output;
my $output_num;
GetOptions('o=s'    => \$output,
           'i=s'    => \$file,
           'n=f'    => \$output_num,
           'help|?' => sub{Usage()}
           );  

&Usage unless ($output && $file && $output_num);

my $total_seq_num=&seq_count($file);
die "The extract number($output_num) is bigger/equal than sequence total number($total_seq_num).\n" if ($output_num >= $total_seq_num); 

if ($output_num>0 and $output_num<1)
{
    $output_num = int ($output_num * $total_seq_num);
}

$output_num = int ($output_num);
print "Extracting random $output_num sequences out of total $total_seq_num\n";

my %get;
while (1)
{
    $get{int(rand($total_seq_num))+1}=1;
    my $num_of_random=scalar (keys %get);
    last if ($num_of_random == $output_num);
}
print STDERR "Done choosing $output_num sequences.\n";

my $name;
my $fastq;
my $count;

open (OUT, ">$output") or die "Cannot write $output\n";
open (IN, "$file") or die "Cannot open $file\n";
while(<IN>)
{
   chomp;
   if ($_=~ /^@/)
   {
           $count++;
           $name=$_;
           my $seq=<IN>;
           $seq =~ s/\n//g;
           while ($seq !~ /\+/)
           {
             $seq .= <IN>;
             $seq =~ s/\n//g;
           }
           my $q_name_pos=index($seq,"+");
           my $qual_id = substr($seq,$q_name_pos);
           $seq = substr($seq, 0, $q_name_pos);
           my $seq_len = length $seq;                                                                                                                    my $qual_seq=<IN>;
           $qual_seq =~ s/\n//g;
           my $qual_seq_len = length $qual_seq;
           while ( $qual_seq_len < $seq_len )
           {
              last if ( $qual_seq_len == $seq_len);
              $qual_seq .= <IN>;
              $qual_seq =~ s/\n//g;
              $qual_seq_len = length $qual_seq;
           }
           my $print_out= "$name\n$seq\n$qual_id\n$qual_seq\n";

           if ($get{$count})
           {
              print OUT $print_out;
           }
           
   }
   elsif ($_=~ /^>/)
   {
       $count++;
       print OUT $_."\n" if ($get{$count});
       
   }
   else 
   { 
       print OUT $_."\n" if ($get{$count});
   }
}
close IN;
close OUT;

sub Usage 
{
  print "perl $0 -o <output> -i <fasta/q> -n # \n";
  print "     -i    input fasta/q file\n";
  print "     -o    output fasta/q file\n";
  print "     -n    number of sequence to extract [int]\n";
  print "           it can be a fraction of total sequence [0.0 to 1.0]\n";
  print "     -h    print usage\n";
  exit;
}

sub seq_count
{
  my $in=$_[0];
  my $count;
  my $mode = "fasta";
  open (IN, "$in") or die "open cannot open $in\n";
  while(<IN>)
  {
   chomp;
   $mode = "fastq" if ($_ =~ /^@/);
  }
  close IN;

  my $lc = `wc -l $in`;
  my ($count) = $lc =~ /^(\d+)/;

  $count = $count/4 if $mode eq "fastq";

 return ($count);
}
