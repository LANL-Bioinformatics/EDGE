#!/usr/bin/env perl
# require samtools 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin qw($RealBin);
use lib $RealBin;
use gi2lineage;

my $file1=$ARGV[0]; #list two columns, readsID, taxID
my $file2=$ARGV[1]; #fastq file
my $target_taxID=$ARGV[2];

if (scalar(@ARGV)!=3){
   print "Usage: perl $0 <list file> <fastq> <target_taxID>\n";
   print "\n         <list file> two columns, readsID<tab>taxID\n";
   exit;
}

my %hash;
loadTaxonomy();
my $target_lineage = taxid2lineage($target_taxID);
$target_lineage =~ s/\W/_/g;
#print $target_lineage,"\n";
if ( -f $file1 || $file1 eq "-") {
  my $fh;
  if (-r $file1){
        open ($fh,"$file1");
  }else{ $fh= *STDIN; }
  while (<$fh>)
  {
    chomp;
    my ($readsID,$taxID) = split /\t/;
    $taxID =~ s/ //g;
    my $lineage=taxid2lineage($taxID);
    $lineage =~ s/\W/_/g;
    if ($lineage =~ /$target_lineage/){
    #  print join("\t",$readsID,$taxID,$lineage),"\n";
      $hash{$readsID}=1;
      $hash{"$readsID/1"}=1;
      $hash{"$readsID/2"}=1;
    }
  }
  close $fh;
}


open (IN2,"$file2");
my $name;
my $name_with_head;
my $name_first_no_space_words;
my $name_first_no_space_words_with_head;;
my $fastq;
while (<IN2>)
{
   chomp;
   if ($_=~ /^(@\S+)/)
   {
           $name=$_;
           $name_with_head=$name;
           $name_first_no_space_words=$1;
           $name =~ s/\@//;
           $name_first_no_space_words_with_head=$name_first_no_space_words;
           $name_first_no_space_words =~ s/\@//;
           my $seq=<IN2>;
           my $qual_id=<IN2>;
           my $qual_seq=<IN2>;
           if ($hash{$name} || $hash{$name_with_head} || $hash{$name_first_no_space_words} || $hash{$name_first_no_space_words_with_head})
           {
              print "\@$name\n$seq$qual_id$qual_seq";
           }
   }
   elsif ($_=~ /^(>\S+)/)
   {
           $name=$_;
           $name_with_head=$name;
           $name_first_no_space_words=$1;
           $name =~ s/^>//;
           $name_first_no_space_words_with_head=$name_first_no_space_words;
           $name_first_no_space_words =~ s/^>//;
           print $_."\n" if ($hash{$name} || $hash{$name_with_head} || $hash{$name_first_no_space_words} || $hash{$name_first_no_space_words_with_head});
   }
   else
   {
           print $_."\n" if ($hash{$name} || $hash{$name_with_head} || $hash{$name_first_no_space_words} || $hash{$name_first_no_space_words_with_head});
   }
}

close IN2;


