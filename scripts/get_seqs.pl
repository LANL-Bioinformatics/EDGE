#! /usr/bin/perl
# retrieve the sequences of list(name) from fasta/q
use strict;
$|=1;
my $file1=$ARGV[0]; #list
my $file2=$ARGV[1]; #fasta file

if (scalar(@ARGV)!=2){
   print "Usage: perl $0 <list> <fasta/q>\n";
   exit;
}

my %hash;

if ( -f $file1 || $file1 eq "-") {
  my $fh;
  if (-r $file1){ 
        open ($fh,"$file1");
  }else{ $fh= *STDIN; }
  while (<$fh>)
  {
    chomp;
    $_ =~ s/ //g;
    $hash{$_}=1;
    $hash{"$_/1"}=1;
    $hash{"$_/2"}=1;
  }
  close $fh;
}else{
    $hash{$file1}=1;
    $hash{"${file1}/1"}=1;
    $hash{"${file1}/2"}=1;
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



