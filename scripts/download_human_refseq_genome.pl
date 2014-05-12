#! /usr/bin/env perl
use strict;
use LWP::Simple;

my $output=$ARGV[0];
if (!$output) {print "perl $0 output_dir\n";exit;}
system("mkdir -p $output");
while(<DATA>)
{
    chomp;
    my ($id,$acc)=split/\t/,$_;
    my $ref_fasta = &get_genbank_by_acc($acc);
    open (OUT,">$output/Human_Chromosome_$id.fasta") or die "Cannot write $output/Human_Chromosome_$id.fasta\n";
    print OUT $ref_fasta;
    close OUT;
}


sub get_genbank_by_acc{
  ## get reference genome genbank from NCBI.
  my $accession=shift;
  print STDERR "Getting reference genome ($accession) Fasta from NCBI...\n";
  my $utils = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils";
  my $efetch = "$utils/efetch.fcgi?"."db=nuccore&id=$accession&rettype=fasta&retmode=text";
  my $efetch_result = get($efetch);
  if ($efetch_result=~ /^\n/)
  {
     print "no \"$accession\" reference genbank file from NCBI genome.\n";
     exit(1);
  }
  return $efetch_result;
}

__DATA__
1	NC_000001
2	NC_000002
3	NC_000003
4	NC_000004
5	NC_000005
6	NC_000006
7	NC_000007
8	NC_000008
9	NC_000009
10	NC_000010
11	NC_000011
12	NC_000012
13	NC_000013
14	NC_000014
15	NC_000015
16	NC_000016
17	NC_000017
18	NC_000018
19	NC_000019
20	NC_000020
21	NC_000021
22	NC_000022
X	NC_000023
Y	NC_000024
MT	NC_012920
