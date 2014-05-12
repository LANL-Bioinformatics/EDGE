#!/usr/bin/env perl
# Purpose: prepare files for phage finder. 
# This script takes a GenBank file from Prokka as input, and produces a
# phage_finder_info.txt (protein table), a prefix.con nucleotide sequence and prefix.faa protein sequences. 
# A PTT file is a line based, tab separated format with fixed column types.
#
# Written by Chien-Chi Lo
# 6 September 2013

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

my $outDir;
my $version=0.1;
my $prefix="Assembly";
GetOptions(
            "o=s"              => \$outDir,
            "p=s"              => \$prefix,
            "version"          => sub{print "Version $version\n";exit;},
            "help|?"           => sub{Usage()} );


if (@ARGV != 1) {&Usage();} 
unless ( -e $ARGV[0] ) { print "File not exist\n"; &Usage();}

my $gbk = Bio::SeqIO->new(-file=>$ARGV[0], -format=>'genbank');
my ($file_name, $path, $suffix)=fileparse("$ARGV[0]", qr/\.[^.]*/);

open (my $nt_fh, ">$outDir/$prefix.con") or die "Cannot write $outDir/$prefix.con\n";
open (my $aa_fh, ">$outDir/$prefix.pep") or die "Cannot write $outDir/$prefix.pep\n";
open (my $ph_fh, ">$outDir/phage_finder_info.txt") or die "Cannot write $outDir/phage_finder_info.txt\n";
while (my $seq = $gbk->next_seq)  # each LOCUS
{
    my $Locus_id = $seq->display_id();
    my $nt_seq = $seq->seq();
    
 #   system ("mkdir -p $outDir/$Locus_id");

    print $nt_fh ">".$Locus_id." Len=".$seq->length."\n".&fold_seq($nt_seq)."\n";

  
    my @cds = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
    #print $ph_fh $seq->display_id(), " - 0..",$seq->length,"\n";
    #print $ph_fh scalar(@cds)," proteins\n";
    #print $ph_fh join("\t", qw(Location Strand Length PID Gene Synonym Code COG Product)),"\n";

    my $cds_count=0;
    for my $f (@cds) {
       my $gi = '-';
       $gi = $1 if tag($f, 'db_xref') =~ m/\bGI:(\d+)\b/;
       my $cog = '-';
       $cog = $1 if tag($f, 'product') =~ m/^(COG\S+)/;
       my $locus_tag = tag($f, 'locus_tag');
       my $aa_seq = tag($f, 'translation');
    # ptt format
    #   my @col = (
     #    $f->start.'..'.$f->end,
     #    $f->strand >= 0 ? '+' : '-',
      #   ($f->length/3)-1,
      #   $gi eq "-"?$locus_tag:$gi,
      #   tag($f, 'gene'),
      #   $locus_tag,
      #   $cog,
      #   tag($f, 'product'),
      # );
       my @col =(
          $Locus_id,
          $seq->length,
          $locus_tag,
          $f->strand >= 0 ?$f->start:$f->end,
          $f->strand >= 0 ?$f->end:$f->start,
          tag($f, 'product'),
       );
       print $ph_fh join("\t", @col), "\n";
       my ($aa_id) = ($gi eq "-")?$locus_tag:"gi\|$gi\|$locus_tag";
       print $aa_fh ">".$aa_id." ".tag($f, 'product')."\n". &fold_seq($aa_seq)."\n";
       $cds_count++;
    }
  #  print $Locus_id, "\t", $cds_count,"\n";
#    if ( -z "$outDir/$Locus_id/$Locus_id.pep")
 #   {
       # no protein encoding in the contig.
  #     system("rm -rf $outDir/$Locus_id/");
  #  }
}
close $ph_fh;
close $aa_fh;
close $nt_fh;

sub tag {
   my($f, $tag) = @_;
   return '-' unless $f->has_tag($tag);
   return join(' ', $f->get_tag_values($tag));
}


sub fold_seq
{
    my $seq=shift;
    $seq =~ s/ //g;
    $seq =~ s/\n//g;
    $seq =~ s/\r//g;
    $seq =~ s/(.{70})/$1\n/g; 
    chomp $seq;
    return $seq;
}

sub Usage 
{
    print <<"END";
    Usage: perl $0 -o outDir GenBank_File
    Version $version
    -o      Output directory.
    -p      Prefix [default: Assembly]
END
    exit;
}
