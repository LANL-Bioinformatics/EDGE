#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use File::Basename;

if (@ARGV != 1) {    die "USAGE: genbank2embl.pl  Genbank_Input\n"; }
my $genbank=$ARGV[0];
my $seqio = Bio::SeqIO->new('-format' => 'genbank', '-file' => "$genbank");
my ($file_name, $file_path, $file_suffix)=fileparse("$genbank", qr/\.[^.]*/);
while( my $seq = $seqio->next_seq) {
  my $Locus_id = $seq->display_id();
  my $seqout = new Bio::SeqIO('-format' => 'embl', '-file' => ">$file_path/$Locus_id.embl");
  $seqout->write_seq($seq)
}
