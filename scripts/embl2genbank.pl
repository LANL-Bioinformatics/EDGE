#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

if (@ARGV != 2) {    die "USAGE: embl2genbank.pl  embl_Iutput Genbank_Onput \n"; }

my $seqio = Bio::SeqIO->new('-format' => 'embl', '-file' => "$ARGV[0]");
my $seqout = new Bio::SeqIO('-format' => 'genbank', '-file' => ">>$ARGV[1]");
while( my $seq = $seqio->next_seq) {
  my $locus = $seq->display_id;
  $locus =~ s/\.final$//;
  $locus =~ s/^(\S+?)\.//;
  $seq->accession_number($locus);
  $seq->display_id($locus);
  $seqout->write_seq($seq)
}
