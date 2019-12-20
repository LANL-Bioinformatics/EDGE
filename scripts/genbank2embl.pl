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
  $seq = clean_up_comment($seq);
  $seqout->write_seq($seq);
}

sub clean_up_comment{
        my $inseq_obj = shift;
        my $anno_collection = $inseq_obj->annotation;
        my @annotations = $anno_collection->get_Annotations('comment');
        my $comment = Bio::Annotation::Comment->new;
        my $old_comment=$annotations[0]->display_text;
        $old_comment =~ s/\n/ /g;
        $comment->text("$old_comment");

        my $coll = new Bio::Annotation::Collection;
        $coll->add_Annotation('comment',$comment);
        $inseq_obj->annotation($coll);
        return $inseq_obj;
}                   
