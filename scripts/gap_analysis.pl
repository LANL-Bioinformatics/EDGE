#! /usr/bin/perl

use Bio::SeqIO;
use Getopt::Long;
use File::Basename;
use strict;

my $usage = qq{
Usage: $0 -genbank <Genbankfile> -gap <Gap_file>
          the gap file is four columns tab-delimited file.
             gap_start, gap_end, gap_length, genome_id
             The genome_id should match with genbank LOCUS.
            
};
my $debug=0;
my ($Genbankfile, $Gap_file);
GetOptions('genbank=s'=>\$Genbankfile,
	       'gap=s'=>\$Gap_file,
           'verbose' =>\$debug,
	   )||die $usage;
die $usage if (!$Genbankfile or !$Gap_file);

my ($basename,$dir,$ext)= fileparse($Gap_file, qr/\.[^.]*/);

open (IN, $Gap_file) or die "$!";
my @Gap_file=<IN>;
close IN;

my %scan_gap;
my %CDS;
my $inseq = Bio::SeqIO->new(-file => $Genbankfile);
my $n;
while (my $seq = $inseq->next_seq){ 
  my $Locus_id = $seq->display_id();
  my $genome_seq = $seq->seq();
  my $genome_len = $seq->length;
  my ($output_aa, $output_contig, $fig_id,
      $genome_id, $project, $fig_contig_id, $contig_id);
  for my $feat_object ($seq->get_SeqFeatures){
    my $fig_id;
    my $product=""; 
    my $locus_tag="";
    my $aa_Seq="";
    my $note="";
    my @start;
    my @end;
    if ($feat_object->primary_tag eq "CDS" or $feat_object->primary_tag =~ /RNA/i){
      $n++;
      if ( $feat_object->location->isa('Bio::Location::SplitLocationI'))
      {
         for my $location ( $feat_object->location->sub_Location ) {
             push @start, $location->start;
             push @end, $location->end; 
         }
      }
      else
      {
          push  @start, $feat_object->location->start;       
          push  @end,  $feat_object->location->end;
      }
      foreach my $index(0..$#start)
      {
         my $start = $start[$index];
         my $end = $end[$index];
         my $nt_Seq = $feat_object->seq()->seq;
         my $strand = $feat_object->location->strand;
          #print $start,"\t",$end,"\t",$strand,"\n";
         eval { $product = join('',$feat_object->get_tag_values("product")); }; warn $@ if ($@ && $debug);
         #$aa_Seq = join('',$feat_object->get_tag_values("translation"));
         eval { $note = join('',$feat_object->get_tag_values("note")); }; warn $@ if ($@ && $debug);
         eval { $locus_tag = join('',$feat_object->get_tag_values("locus_tag")); }; warn $@ if ($@ && $debug);
         $product = $note if (!$product);
         $product = "Unknown" if (!$product);
         $product = "$locus_tag:$product" if ($locus_tag);
       #  print $start,"\t",$end,"\t",$strand,"\n";
         $CDS{$Locus_id}->{start}=$start;
         $CDS{$Locus_id}->{end}=$end;
         $CDS{$Locus_id}->{strand}=$strand;
   
         foreach my $gap_line(@Gap_file)
         {
            chomp $gap_line;
            my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;
            my $scan_key="$gap_s"."$gap_end";
            my $print_string;
            next if ($gap_s !~ /\d+/);
	    #next if ($scan_gap{$scan_key} =~ /Complete/ or $scan_gap{$scan_key} =~ /Partial/);
            if ($ref_id =~ /$Locus_id/i or $basename =~ /$Locus_id/i){
              if ($gap_s <= $start && $gap_end >= $end)
              {
                $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tComplete\t".$start."\t".$end."\t".$product."\n";
                $scan_gap{$scan_key}.=$print_string;
              }
              elsif($gap_s >= $start && $gap_end <= $end)
              {
                $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
                $scan_gap{$scan_key}.=$print_string;
              }
              elsif($gap_s >= $start && $gap_end > $end && $gap_s <= $end )
              {
                $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
                $scan_gap{$scan_key}.=$print_string;
              }
              elsif($gap_s < $start && $gap_end >= $start && $gap_end <= $end)
              {
                $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
                $scan_gap{$scan_key}.=$print_string;
              }
              else
              {
	      #  next if ($scan_gap{$scan_key} =~ /Complete/ or $scan_gap{$scan_key} =~ /Partial/);
              #  $print_string= $Locus_id."\t".$gap_s. "\t".$gap_end."\t \t \t \t \n";
              #  $scan_gap{$scan_key}=$print_string;
              }
            }
         }
      } # end foreach $start
    } # end if CDS/RNA
    if ($feat_object->primary_tag eq "source"){
      #$genome_id = join('',$feat_object->get_tag_values("genome_id"));
      #$project = join('',$feat_object->get_tag_values("project"));
      #($contig_id) = $Locus_id =~ /(\S+)/;
      #$fig_contig_id = "$genome_id.contig.$contig_id";
      #$output_contig = ">$fig_contig_id $project\n$contig_seq\n";
    }
  } # end for  my $feat_object
} #end while my $seq

close OUT;

print  "Chromosome\tgap_start\tgap_end\tgap_len\tMissing\tcds_start\tcds_end\tcds_product\n";
foreach my $gap_line(@Gap_file)
{
   chomp $gap_line;
   my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;  
   my $scan_key="$gap_s"."$gap_end";
   next if ($gap_s !~ /\d+/);
   if ($scan_gap{$scan_key})
   {
      print $scan_gap{$scan_key};
   }
   else
   {
     print $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\t \t \t \t \n";
   }
}



