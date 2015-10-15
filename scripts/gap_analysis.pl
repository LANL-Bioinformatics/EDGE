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
          Options:   
          -gff <gff3_file>    GFF is an alternative to genbank. Use either one.   
            
};
my $debug=0;
my $time = time;
my ($Genbankfile, $Gap_file,$gff_file);
GetOptions('genbank=s'=>\$Genbankfile,
           'gff=s' => \$gff_file,
	       'gap=s'=>\$Gap_file,
           'verbose' =>\$debug,
	   )||die $usage;
die $usage if ((!$Genbankfile and !$gff_file ) or !$Gap_file);

my ($basename,$dir,$ext)= fileparse($Gap_file, qr/\.[^.]*/);

&print_timeInterval($time,"Loading Gap file");
open (IN, $Gap_file) or die "$!";
my @Gap_file=<IN>;
close IN;

if (!@Gap_file)
{
	print STDERR "No gap region found in file \"$Gap_file\"\n";
	exit;
}

my $inseq;
my $gffio;
my $scan_gap_r;
if ($Genbankfile){
    &print_timeInterval($time,"Loading Genbank");
    $inseq = Bio::SeqIO->new(-file => $Genbankfile);
    $scan_gap_r=&process_with_Genbank($inseq);
    $inseq->close();
}
elsif($gff_file)
{
    $scan_gap_r = &process_with_GFF($gff_file);
    #$scan_gap_r=&process_with_GFF($gffio);

}


print  "Chromosome\tgap_start\tgap_end\tgap_len\tMissing\tcds_start\tcds_end\tcds_product\n";
foreach my $gap_line(@Gap_file)
{
   chomp $gap_line;
   my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;  
   my $scan_key="$ref_id"."$gap_s"."$gap_end";
   next if ($gap_s !~ /\d+/);
   if ($scan_gap_r->{$scan_key})
   {
      print $scan_gap_r->{$scan_key};
   }
   else
   {
     print $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\t \t \t \t \n";
   }
}

&print_timeInterval($time,"Finished"); 

## END MAIN ##

sub process_with_GFF
{
    my $gff = shift;
    my $scan_gap;
    my $id;
    my $start;
    my $stop;
    my $old_id;
    my %seq_hash;
    my %gap_pos;
    open (my $fh, $gff) or die "Cannot open $gff\n";
    while(<$fh>)
    {
        chomp;
        if (/^#/)
        {
            if ($_ =~ /sequence-region\s+(\S+)\s+(\S+)\s+(\S+)/)
            {
                $id=$1;  $start=$2; $stop=$3;
               # &print_timeInterval($time,"Processing $id"); 
            }
            if (/FASTA/)
            {
                $/ = ">";
                while (my $line=<$fh>)
                {
                    $line =~ s/\>//g;
                    my ($id, @seq) = split /\n/, $line;
                    next if (!$id);
                    ($id) =~ s/^(\S+).*/$1/;
                    my $seq = join "", @seq;
                    my $len = length($seq);
                    $seq_hash{$id}=$len;
                }    
                $/="\n";
            }
        }
        else
        {
            my ($id,$source,$type,$start,$end,$score,$strand,$phase,$Attributes)=split /\t/,$_;
            
            if ($old_id ne $id){
            	&print_timeInterval($time,"Processing Gaps on $id");
            	delete $gap_pos{$old_id} if ($gap_pos{$old_id});
            	my @gaps = grep {$id} @Gap_file;
            	foreach my $gap_line(@gaps){
            		chomp $gap_line;
			my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;
			for($gap_s..$gap_end){$gap_pos{$_}=$gap_line;}
            	}
            }
            if (($type eq "CDS"  or $type =~ /RNA/) and $type ne "mRNA")
            {
                my %annotations=map { split /=/;} split /;/,$Attributes;
                my $product = $annotations{"product"} || $annotations{"Note"} ||  $annotations{"function"} || "Unknown" ;
                my $locus_tag = $annotations{"locus_tag"} || $annotations{"Name"} || "";
                $product = "$locus_tag:$product" if ($locus_tag);
                $product =~ s/\%2C/,/g;
                $product =~ s/\%3B/;/g;
                $scan_gap=&process_gap_file($product,$start,$end,\%gap_pos,$scan_gap);
            }
         
            $old_id=$id;    
        }
    }   
    close $fh;
    
    return ($scan_gap);
}


sub process_with_Genbank
{
my $inseq =shift;
my $n;
my $scan_gap;
my %gap_pos;
my $old_id;
while (my $seq = $inseq->next_seq){ 
  my $Locus_id = $seq->display_id();
  my $genome_seq = $seq->seq();
  my $genome_len = $seq->length;
  if ($old_id ne $Locus_id){
        &print_timeInterval($time,"Processing Gaps on $Locus_id");
        delete $gap_pos{$old_id} if ($gap_pos{$old_id});
        my @gaps = grep {$Locus_id} @Gap_file;
        foreach my $gap_line(@gaps){
		my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;
		for($gap_s..$gap_end){$gap_pos{$_}=$gap_line;}
        }
  }
  for my $feat_object ($seq->get_SeqFeatures){
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
         #my $nt_Seq = $feat_object->seq()->seq;
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
         $scan_gap=&process_gap_file($product,$start,$end,\%gap_pos,$scan_gap);
        # $scan_gap=&process_gap_file($product,$start,$end,$Locus_id,$scan_gap);
        
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
  	 $old_id=$Locus_id; 
} #end while my $seq
    return $scan_gap;
}

sub process_gap_file
{
    my $product = shift;
    my $feature_start = shift;
    my $feature_end = shift;
   # my $Locus_id = shift;
    my $gap_pos=shift;
    my $scan_gap = shift;
    my $cover;
    my $scan_key;
    my $print_string;
    my ($gap_s, $gap_end, $gap_len,$ref_id);
    my %features_pos;
    
    if ($gap_pos->{$feature_start} && ($gap_pos->{$feature_start} eq $gap_pos->{$feature_end})){
    	($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_pos->{$feature_start};
    	$cover = "Complete";
    	$scan_key="$ref_id"."$gap_s"."$gap_end";
   	$print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\t".$cover."\t".$feature_start."\t".$feature_end."\t".$product."\n";
   	$scan_gap->{$scan_key}.=$print_string;
    }else{
        $cover = "Partial";
    	my $gap_pos_record;
    	my $old_gap_pos_record;
    	my @covered_gaps;
    	for ($feature_start..$feature_end){
    		if ($gap_pos->{$_} ne $old_gap_pos_record){
    			$gap_pos_record = $gap_pos->{$_};
    			push @covered_gaps, $gap_pos_record;
    		}
    		$old_gap_pos_record=$gap_pos_record;
    	}
    	foreach my $gap (@covered_gaps){
    		($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap;
    		$scan_key="$ref_id"."$gap_s"."$gap_end";
   		$print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\t".$cover."\t".$feature_start."\t".$feature_end."\t".$product."\n";
   		$scan_gap->{$scan_key}.=$print_string;
    	}
    }

    
    #foreach my $gap_line(@Gap_file)
    #{
    #   chomp $gap_line;
    #   my ($gap_s, $gap_end, $gap_len,$ref_id) = split /\s+/,$gap_line;
    #   my $scan_key="$ref_id"."$gap_s"."$gap_end";
    #   my $print_string;
    #   next if ($gap_s !~ /\d+/);
	    #next if ($scan_gap{$scan_key} =~ /Complete/ or $scan_gap{$scan_key} =~ /Partial/);
    #   if ($ref_id =~ /$Locus_id/i or $basename =~ /$Locus_id/i){
    #     if ($gap_s <= $start && $gap_end >= $end)
    #     {
    #       $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tComplete\t".$start."\t".$end."\t".$product."\n";
    #       $scan_gap->{$scan_key}.=$print_string;
    #     }
    #     elsif($gap_s >= $start && $gap_end <= $end)
    #     {
    #       $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
    #       $scan_gap->{$scan_key}.=$print_string;
    #     }
    #     elsif($gap_s >= $start && $gap_end > $end && $gap_s <= $end )
    #     {
    #       $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
    #       $scan_gap->{$scan_key}.=$print_string;
    #     }
    #     elsif($gap_s < $start && $gap_end >= $start && $gap_end <= $end)
    #     {
    #       $print_string= $ref_id."\t".$gap_s. "\t".$gap_end."\t".$gap_len."\tPartial\t".$start."\t".$end."\t".$product."\n";
    #       $scan_gap->{$scan_key}.=$print_string;
    #     }
    #     else
    #     {
	      #  next if ($scan_gap{$scan_key} =~ /Complete/ or $scan_gap{$scan_key} =~ /Partial/);
         #  $print_string= $Locus_id."\t".$gap_s. "\t".$gap_end."\t \t \t \t \n";
         #  $scan_gap{$scan_key}=$print_string;
    #     }
    #   }
    #}
    return $scan_gap;
}

sub print_timeInterval{
    my $now = shift;
    my $msg = shift;
    $now = time - $now;
    my $string=sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
    print STDERR "[$string]  $msg\n";
}

