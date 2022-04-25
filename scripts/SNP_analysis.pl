#! /usr/bin/perl

use Bio::SeqIO;
use Bio::Tools::CodonTable;
use Getopt::Long;
use File::Basename;
use Cwd;
use strict;
use Data::Dumper;

my $workingDir=getcwd();
my $usage = qq{
Usage: $0 -genbank <ref_Genbankfile> -SNP <SNPs_file> -format [nucmer | vcf | changelog] -output <outputDir>
          # format "nucmer" is generated by show-snps -CT
          # format "vcf is Variant Call Format 4.1
          # format "changelog is consensus workflow output
          # note: make sure the genbank LOCUS match with the reference id in snps file. 
          Options:
          -gff <gff3_with_sequence_file>   GFF is an alternative to genbank. Use either one.   
             OR
          -gff <gff3_without_sequence_file>   -fasta  <genome_fasta_file>
          
};

my $time = time;
my $msg;
my $debug=0;
my ($Genbankfile, $SNPs_file, $format,$gff_file,$fasta_file);
my $outDir=$workingDir;
GetOptions('genbank=s'=>\$Genbankfile,
           'gff=s' => \$gff_file,
           'fasta=s' => \$fasta_file,
	       'SNP=s'=>\$SNPs_file,
	       'format=s'=>\$format,
           'output=s' =>\$outDir,
           'verbose' =>\$debug,
	   )||die $usage;
die $usage if ((!$Genbankfile and !$gff_file ) or !$SNPs_file );

my @SNPs_lines;
my $variant_count=&load_SNP();

if (!$variant_count)
{
	print STDERR "No variants found in file \"$SNPs_file\"\n";
	exit;
}

my $result;  # $result->{Locus_id}->{snp}, $result->{Locus_id}->{indel}
my $inseq;
my $myCodonTable   = Bio::Tools::CodonTable->new();
my $ambiguous_mode = ($SNPs_file =~ /w_ambiguous/)? 1 : 0;

if ($Genbankfile){
    &print_timeInterval($time,"Loading Genbank");
    $inseq = Bio::SeqIO->new(-file => $Genbankfile);
    $result = &process_with_Genbank($inseq);
    $inseq->close();
}
elsif($gff_file)
{
    $result = &process_with_GFF($gff_file,$fasta_file);
}


#### Output ####
my ($basename,$dir,$ext)= fileparse($SNPs_file, qr/\.[^.]*/);
&print_timeInterval($time,"Printing result to $basename.SNPs_report.txt and $basename.Indels_report.txt"); 
open (OUT ,">$outDir/$basename.SNPs_report.txt") or die "$!";
if ($format =~ /vcf/i){
    print OUT "Chromosome\tSNP_position\tRef_codon\tSub_codon\taa_Ref\taa_Sub\tSynonymous\tProduct\tCDS_start\tCDS_end\tCDS_strand\tCount\tDP\tRatio\tF:R\tRoot-mean-square-mapping_quality\n";
}else{
    print OUT "Chromosome\tSNP_position\tRef_codon\tSub_codon\taa_Ref\taa_Sub\tSynonymous\tProduct\tCDS_start\tCDS_end\tCDS_strand\n";
}

open (OUT2, ">$outDir/$basename.Indels_report.txt") or die "$!";
if ($format =~ /vcf/i)
{
   print  OUT2 "Chromosome\tINDEL_position\tRef_Seq\tIndel_seq\tLength\tType\tProduct\tCDS_start\tCDS_end\tCDS_strand\tFrameshift\tCount\tDP\tRatio\tF:R\tRoot-mean-square-mapping_quality\n";
}
else
{
   print  OUT2 "Chromosome\tINDEL_position\tSequence\tLength\tType\tProduct\tCDS_start\tCDS_end\tCDS_strand\tFrameshift\n";
}


foreach my $locus_id (keys %{$result})
{
  my %snps = %{$result->{$locus_id}->{snp}};
  my %indel = %{$result->{$locus_id}->{indel}};

  foreach my $pos (sort {$a <=> $b} keys %snps)
  {
       my $locus = $snps{$pos}->{Locus};
       my $ref_codon = $snps{$pos}->{ref_codon};
       my $snp_codon = $snps{$pos}->{snp_codon};
       my $product = $snps{$pos}->{product};
       my $ref_aa = $snps{$pos}->{ref_aa};
       my $snp_aa = $snps{$pos}->{snp_aa};
       my $synonymous = $snps{$pos}->{synonymous};
       my $strand = $snps{$pos}->{strand};
       $snp_aa =~ s/,$//;
       $synonymous =~ s/,$//;
       $snp_codon =~ s/,$//;
       
       print OUT  $locus,"\t",
                  $pos,"\t",
                  $ref_codon,"\t",
                  $snp_codon,"\t",
                  $ref_aa,"\t",
                  $snp_aa,"\t",
                  $synonymous,"\t",
                  $product,"\t",
                  $snps{$pos}->{cds_s},"\t",
                  $snps{$pos}->{cds_e},"\t",
                  $strand;
       if ($format =~ /vcf/i){
           print OUT "\t",
                  $snps{$pos}->{dpalt},"\t",
                  $snps{$pos}->{depth},"\t",
                  $snps{$pos}->{alt_ratio},"\t",
                  $snps{$pos}->{alt_f_r},"\t",
                  $snps{$pos}->{MQ};
       }
       print OUT "\n";
  }
  
  
  foreach my $pos(sort {$a <=> $b} keys %indel)
  {
      next if ($indel{$pos}->{skip});
      my $seq = ($indel{$pos}->{direction} == -1)?reverse($indel{$pos}->{seq}):$indel{$pos}->{seq};
      my $length = length($seq);
      my $product = $indel{$pos}->{prodcut};
      print  OUT2 $indel{$pos}->{Locus},"\t",
                 $pos, "\t";
      if ($format =~ /vcf/i)
      {
         my $ref_seq = $indel{$pos}->{ref_seq};
         my $indel_seq = $indel{$pos}->{indel_seq};
	 $length = abs (length($indel_seq)-length($ref_seq));
         print OUT2 $ref_seq,"\t",
                   $indel_seq,"\t",
                   $length,"\t";     
      }
      else
      {
         print OUT2 $seq,"\t",
                   $length,"\t";
      }       
      my $frameshift = ($length % 3 and $product !~ /Intergenic/i)? "Yes":"No";
	      
      print  OUT2 $indel{$pos}->{type},"\t",
             $product,"\t",
             $indel{$pos}->{cds_s},"\t",
             $indel{$pos}->{cds_e}."\t",
             $indel{$pos}->{strand}."\t",
	     $frameshift;
      if ($format =~ /vcf/i){
           print OUT2 "\t",
                  $indel{$pos}->{dpalt},"\t",
                  $indel{$pos}->{depth},"\t",
                  $indel{$pos}->{alt_ratio},"\t",
                  $indel{$pos}->{alt_f_r},"\t",
                  $indel{$pos}->{MQ};
      }
      print OUT2 "\n";
  }
}
close OUT;
close OUT2;

&print_timeInterval($time,"Finished"); 

### END MAIN ###

sub load_SNP {

	&print_timeInterval($time,"Loading SNP file");
	open (my $fh, $SNPs_file) or die "$!";
	my $header = <$fh>;
	$format = ($header =~ /VCF/i)? "vcf":"nucmer" if (! $format);
	my $count=0;
	while(<$fh>)
	{
		chomp;
		next if ($format =~ /nucmer/i && $_ !~ /^\d+/);
		next if ($format =~ /vcf|changelog/i && $_ =~ /^#/);
		push @SNPs_lines, $_;
		$count++;
	}
	close IN;
	return $count;
}


sub process_with_GFF 
{
    my $gff = shift;
    my $fasta = shift;
    my %result;
    my %seq_hash;
    my %coding_location;
    open (my $fh, $gff) or die "Cannot open $gff\n";
    %seq_hash=&readFastaFile(\%seq_hash,$fasta) if ($fasta);
    while(<$fh>)
    {
        chomp;
        if (/^#/)
        {
            if ($_ =~ /sequence-region\s+(\S+)\s+(\S+)\s+(\S+)/)
            {
                my $id=$1;  my $start=$2;  my $stop=$3;
                &print_timeInterval($time,"Loading $id GFF"); 
                if ($fasta and !$seq_hash{$id})
                {
                    foreach(keys %seq_hash)
                    {
                        $seq_hash{$id} = $seq_hash{$_} if ($_ =~ $id);
                    }
                }
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
                    $seq_hash{$id}=$seq;
                }    
                $/="\n";
            }
        }
        else
        {
            my ($id,$source,$type,$start,$end,$score,$strand,$phase,$Attributes)=split /\t/,$_;
            if ($type eq "CDS" or $type =~ /pseudogenic/)
            {
                my %annotations=map { split /=/;} split /;/,$Attributes;
                my $product = $annotations{"product"} || $annotations{"Note"} ||  $annotations{"function"} || "Unknown" ;
                my $locus_tag = $annotations{"locus_tag"} || $annotations{"Name"} || "";
                $product = "$locus_tag:$product" if ($locus_tag);
                $product =~ s/\%2C/,/g;
                $product =~ s/\%3B/;/g;
                for ($start..$end){$coding_location{$id}->{$_}="$start\t$end\t$strand\t$product";}
            }
        }
    }
    foreach my $Locus_id(keys %coding_location)
    {
        &print_timeInterval($time,"Processing SNP files"); 
        my $genome_seq = $seq_hash{$Locus_id};
        (my $snps,my $indel)=process_snp_file(\%coding_location,$genome_seq,$Locus_id);
       # print Dumper $indel;
        $result{$Locus_id}->{snp} = $snps;
        $result{$Locus_id}->{indel} = $indel;
    }
    return \%result;
}

sub process_with_Genbank
{
    my $inseq = shift;
    my %result;
    while (my $seq = $inseq->next_seq){  # read Genbank
    
      my $Locus_id = $seq->display_id();
      my %coding_location;
      my $genome_seq = $seq->seq();
      my $genome_len = $seq->length;
      &print_timeInterval($time,"Processing $Locus_id"); 
      for my $feat_object ($seq->get_SeqFeatures){
        my $locus_tag="";
        my $product=""; 
        my $aa_Seq="";
        my $note="";
        my @start;
        my @end;
        if ($feat_object->primary_tag eq "CDS"){
          if ( $feat_object->location->isa('Bio::Location::SplitLocationI'))
          {
             for my $location ( $feat_object->location->sub_Location() ) {
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
             # my $nt_Seq = $feat_object->seq()->seq;
              my $strand = $feat_object->location->strand;
              #my $aa_Seq_tmp = $feat_object->get_tag_values("translation");
              #$aa_Seq = join('',$feat_object->get_tag_values("translation")) if ($aa_Seq_tmp);
              eval { $product = join('',$feat_object->get_tag_values("product")); }; warn $@ if ($@ && $debug);
              eval { $note = join('',$feat_object->get_tag_values("note")); }; warn $@ if ($@ && $debug);
              eval { $locus_tag = join('',$feat_object->get_tag_values("locus_tag")); }; warn $@ if ($@ && $debug);
              $product = $note if (!$product);
              $product = "Unknown" if (!$product);
              $product = "$locus_tag:$product" if ($locus_tag);
              for ($start..$end){$coding_location{$Locus_id}->{$_}="$start\t$end\t$strand\t$product";}
          } #  foreach my $index(0..$#start)
        }#end if CDS
      }# end for  my $feat_object
      &print_timeInterval($time,"Processing SNP files"); 

      (my $snps,my $indel)=process_snp_file(\%coding_location,$genome_seq,$Locus_id);
      $result{$Locus_id}->{snp} = $snps;
      $result{$Locus_id}->{indel} = $indel;
    
    
    } #end while () Genbank
    return \%result;
}

sub process_snp_file 
{
   my $coding_location_r=shift;
   my $genome_seq=shift;
   my $Locus_id= shift;
   my %indel;
   my %snps;
   my $line_num=0;
   my %insertion_flag;
   my $insert_base=0;
   my $delet_base=1;
   foreach my $snps_line (@SNPs_lines)
   {
      $line_num++;
      chomp $snps_line;
      my $codon;
      my $ref_aa;
      my $ref_codon;
      my $snp_aa;
      my $snp_codon;
      my $tmp;
      my $snp_quality;
      my $vcf_info;
      my $vcf_info2;
      my $indel_flag=0;
      my ($dp, $dp_alt, $dp_alt_ratio, $dp_f_r, $mean_quality);
      my ($ref_pos,$ref_base,$snp,$snp_pos,$buff,$dist,$ref_direction,$snp_direction,$ref_id,$query_id);
      if ($format =~ /nucmer/i) 
      {
	 
          ($ref_pos,$ref_base,$snp,$snp_pos,$buff,$dist,$ref_direction,$snp_direction,$ref_id,$query_id)=split /\t/,$snps_line;
      }
      elsif ($format =~ /vcf/i)
      {
          ($ref_id,$ref_pos,$tmp, $ref_base,$snp,$snp_quality,$tmp,$vcf_info, $vcf_info2,$tmp)=split /\t/,$snps_line;
	  if ( $vcf_info =~ /DP4=(\d+),(\d+),(\d+),(\d+)/i) {
              $dp = $1 + $2 + $3 + $4;
              $dp_alt = $3 + $4;
              $dp_alt_ratio = sprintf("%.2f",$dp_alt/$dp);
              $dp_f_r = "$3:$4";
          }
	  my @SNPs_events = split/,/, $snp;
	  map { if (length($_) > 1) {$indel_flag = 1;} } @SNPs_events;
	  $indel_flag = 1 if length($ref_base) > 1 or $vcf_info =~ /INDEL/;
          if ($vcf_info =~ /MQ=(\d+)/){
		  $mean_quality = $1
	  }else{
		  $mean_quality = 'NA'
	  }
      }
      elsif ($format =~ /changelog/){
	  my ($RATIO_OF_REF,$HQ_DOC,$HQ_RATIO_OF_REF,$HQ_RATIO_OF_CON);
          my $amb_base="";
	  ($ref_id,$ref_pos,$ref_base,$snp,$dp,$RATIO_OF_REF,$dp_alt_ratio,$HQ_DOC,$HQ_RATIO_OF_REF,$HQ_RATIO_OF_CON)=split /\t/,$snps_line;
          if ($ambiguous_mode){
             $snp =~ s/\w\((\w)\)/$1/;
          }else{
             if ($snp =~ /(\w)\((\w)\)/){ $snp = $1; $amb_base = $2;} 
             $snp = $amb_base if ($amb_base and $ref_base eq $snp);
          }
          next if ($ref_base eq $snp);
	  $ref_base = "." if $ref_base eq '-';
	  $snp = '.' if $snp eq '-';
      }else
      { 
          print STDERR "Format is required and either nucmer(show-snps -CT) or vcf.\n";
          exit;
      }
      next if ($ref_id !~ /$Locus_id/i);
      next if $snp =~ /N/i and $format =~ /changelog/;
      if (defined $coding_location_r->{$Locus_id}->{$ref_pos}) # in CDS
      {
          my ($start,$end,$strand,$product)= split /\t/, $coding_location_r->{$Locus_id}->{$ref_pos};
          if ($ref_base eq "." or $snp eq "." or $indel_flag)  #insertion/deletion
          #if ($ref_base eq "." or $snp eq "." or $vcf_info =~ /INDEL/ or length ($snp) != length($ref_base))  #insertion/deletion
          {
             my $type;
             
             if ($ref_base eq '.' or length ($snp)>length($ref_base)){$type ="Insertion"};
             if ($snp eq '.' or length ($snp)<length($ref_base)){$type ="Deletion"};
             
             if ($format =~ /vcf/i)
             {
                 $indel{$ref_pos}->{ref_seq}=$ref_base;
                 $indel{$ref_pos}->{indel_seq}=$snp;
             }
             else
             { 
               if ($type eq "Deletion")
               {
                 if ($indel{$ref_pos})
                 {
                 }else{
                    if ($indel{$ref_pos-1})
                    {
                       $indel{$ref_pos}->{skip}=1;
                       $ref_pos = $ref_pos - $delet_base;
                       $delet_base++;
                    }else 
                    {
                       $delet_base=1;
                    }
                    $indel{$ref_pos}->{seq}.= $ref_base;
                    $indel{$ref_pos}->{direction}=$ref_direction;
                 }
               }
               else  # INSERTION
               { 
                # if (!$nucmerInsertSkip_r->{$snp_pos})
                # {
                   $indel{$ref_pos}->{seq}.=$snp;
                # }
                # $nucmerInsertSkip_r->{$snp_pos}=1;
                 $indel{$ref_pos}->{direction}=$snp_direction;
                 
               }
               
             }
             $indel{$ref_pos}->{type}=$type;
             $indel{$ref_pos}->{prodcut}=$product;
             $indel{$ref_pos}->{cds_s}=$start;
             $indel{$ref_pos}->{cds_e}=$end;
             $indel{$ref_pos}->{strand}=$strand;
             $indel{$ref_pos}->{Locus}=$ref_id;
             $indel{$ref_pos}->{depth}=$dp;
             $indel{$ref_pos}->{dpalt}=$dp_alt;
             $indel{$ref_pos}->{alt_ratio}=$dp_alt_ratio;
             $indel{$ref_pos}->{alt_f_r}=$dp_f_r;
             $indel{$ref_pos}->{MQ}=$mean_quality;
          }
          else #SNPs
          { 
             my @snps;
             $snps{$ref_pos}->{snp_codon}="";
             if (length($snp)>1) {@snps=split /,/,$snp;}else {@snps=$snp;}
             foreach my $snp (@snps){
               $snps{$ref_pos}->{snp_base}=$snp;
               $snps{$ref_pos}->{depth}=$dp;
               $snps{$ref_pos}->{dpalt}=$dp_alt;
               $snps{$ref_pos}->{alt_ratio}=$dp_alt_ratio;
               $snps{$ref_pos}->{alt_f_r}=$dp_f_r;
               $snps{$ref_pos}->{MQ}=$mean_quality;
               my $mod = ($ref_pos-$start+1) % 3;
               if ($mod % 3 == 0) # third base
               {
                  $codon=substr($genome_seq,$ref_pos-3,3);
                  $ref_codon=$codon;
                  if ($snps{$ref_pos-2})
                  {
                       if ($snps{$ref_pos-1})
                       {
                         $snp_codon=$snps{$ref_pos-2}->{snp_base}.$snps{$ref_pos-1}->{snp_base}.$snp;
                       }
                       else
                       {
                         $snp_codon=$snps{$ref_pos-2}->{snp_base}.substr($ref_codon,1,1).$snp;
                       }
                       #delete $snps{$ref_pos};
                       $snps{$ref_pos}->{ref_aa}="Merged with SNP ". ($ref_pos-2);
                       $snps{$ref_pos}->{product}="";
                       $snps{$ref_pos}->{ref_codon}=$ref_base;
                       $snps{$ref_pos}->{snp_codon}=$snp;
                       $snps{$ref_pos}->{Locus}=$ref_id;
                       $ref_pos = $ref_pos-2;
                  }elsif ($snps{$ref_pos-1})
                  {
                       #$snp_codon=$snps{$ref_pos-1}->{snp_base}.$snp.substr($ref_codon,2,1);
                       $snp_codon=substr($ref_codon,0,1).$snps{$ref_pos-1}->{snp_base}.$snp;
                       #delete $snps{$ref_pos};
                       $snps{$ref_pos}->{ref_aa}="Merged with SNP ".($ref_pos-1);
                       $snps{$ref_pos}->{product}="";
                       $snps{$ref_pos}->{ref_codon}=$ref_base;
                       $snps{$ref_pos}->{snp_codon}=$snp;
                       $snps{$ref_pos}->{Locus}=$ref_id;
                       $ref_pos = $ref_pos-1;
                  }
                  else
                  {
                    substr($codon,2,1,$snp);
                    $snp_codon=$codon;
                  }
                 # print $ref_codon,"\t",$snp_codon,"\n";
                }
                elsif ($mod % 3 == 1) # first base
                {
                  $codon=substr($genome_seq,$ref_pos-1,3);
                  $ref_codon=$codon;
                  substr($codon,0,1,$snp);
                  $snp_codon=$codon;
                 #print $ref_codon,"\t",$snp_codon,"\n";
                }
                elsif ($mod % 3 == 2) # sedond base
                {
                  $codon=substr($genome_seq,$ref_pos-2,3);
                  $ref_codon=$codon;
                  if ($snps{$ref_pos-1})
                  {
                       $snp_codon=$snps{$ref_pos-1}->{snp_base}.$snp.substr($ref_codon,2,1);
                       #delete $snps{$ref_pos};
                       $snps{$ref_pos}->{ref_aa}="Merged with SNP ".($ref_pos-1);
                       $snps{$ref_pos}->{product}="";
                       $snps{$ref_pos}->{ref_codon}=$ref_base;
                       $snps{$ref_pos}->{snp_codon}=$snp;
                       $snps{$ref_pos}->{Locus}=$ref_id;
                       $ref_pos = $ref_pos-1;
                  }
                  else
                  {
                       $snp_codon=substr($ref_codon,0,1).$snp.substr($ref_codon,2,1);
                  }
                  #print $ref_codon,"\t",$snp_codon,"\n";
                }
             
             if ($strand eq "-1" or $strand eq "-")
             { 
              $ref_codon=ReverseComplement($ref_codon);
              $snp_codon=ReverseComplement($snp_codon);
             } 
             $ref_aa= $myCodonTable->translate($ref_codon);
             $snp_aa= $myCodonTable->translate($snp_codon);
             my $codon_index =  int(($ref_pos - $start)/3) + 1;
             my $synonymous = ($ref_aa eq $snp_aa)? "Yes":"$ref_aa$codon_index$snp_aa";
             
             if (scalar (@snps)>1){
               $snps{$ref_pos}->{snp_aa}.=$snp_aa.",";
               $snps{$ref_pos}->{synonymous}.=$synonymous.",";
               $snps{$ref_pos}->{snp_codon}.=$snp_codon.",";
             }
             else
             {
               $snps{$ref_pos}->{snp_aa}=$snp_aa;
               $snps{$ref_pos}->{synonymous}=$synonymous;
               $snps{$ref_pos}->{snp_codon}=$snp_codon;
             }
             $snps{$ref_pos}->{ref_codon}=$ref_codon;
             $snps{$ref_pos}->{ref_aa}=$ref_aa;
             $snps{$ref_pos}->{product}=$product;
             $snps{$ref_pos}->{cds_s}=$start;
             $snps{$ref_pos}->{cds_e}=$end;
             $snps{$ref_pos}->{Locus}=$ref_id;
             $snps{$ref_pos}->{strand}=$strand;
            # print OUT $Locus_id,"\t",$ref_pos,"\t",$ref_codon,"\t",$snp_codon,"\t",$ref_aa,"\t",$snp_aa,"\t",$synonymous,"\t",$product,"\t",$start, "\t",$end,"\n";
            } # foreach my $snp (@snps){
          }
      } # end in CDS
      else   # not in CDS
      {
           if ($ref_base eq "." or $snp eq "." or $indel_flag)  #insertion/deletion
         # if ($ref_base eq "." or $snp eq "." or $vcf_info =~ /INDEL/ or length ($snp) != length($ref_base))  #insertion/deletion
          {
             my $type;
             
             if ($ref_base eq '.' or length ($snp)>length($ref_base)){$type ="Insertion"};
             if ($snp eq '.' or length ($snp)<length($ref_base)){$type ="Deletion"};
      
             if ($format =~ /vcf/i)
             {
                 $indel{$ref_pos}->{ref_seq}=$ref_base;
                 $indel{$ref_pos}->{indel_seq}=$snp;
             }
             else
             { 
               if ($type eq "Deletion")
               {
                 if ($indel{$ref_pos-1})
                 {
                   $indel{$ref_pos}->{skip}=1;
                   $ref_pos = $ref_pos - $delet_base;
                   $delet_base++;
                 }else 
                 {
                   $delet_base=1;
                 } 
                 $indel{$ref_pos}->{seq}.= $ref_base;
                 $indel{$ref_pos}->{direction}=$ref_direction;
               }
               else  # INSERTION
               {
                # next if ($nucmerInsertSkip_r->{$snp_pos});
                # print join ("\t",$start,$end,$ref_pos,$product),"\n";
                # $nucmerInsertSkip_r->{$snp_pos}=1;
                 $indel{$ref_pos}->{seq}.=$snp;
                 $indel{$ref_pos}->{direction}=$snp_direction;
               }
             }
             $indel{$ref_pos}->{type}=$type;
             $indel{$ref_pos}->{prodcut}="Intergenic region";
             $indel{$ref_pos}->{cds_s}="";
             $indel{$ref_pos}->{cds_e}="";
             $indel{$ref_pos}->{Locus}=$ref_id;
             $indel{$ref_pos}->{depth}=$dp;
             $indel{$ref_pos}->{dpalt}=$dp_alt;
             $indel{$ref_pos}->{alt_ratio}=$dp_alt_ratio;
             $indel{$ref_pos}->{alt_f_r}=$dp_f_r;
             $indel{$ref_pos}->{MQ}=$mean_quality;
          }
          else # not in CDS SNPs
          {  
        
             $snps{$ref_pos}->{snp_codon}=$snp;
             $snps{$ref_pos}->{snp_aa}="";
             $snps{$ref_pos}->{synonymous}="";
             $snps{$ref_pos}->{ref_codon}=$ref_base;
             $snps{$ref_pos}->{ref_aa}="";
             $snps{$ref_pos}->{product}="Intergenic region";
             $snps{$ref_pos}->{cds_s}="";
             $snps{$ref_pos}->{cds_e}="";
             $snps{$ref_pos}->{Locus}=$ref_id;
             $snps{$ref_pos}->{strand}="";
             $snps{$ref_pos}->{depth}=$dp;
             $snps{$ref_pos}->{dpalt}=$dp_alt;
             $snps{$ref_pos}->{alt_ratio}=$dp_alt_ratio;
             $snps{$ref_pos}->{alt_f_r}=$dp_f_r;
             $snps{$ref_pos}->{MQ}=$mean_quality;
          }
          
      }# not in CDS
     
   } # end  foreach (my $snps_line= <IN>)
             
    return (\%snps,\%indel);
}

sub readFastaFile
{
    my $seq_r = shift;
    my $file = shift;
    open (my $fh, $file) or die "Cannot open $file\n";
    $/ = ">";
    while (my $line=<$fh>)
    {
         $line =~ s/\>//g;
         my ($id, @seq) = split /\n/, $line;
         next if (!$id);
         ($id) =~ s/^(\S+).*/$1/;
         my $seq = join "", @seq;
         my $len = length($seq);
         $seq_r->{$id}=$seq;
    }    
    $/="\n";
    return %{$seq_r};
}

sub ReverseComplement{
        my $dna = $_[0];
        my $ReverseCompSeq = reverse ($dna);
        $ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
        return($ReverseCompSeq);
}

sub print_timeInterval{
    my $now = shift;
    my $msg = shift;
    $now = time - $now;
    my $string=sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
    print STDERR "[$string]  $msg\n";
}
