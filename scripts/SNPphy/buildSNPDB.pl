#!/usr/bin/perl -w
################################################################
# Written by Sanaa Ahmed
# Jan. 03, 2013

# Given a directory containing snp files creates SNP alignment 
# SNP files can be:  
#  vcf files from read mapping (suffix: .vcf)
#  show-snps file from nucmer  (suffix .snps)
# Requires:
#  Gap tab-delimited file with the following template
#   Fasta header, start coord, end coord
#  Reference fasta file used for read mapping and nucmer
#  List of all genomes to report SNP alignment for.
################################################################

use strict;
use Getopt::Long;
use diagnostics;
use FileHandle;
use File::Basename;

my ($in_dir,$snp_file,$ref_file,$ref_header);
my ($ref,$rpos,$rbase,$query,$qpos,$qbase,$buff,$dist,$frm);
my (%snp_location,@header_list,$ref_sequence,%gap_location);
my ($reference,$gap_fasta,@fragments);
my $contig=0;
my %positions;
my %noncoding_location;
my %coding_location;
my %coding_snp;
my $outfile;
my $genbank;
my $cdsSNPoutfile;
my $cdsSNPstats;
my $CDScoords;
my $noncoding;

GetOptions(
   'i=s'      => \$in_dir,
   'r=s'      => \$ref_file,
   'l=s'      => \$ref_header,
   'o=s'      => \$outfile,
   'g=i'      => \$genbank,
   'help|h'   => sub{usage()},
);

if (!$in_dir && !$ref_file && !$ref_header){die &usage;}
if ($in_dir=~ /.+\/$/){my $tmp= chop($in_dir);}
my ($name,$path,$suffix)=fileparse("$ref_file",qr/\.[^.]*/);

my $allSNPoutfile= $in_dir."/${outfile}_all.aln.fasta";
my $allSNPstats= $in_dir.'/snp_stats_all.txt';
#my $allSNPstats= $in_dir.'/snp_stats_all_'.$name.'.txt';
my $allSNPfile= $in_dir.'/snp_comparison_all.txt';
if ($genbank==1){
   $cdsSNPoutfile= $in_dir."/${outfile}_cds.aln.fasta";
   $CDScoords=$in_dir.'/CDScoords.txt';
   $noncoding= $in_dir.'/noncoding.txt';
   $cdsSNPstats=$in_dir.'/snp_stats_cds.txt';
}
my $ambiguous= $in_dir.'/ambiguousSNPpositions.txt';
my $gapfile= $in_dir.'/gaps.txt';
my $basefile=$in_dir.'/basesUsed.txt';
my $gap_file= $in_dir.'/all_gaps.txt';

open (OUT, ">$allSNPoutfile")|| die "$!";
open (STAT, ">$allSNPstats")|| die "$!";
open (AMB, ">$ambiguous")|| die "$!";
#open (COMP, ">$allSNPfile")|| die "$!";
open (GAPF, ">$gapfile")|| die "$!";
open (BASE, ">$basefile")|| die "$!";
if ($genbank==1){
   open (CDSOUT, ">$cdsSNPoutfile")|| die "$!";
   open (CDSSTAT,">$cdsSNPstats")||die "$!";
}

print AMB "Reference\tQuery\tSNP\tSNP\tref Pos.\n";
if ($genbank==1){print STAT "Strain1\tStrain2\tType\tSNP1 pos\tSNP2 pos\tSNP1\tSNP2\tCDS start\tCDS end\n";}
else{print STAT "Strain1\tStrain2\tSNP1 pos\tSNP2 pos\tSNP1\tSNP2\n";}
#print COMP "SNP Pos\tGenomes\n";
#if ($genbank==1){print CDSSTAT "SNP Pos\tGenomes\n";}
if ($genbank==1){print CDSSTAT "Strain1\tStrain2\tSNP1 pos\tSNP2 pos\tSNP1\tSNP2\tCDS start\tCDS end\n";}
if ($genbank==1){read_coords();}
read_reference($ref_file);
print "Reading Gaps file.\n";
my $skip_query_ref=read_gap($gap_file);
print "Printing Gap coords.\n";
read_directory($in_dir);
print_gap_coords();
#check_positions(%positions);
create_snp_array();
#create_comparison();
print "SNP alignment complete.\n";

sub read_reference
{
my ($header,@seq);
my $fh= FileHandle->new($ref_file)|| die "$!";
if ($fh->open("< $ref_file")){
   $/=">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless($_){next;};
      ($header,@seq)=split /\n/,$_;
      $ref_sequence= join "",@seq;
      print length $ref_sequence,"\n";
      $reference= $header;
   }
   $/="\n";
   $fh->close;
}

open (IN,"$ref_header")||die "$!";
while (<IN>){
   chomp;
   $_ =~ s/(.+_read)\.\w+/$1/;
   push(@header_list,"$reference:$_");
}
close IN;
}

sub read_coords
{

open (CDS, "$CDScoords")||die "$!";
while (<CDS>){
   chomp;
   my ($id,$start,$end,$product)= split /\s+/,$_;
   for ($start..$end){$coding_location{$_}="$start,$end,$product";}
  # if ($snp>=$start && $snp<=$end){return ($start,$end);}
}
}

sub read_directory
{
my $snp_dir= $in_dir.'/snps';
my $snp_count;
if (-e '$in_dir/*.gaps'){`mv $in_dir/*.gaps $in_dir/gaps`;}
if (-e '$in_dir/*.vcf'){`mv $in_dir/*.vcf $snp_dir`;}
opendir(DIR,$snp_dir);
while (my $files= readdir(DIR)){
   next if ($files=~ /^..?$/);
   if ($files=~ /^${name}_(.+)\.snps$/){
      my $query=$1;
      if ($files=~ /contigs.snps$/){
         $snp_file= $snp_dir.'/'.$files;
         contig_nucmer_snp($snp_file);
      }
      else{
         $snp_file= $snp_dir.'/'.$files;
         $snp_count=read_nucmer_snp($snp_file);
         #if ($snp_count==0){
             #$skip_query{$query}=1;
             #print "No snps $query to $name. Will not use $query to build SNP phylogeny\n";
         #}
      }
   }
   elsif ($files=~/.+\.vcf/){
      $snp_file= $snp_dir.'/'.$files;
      read_vcf_snp($snp_file);
   }
}
}

sub read_gap
{
my $count=0;
my %skip_query;
my ($gap_start,$gap_end,$gap_length);
open (IN, "$gap_file")|| die "$!";
while (<IN>){
   chomp;
   if (/^$reference\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)$/){
      ($gap_start,$gap_end,$gap_length,$query)=($1,$2,$3,$4);
      if ( ($gap_length/(length $ref_sequence)) > 0.75)
      {
          $skip_query{$query}=1;
          print "$query are too different >75% to $reference. Will not use in building SNP tree\n";
          next;
      }
      for (my $i=$gap_start;$i<=$gap_end;$i++){$gap_location{$i}++;}
   }
}
close IN;

#if (-e $in_dir.'/noncoding.txt'){
#   open (IN, "$noncoding")|| die "$!";
#   while (<IN>){
#      chomp;
#      if (/^$reference\s+(\d+)\s+(\d+)\s+\S+$/){
#         my ($start,$end)=($1,$2);
#         for (my $i=$start;$i<=$end;$i++){$noncoding_location{$i}++;}
#      }
#   }
#}
#close IN;
return \%skip_query;
}

sub print_gap_coords
{
my $last=0;
my $first=1;
my $gap_total=0;

my ($start,$end);
foreach my $keys(sort{$a<=>$b}keys %gap_location){
#   print "$keys\t$start\t$last\n";
   if ($first){$start=$keys;$last=$keys;$first=0;}
   elsif ($last!=$keys-1){
#      print "$start\t$last\n";
      print GAPF "$start\t$last\n";
      my $gap_length=$last-$start+1;
      $gap_total+=$gap_length; 
      $start=$keys;
   }
   $last=$keys;
}
my @locations=keys %gap_location;
$gap_total = $#locations + 1;
my $base_total= (length $ref_sequence)-$gap_total;
print BASE "Total sequence length:\t",length $ref_sequence,
"\nTotal gap length:\t$gap_total
Total bases used in SNP determination:\t$base_total\n";

die "ERROR: No core SNP region to build SNP alignment\n" if ($gap_total >= length ($ref_sequence));

close GAPF;
}

sub contig_nucmer_snp
{
my $count =0;
open (IN,"$snp_file")||die "$!";
while (<IN>){
   chomp;
   if (/^\d/){
      $count=1;
      ($rpos,$rbase,$qbase,$qpos,$buff,$dist,$frm,$frm,$ref,$query)= split ("\t",$_);
      if ($query=~/(\S+)_\d+$/){$query= $1.'_contig';}
      if ($qbase!~ /\./ && $rbase!~ /\./ && $ref=~ /$reference/){
         $snp_location{$rpos}{"$ref:$query"}= $qbase;
         $snp_location{$rpos}{"$ref:$ref"}=$rbase;
         $positions{$rpos}{"$ref:$query"}=$qpos;
      }
      else {$gap_location{$rpos}++;}
   }
}
close IN;
}
sub read_nucmer_snp
{
my $count =0;
open (IN,"$snp_file")||die "$!";
while (<IN>){
   chomp;
   if (/^\d/){
      $count=1;
      ($rpos,$rbase,$qbase,$qpos,$buff,$dist,$frm,$frm,$ref,$query)= split ("\t",$_);
      if ($query=~ /(.*)\_(\d+)\_\d+$/){
         $query=$1;
         $qpos=$qpos+$2-1;
      }
      if ($ref=~ /(.*)\_(\d+)\_\d+$/){
         $ref=$1;
         $rpos=$rpos+$2-1;
      }
      if ($qbase!~ /\./ && $rbase!~ /\./ && $ref=~ /$reference/){
         $snp_location{$rpos}{"$ref:$query"}= $qbase;
         $positions{$rpos}{"$ref:$query"}=$qpos;
      }
      else {$gap_location{$rpos}++;}
   }
}
close IN;
return $count;
}

sub read_vcf_snp
{
my ($id,$qual,$filter,$info,$info2,$tmp);
my %mapped;
open (IN,"$snp_file")||die "$!";
while (<IN>){
   chomp;
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /Users/218819/Projects/edge_run/testData/output/SNP_Phylogeny/Ecoli/results/Ecoli_55989_output.sort.bam
   if (/^#CHROM.+\/\S+\/$reference\_(\S+)\.sort\.bam/){$query=$1.'_read';}
   if (!/^\#/){
      ($ref,$rpos,$id,$rbase,$qbase,$qual,$filter,$info,$info2,$tmp)= split ("\t",$_);
      if ($qbase!~ /\./ && $rbase!~ /\./ && $info!~ /INDEL/){
         if ($ref=~ /(.+)_(\d+)_(\d+)$/){
            $ref= $1;
            my $start=$2;
            if ($start>1){$rpos= $rpos+$start-1;}
         }
         if ($qbase=~ /(\D),(\D)/){
            $qbase=$1; my $temp=$2;
            print AMB "$reference\t$query\t$qbase\t$temp\t$rpos\n";
         }
         else{$snp_location{$rpos}{"$reference:$query"}= $qbase;}
         $positions{$rpos}{"$ref:$query"}='read';
      }
   }
}
close IN;
}

sub check_positions
{
my $prev=0;
my $count=0;
my $exclude=0;
foreach my $position (sort {$a<=>$b}keys %positions){
   my $difference=$position-$prev;
   if ($difference<10){
      $exclude++;
      $gap_location{$position}++;
      $prev= $position;
   }
   else {$prev=$position;}
   $count++;
}
print "Excluding $exclude out of $count SNPs\nWarning: This number includes SNPs found in gap regions\n";
}

sub create_snp_array
{
print "Creating SNP alignment.\n";
my $ref=0;
my $current_snp=0;
my $first=0;
my $second=0;

foreach my $comparison(@header_list){
   if ($comparison=~/(.+):(.+)/){
      ($first,$second)=($1,$2);
      if ($skip_query_ref->{$second}) {print $second,"\n"; next;}
      print OUT ">$2\n"; 
      if ($genbank==1){print CDSOUT ">$2\n";}
   }
   #print $second,"\n";
   foreach (sort {$a<=>$b}keys %snp_location){
      $current_snp=$_ -1;
      if (!defined $gap_location{$_}){
         if (defined $snp_location{$_}{$comparison}){
            print OUT $snp_location{$_}{$comparison};
            if ($genbank==1){
               if (defined $coding_location{$_}){
                  print CDSOUT $snp_location{$_}{$comparison};
               }
            }
            if ($first ne $second){
               $ref= substr($ref_sequence,$current_snp,1);
#               print "$current_snp\n";
               if ($genbank==1){
                  if (!defined $coding_location{$_}){
                     print STAT "$first\t$second\tnoncoding SNP\t$_\t$positions{$_}{$comparison}\t$ref\t$snp_location{$_}{$comparison}\n";
                  }
                  elsif (defined $coding_location{$_}){
                     my $snp=$_;
                     my ($start,$end,$product)= split /,/, $coding_location{$_};
                     print STAT "$first\t$second\tcoding SNP\t$snp\t$positions{$snp}{$comparison}\t$ref\t$snp_location{$snp}{$comparison}\t$start\t$end\n";
                     print CDSSTAT "$first\t$second\t$snp\t$positions{$snp}{$comparison}\t$ref\t$snp_location{$snp}{$comparison}\t$start\t$end\n";
                     
                  }
               }
               else{print STAT "$first\t$second\t$_\t$positions{$_}{$comparison}\t$ref\t$snp_location{$_}{$comparison}\n";}
            }
         }
         elsif(!defined $snp_location{$_}{$comparison}){
            $ref= substr($ref_sequence,$current_snp,1);
#            print "$current_snp\n";
            print OUT $ref;
            if ($genbank==1){if (defined $coding_location{$_}){print CDSOUT $ref;}}
         }
      }
   }
   print OUT "\n";
   if ($genbank==1){print CDSOUT "\n";}
}
}

sub create_comparison
{
my $ref=0;
my $current=0;

print "Creating SNP stats file. \n";
foreach (@header_list){#   print ">$comparison\n";
   if (/(.+):(.+)/){
      #print COMP "\t$2";
      if ($genbank==1){print CDSSTAT "\t$2"}  
   }
}
#print COMP "\n";
if ($genbank==1){print CDSSTAT "\n"}

foreach (sort {$a<=>$b}keys %snp_location){
   $current=$_ -1;
   if (!defined $gap_location{$_}){
 #     print COMP "$_\t";
      foreach my $comparison(@header_list){
         if (defined $snp_location{$_}{$comparison}){print COMP "$snp_location{$_}{$comparison}\t";}
         elsif(!defined $snp_location{$_}{$comparison}){
            $ref= substr($ref_sequence,$current,1);
            print COMP "$ref\t";
         }
      }
   #   print COMP "\n";
   }
   if ($genbank==1){
      if (defined $coding_location{$_}){
         print CDSSTAT "$_\t";
         foreach my $comparison(@header_list){
            if (defined $snp_location{$_}{$comparison}){print CDSSTAT "$snp_location{$_}{$comparison}\t";}
            elsif(!defined $snp_location{$_}{$comparison}){
               $ref= substr($ref_sequence,$current,1);
               print CDSSTAT "$ref\t";
            }
         }
         print CDSSTAT "\n";
      }
   }
}
}

close OUT;
close STAT;
close AMB;
close COMP;
close CDSSTAT;
 
sub usage
{
print STDERR"
Usage:
   $0 -i <in directory> -r <reference filename> -l <query list file> -o <output directory>

   -i          in directory name
   -r          reference file
   -l 	       header list
   -g          (1|0) with genbank CDS coords
   -help|h|?   print usage/help
";
exit;
}

