#! /usr/bin/perl
use strict;
use File::Basename;

my $file1=$ARGV[0]; #  output from nucmer | show-snps -CT 

if (scalar(@ARGV)<1){              
   print STDERR "Usage:\nperl $0 nucmer_show_snps_output(show-snps -CT)\n\n" ;
   exit;
}

open (IN, $file1);
my %ref_hash;
my %query_hash;
my $ref_query=<IN>;
my ($ref_full_path,$query_full_path)=split /\s+/, $ref_query;
(my $r_file_name, undef, undef)=fileparse("$ref_full_path", qr/\.[^.]*/);
(my $q_file_name, undef, undef)=fileparse("$query_full_path", qr/\.[^.]*/);
$r_file_name =~ s/\.fna//;
$q_file_name =~ s/\.fna//;
$r_file_name =~ s/\.fasta//;
$q_file_name =~ s/\.fasta//;

%ref_hash=&fastalength($ref_full_path);
%query_hash=&fastalength($query_full_path);

my $SNPs_count=0;
my $Indels_count=0;

#print "\t",$r_file_name,"\t",$q_file_name,"\n";
while (my $line=<IN>){
   chomp $line;
   next if ($line !~ /^\d+/);
   my @array=split /\s+/, $line;
   my $r_position = $array[0];
   my $r_base = $array[1];
   my $q_base = $array[2];
   my $q_position = $array[3];
   my $r_tag = $array[8];
   my $q_tag = $array[9];
    
   if ($r_base =~ /\w/ and $q_base =~ /\w/){  # SNPs
      #print $r_tag,"_",$r_position,"\t",$r_base,"\t",$q_base,"\n";
      $ref_hash{$r_tag}->{SNPpos}=$r_position;
      $query_hash{$q_tag}->{SNPpos}=$q_position;
      $SNPs_count++;
   }
   else{  # Indels 
      #print $r_tag,"_",$r_position,"\t",$r_base,"\t",$q_base,"\n";
      $ref_hash{$r_tag}->{INDELpos}=$r_position;
      $query_hash{$q_tag}->{INDELpos}=$q_position;
      $Indels_count++;
   } 
}
close IN;

print $SNPs_count,"\t",$Indels_count,"\n";


sub fastalength 
{
my ($file)= @_;
open (IN1, $file);
my ($id, $seq, %hash);
while (<IN1>){
   chomp;
   if (/>(\S+)/){   
      if ($seq){$hash{$id}->{len}=length($seq);}
      $id=$1; 
      $seq="";
   }
   else{$seq .=$_;}
}
if ($seq){$hash{$id}->{len}=length($seq);}
close IN1;
return (%hash); 
}

