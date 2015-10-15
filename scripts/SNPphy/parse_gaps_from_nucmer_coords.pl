#! /usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;

$ENV{PATH}  = "/sardis1/chienchi/bin:/usr/local/bin:/usr/bin:/gorgonzola2/bin:/bin:/opt/apps/bin";

my $len_cutoff = "0";
GetOptions( "l|len_cutoff=i"   => \$len_cutoff);

my $file1=$ARGV[0]; #  output from nucmer | show-coords -clrT 

if (scalar(@ARGV)<1){              
   print STDERR "Usage:\nperl $0 [-len_cutoff #] nucmer_show_coords_output(show-coords -clrT)\n\n" ;
   exit;
}


my %seq_hash;
my %query_hash;
open (IN,"$file1") or die "Can't open $file1:$!";
my $ref_query=<IN>;
my ($ref_full_path,$query_full_path)=split /\s+/, $ref_query;
(my $r_file_name, undef, undef)=fileparse("$ref_full_path", qr/\.[^.]*/);
(my $q_file_name, undef, undef)=fileparse("$query_full_path", qr/\.[^.]*/);
$r_file_name =~ s/\.fna//;
$q_file_name =~ s/\.fna//;
$r_file_name =~ s/\.fasta//;
$q_file_name =~ s/\.fasta//;

%seq_hash=&fastalength($ref_full_path);
%query_hash=&fastalength($query_full_path);

while (<IN>)
{
    if ($_=~ /^\d+/)
    {
       my @array=split;
       my $r_start=$array[0];
       my $r_end=$array[1];
       my $q_start=$array[2];
       my $q_end=$array[3];
       my $r_hit_len=$array[4];
       my $q_hit_len=$array[5];
       my $idy=$array[6];
       my $ref_len=$array[7];
       my $query_len=$array[8];
       my $ref_coverage=$array[9];
       my $query_coverage=$array[10];
       my $ref_tag=$array[11];
       my $query_tag=$array[12];
       my $r_loopstart=$r_start;
       my $r_loopend=$r_end;
       my $q_loopstart=$q_start;
       my $q_loopend=$q_end;
       #$seq_hash{$ref_tag}->{len}=$ref_len;
       #$query_hash{$query_tag}->{len}=$query_len;
       if ($r_start>$r_end)
       {
           $r_loopstart=$r_end;
           $r_loopend=$r_start;
       }
       for (my $i=$r_loopstart; $i<= $r_loopend; $i++)
       {
          $seq_hash{$ref_tag}->{$i}++; 
       }
       if ($q_start>$q_end)
       {
           $q_loopstart=$q_end;
           $q_loopend=$q_start;
       }
       for (my $i=$q_loopstart; $i<= $q_loopend; $i++)
       {
          $query_hash{$query_tag}->{$i}++; 
       }
    }
    
}
close IN;


open (OUT,"> ${r_file_name}_${q_file_name}.gaps_coord.txt");
open (OUT2,">${q_file_name}_${r_file_name}.gaps_coord.txt");

#foreach my $ref (sort {$seq_hash{$b}->{len}<=>$seq_hash{$a}->{len}} keys %seq_hash)
my $ref_total_gaps_bases=0;
my $ref_total_gaps_number=0;
foreach my $ref (sort keys %seq_hash)
{
  my ($genome_covered_pos, $gaps, @gap_array, $gap_num,$cov,$total_cov, $gap_len);
  my $ref_len_1=$seq_hash{$ref}->{len};
  #print $ref," ",$ref_len_1,"\n";
  for (1..$ref_len_1){
    if ($seq_hash{$ref}->{$_})
    {
        if (@gap_array){
         $gap_len = $gap_array[-1] - $gap_array[0]+1;
         if ($gap_len > $len_cutoff){
             $gap_num++; 
             
             $ref_total_gaps_bases += $gap_len;
             print OUT $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
         }
         @gap_array=();
        }
        $genome_covered_pos++;
        $cov=$seq_hash{$ref}->{$_};
        $total_cov += $cov;
    }
    else
    {
       push @gap_array, $_;
       $gaps++;
    }
  }

  if (@gap_array){
      $gap_len = $gap_array[-1] - $gap_array[0]+1;
      if ($gap_len > $len_cutoff){
         $gap_num++;
         $ref_total_gaps_bases += $gap_len;
         print OUT $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
      }
  }
  
  #print OUT $ref,"\t",$seq_hash{$ref}->{len},"\t",$seq_hash{$ref}->{GC},"\t",
  #      $total_cov/$ref_len_1,"\t",
  #      $genome_covered_pos/$ref_len_1*100,"\n";
  #$total_ref_len += $ref_len_1;
  #$total_ref_cov += $genome_covered_pos;
  #$total_fold_cov += $total_cov; 
  $ref_total_gaps_number += $gap_num;
  #printf STDERR ("%s\t%d\t%d\n",$ref, $gap_num, $gaps);
} 
   print $ref_total_gaps_number," (",$ref_total_gaps_bases,")","\n";
   #my $avg_cov_fold=$total_fold_cov/$total_ref_len;
   #my $coverage=$total_ref_cov/$total_ref_len * 100;
   #printf STDERR ("Coverage:\t%.4f%%\n",$coverage);

my $query_total_gaps_bases=0;
my $query_total_gaps_number=0;
foreach my $ref (sort keys %query_hash)
{
  my ($genome_covered_pos, $gaps, @gap_array, $gap_num,$cov,$total_cov,$gap_len);
  my $ref_len_1=$query_hash{$ref}->{len};
  for (1..$ref_len_1){
    if ($query_hash{$ref}->{$_})
    {
        if (@gap_array){
         $gap_len = $gap_array[-1] - $gap_array[0]+1;
         if ($gap_len > $len_cutoff){
             $gap_num++; 
             $query_total_gaps_bases += $gap_len;
             print OUT2 $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
         }
         @gap_array=();
        }
        $genome_covered_pos++;
        $cov=$query_hash{$ref}->{$_};
        $total_cov += $cov;
    }
    else
    {
       push @gap_array, $_;
       $gaps++;
    }
  }

  if (@gap_array){
        $gap_len = $gap_array[-1] - $gap_array[0]+1;
        if ($gap_len > $len_cutoff){
         $gap_num++;
         $query_total_gaps_bases += $gap_len;
         print OUT2 $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
      }
  }
  
  $query_total_gaps_number += $gap_num;
  printf STDERR ("%s\t%d\t%d\n",$ref, $gap_num, $gaps);
}
  print $query_total_gaps_number," (",$query_total_gaps_bases,")","\n";
   
close OUT;
close OUT2;

sub fastalength {
    my ($file)= @_;
    open (IN1, $file);
    my ($id, $seq, %hash);
    while (<IN1>)
    {
       chomp;
       if (/>(\S+)/)
       {   
           if ($seq)
           { 
              $hash{$id}->{len}=length($seq);
           }
           $id = $1; 
           $seq="";
       }
       else
       {
           $seq .= $_;
       }
    }
    if ($seq)
    { 
       $hash{$id}->{len}=length($seq);
    }
    close IN1;
    return (%hash); 
}




