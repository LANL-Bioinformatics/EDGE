#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;

my $len_cutoff="0";
GetOptions(
   "l|len_cutoff=i" => \$len_cutoff
);
my $file=$ARGV[0];
if (scalar(@ARGV)<1){&usage();}

my %seq_hash;
my %query_hash;
my $reference;
my $query;
my ($start_r,$end_r,$start_q,$end_q);
my $repeats=0;
my $repeat=0;

open (IN,"$file") or die "Can't open $file:$!";
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

while(<IN>){
   if ($_=~ /^\d+/){
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
      if ($ref_tag=~/.+_(\d+)_(\d+)$/){
         ($start_r,$end_r)=($1,$2);
         $r_loopstart=$start_r+$r_loopstart-1;
         $r_loopend=$r_loopstart+$r_hit_len;
         $repeat=1;
      }
      if ($query_tag=~/.+_(\d+)_(\d+)$/){
         ($start_q,$end_q)=($1,$2);
         $q_loopstart=$start_q+$q_loopstart-1;
         $q_loopend=$q_loopstart+$q_hit_len;
      }
      if ($r_start>$r_end){
         $r_loopstart=$r_end+$start_r-1;
         $r_loopend=$r_loopstart+$r_hit_len;
      }
      for (my $i=$r_loopstart;$i<=$r_loopend;$i++){$seq_hash{$ref_tag}->{$i}++;}
      if ($q_start>$q_end){
         $q_loopstart=$start_q+$q_end-1;
         $q_loopend=$q_loopstart+$q_hit_len;
      }
      for (my $i=$q_loopstart; $i<= $q_loopend; $i++){$query_hash{$query_tag}->{$i}++;}
   }
}
close IN;

open (OUT,"> ${r_file_name}_${q_file_name}.gaps");
open (OUT2,">${q_file_name}_${r_file_name}.gaps");

#foreach my $ref (sort {$seq_hash{$b}->{len}<=>$seq_hash{$a}->{len}} keys %seq_hash)
my $ref_total_gaps_bases=0;
my $ref_total_gaps_number=0;
my $gap_num=0;
foreach my $ref (sort keys %seq_hash){
   my ($genome_covered_pos,$gaps,@gap_array,$cov,$total_cov,$gap_len);
   my $ref_len_1=$seq_hash{$ref}->{len};
   if ($repeat==1){
      if ($ref=~/(.+)_(\d+)_(\d+)$/){($reference,$start_r,$end_r)=($1,$2,$3);}
   }else{ # $repeat==0
      $reference=$ref;
      $start_r=1;
      $end_r=$ref_len_1;
   }
   for ($start_r..$end_r){
       if ($seq_hash{$ref}->{$_}){
          if (@gap_array){
             $gap_len=$gap_array[-1] - $gap_array[0]+1;
             if ($gap_len > $len_cutoff){
                $gap_num++;
                $ref_total_gaps_bases += $gap_len;
#                print $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
                print OUT $reference,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
             }
             @gap_array=();
          }
          $genome_covered_pos++;
          $cov=$seq_hash{$ref}->{$_};
          $total_cov += $cov;
       }
       else{
          push @gap_array, $_;
          $gaps++;
       }
   }
   
   if (@gap_array){
      $gap_len=$gap_array[-1] - $gap_array[0]+1;
      if ($gap_len > $len_cutoff){
         $gap_num++;
         $ref_total_gaps_bases += $gap_len;
#         print $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
         print OUT $reference,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$q_file_name,"\n";
      }
   }
   $ref_total_gaps_number=+$gap_num;
}
print $ref_total_gaps_number," (",$ref_total_gaps_bases,")","\n";

my $query_total_gaps_bases=0;
my $query_total_gaps_number=0;
$gap_num=0;
foreach my $ref (sort keys %query_hash){
   my ($genome_covered_pos,$gaps,@gap_array,$cov,$total_cov,$gap_len);
   my $ref_len_1=$query_hash{$ref}->{len};
   if ($q_file_name =~ /contig/){
       $query = $ref;
       $start_r = 1 ;
       $end_r = $ref_len_1;
   }else{
       if ($ref=~/(.+)_(\d+)_(\d+)$/){($query,$start_r,$end_r)=($1,$2,$3);}
   }
   for ($start_r..$end_r){
      if ($query_hash{$ref}->{$_}){
         if (@gap_array){
            $gap_len = $gap_array[-1] - $gap_array[0]+1;
            if ($gap_len > $len_cutoff){
               $gap_num++;
               $query_total_gaps_bases += $gap_len;
#               print $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
               print OUT2 $query,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
            }
            @gap_array=();
         }
         $genome_covered_pos++;
         $cov=$query_hash{$ref}->{$_};
         $total_cov += $cov;
      }
      else{
         push @gap_array, $_;
         $gaps++;
      }
   }
   if (@gap_array){
      $gap_len = $gap_array[-1] - $gap_array[0]+1;
      if ($gap_len > $len_cutoff){
         $gap_num++;
         $query_total_gaps_bases += $gap_len;
#         print $ref,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
         print OUT2 $query,"\t",$gap_array[0],"\t",$gap_array[-1],"\t",$gap_len,"\t",$r_file_name,"\n";
      }
   }
   $query_total_gaps_number=+$gap_num;
   printf STDERR ("%s\t%d\t%d\n",$ref, $gap_num, $gaps);
}
print $query_total_gaps_number," (",$query_total_gaps_bases,")","\n";

close OUT;
close OUT2;

sub fastalength 
{
my ($file)= @_;
open (IN1,$file);
my ($id,$seq,%hash);
while (<IN1>){
   chomp;
   if (/>(\S+)/){
      if ($seq){$hash{$id}->{len}=length($seq);}
      $id=$1;
      $seq="";
   }
   else{$seq .= $_;}
}
if ($seq){$hash{$id}->{len}=length($seq);}
close IN1;
#print "$id\t$seq\n";
return (%hash);
}

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

sub usage
{
print <<Usage;
$0 -q <query_fasta> -d <output_directory> -t <# threads>

Usage
exit;
}

