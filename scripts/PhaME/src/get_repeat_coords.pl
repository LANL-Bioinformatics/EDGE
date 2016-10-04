#!/usr/bin/perl

use strict;
use Getopt::Long;
use FindBin qw($RealBin);

my $identity=95;
my $len_cutoff=0;
my $output="repeats_coords.txt";
my $stats= "repeats_stats.txt";

# set up environments
$ENV{PATH}="$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

GetOptions(
   'i=i'      => \$identity,
   'l=i'      => \$len_cutoff,
   'o=s'      => \$output,
   's=s'      => \$stats,
   'help|?'   => sub{Usage()},
);

sub Usage
{
   print <<USAGE;
   perl $0 [options] <fasta>
        -i INT      the identity cutoff 0 to 100 (default: 95)
        -l INT      the repeat length cutoff (default:0)
        -o STRING   output filename (default: repeats_coords.txt)
        -s STRING   output stats filename (default: repeats_stats.txt)

USAGE
exit;
}

my $file=$ARGV[0];
&Usage unless (-e $file);

my $command="nucmer --maxmatch --nosimplify --prefix=seq_seq$$ $file $file 2>/dev/null";
print "Running self-nucmer on $file.\n";
if (system ("$command")) {die "$command"}; 
# apply identity cutoff and lenght cutoff and use awk to skip self-hits
my $command="show-coords -r -I $identity -L $len_cutoff -T seq_seq$$.delta | awk \'\$1 != \$3 || \$2 != \$4 && \$8==\$9 {print}\' > seq_seq$$.coords";
if (system ("$command")) {die "$command"};
&get_coords_file("seq_seq$$.coords");
unlink "seq_seq$$.delta";
unlink "seq_seq$$.coords";

sub get_coords_file 
{
my $file=shift;
my %hash;
open (IN,$file);
my $comparison_files=<IN>;
my $tmp=<IN>;
my $field_head=<IN>;
my @comparison_files=split /\s+/,$comparison_files;
my $seq_len = &get_fasta_len($comparison_files[0]);
my %seq_len=%{$seq_len};
my $seq_num=0;
my $total_seq_len=0;
my $total_repeats_len=0;
my $repeats_seg_num=0;
my $repeat_id="r00000";
open (OUT, ">$output");
open (OUT1, ">$stats");
# this is nucmer coordinate output
#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[TAGS]
#15376	16742	607219	608585	1367	1367	99.56	gi|49175990|ref|NC_000913.2|	gi|49175990|ref|NC_000913.2|
# read through each repeat region and record its position of the seq in hash 

while(<IN>){
   chomp;
   my @fields=split /\t/,$_;
   my $seq_id=$fields[7];
   my $start=$fields[0];
   my $end=$fields[1];
   for my $pos ($start..$end){
      $hash{$seq_id}->{$pos}=1;
   }
}
close IN;

foreach my $seq_id (sort keys %seq_len){  
   $seq_num++;
   my @repeats;
   my $seq_total_repeats_len=0;
   my $seq_repeats_seg_num=0;
   my $len = $seq_len{$seq_id};
   for my $pos (1..$len){
      if($hash{$seq_id}->{$pos}==1){ 
         $seq_total_repeats_len++;
         push @repeats, $pos;
      }
      else{
         if (@repeats){
            $seq_repeats_seg_num++;
            $repeat_id++;
            my $r_start= $repeats[0];
            my $r_end= $repeats[-1];
            print OUT $seq_id, "\t", $repeat_id, "\t",$r_start,"\t",$r_end, "\t",$r_end-$r_start+1,"\n";
            undef @repeats;
         }
      }
   }
   if (@repeats){
      $seq_repeats_seg_num++;
      $repeat_id++;
      my $r_start= $repeats[0];
      my $r_end= $repeats[-1];
      print OUT $seq_id, "\t", $repeat_id, "\t",$r_start,"\t",$r_end, "\t",$r_end-$r_start+1,"\n";
      undef @repeats;
   }
   print OUT1 "\n$seq_id size:\t$len\n";
   print OUT1 "Repeats segment #:\t$seq_repeats_seg_num\n";
   printf OUT1 ("Repeats total length:\t%d (%.2f%%)\n", $seq_total_repeats_len, $seq_total_repeats_len/$len*100);
   $total_repeats_len += $seq_total_repeats_len;
   $total_seq_len += $len;
}
close OUT;

if ($seq_num>1){
   print OUT1 "\nOVERALL \nTotal seq #:\t$seq_num\n";
   print OUT1 " Total bases #:\t$total_seq_len\n";
   printf OUT1 (" Repeats total length:\t%d (%.2f%%)\n",$total_repeats_len,$total_repeats_len/$total_seq_len*100);
}
print OUT1;
}

sub get_fasta_len 
{
my $file=shift;
my %len;
my $id;
open (INFASTA, $file);
while(<INFASTA>){
   chomp;
   if (/^>(\S+)/){$id = $1;}
   else{
      $_=~ s/\n//g;
      $_=~ s/ //g;
      $len{$id} += length $_;
   }
}
close INFASTA;
return \%len;
}

