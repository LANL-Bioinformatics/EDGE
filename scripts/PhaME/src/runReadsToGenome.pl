#! /usr/bin/perl
# required: 1. R
#           2. samtools > 1.1 
#           3. bwa 0.6 sampe patched  
#           4. bowtie2
#           5. bcftools  > 1.1 
#           6. vcfutils.pl  (from samtools package)
#           7. snap
#     input: paired reads files: forward.fasta/q and reverse.fasta/q
#            reference genome
#     output: bam file (reads placement from bwa + samtools)
#             aln_stats.txt
#             coverage plots: genome plot and histogram
#             gap coordiates
#             SNP file in variant call format(VCF v4.1)
# chienchi@lanl.gov
# 20100811
# 20110125 updated for samtools and bwa
# 20110617 window size coverage plot
# 20120112 add -aligner
# 20120327 add proper and unproper paired comparision plot and -plot_only flag

use Getopt::Long;
use File::Basename;
use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../lib";
use fastq_utility;

# set up environments
$ENV{PATH}="$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

my $debug=0;

$|=1;
my ($file1, $file2, $paired_files,$prefix, $ref_file, $outDir,$file_long,$singleton,$pacbio, $offset);
my $bwa_options="-t 4 ";
my $bowtie_options="-p 8 -a";
my $snap_options="-t 4 -M ";
my $minimap2_options="-t 4 ";
my $aligner="bwa";
my ($window_size, $step_size)=(1000,200);
my $pacbio_bwa_option="-b5 -q2 -r1 -z10 "; 
my $prefix="ReadsMapping";

GetOptions( 
   'aligner=s'        => \$aligner,
   'p=s'              => \$paired_files,
   'ref=s'            => \$ref_file, # reference/contigs file
   'pre=s'            => \$prefix,
   'long=s'           =>  \$file_long,
   'u=s'              => \$singleton, # illumina singleton 
   'window_size=i'    => \$window_size,  # for coverage plot
   'step_size=i'      => \$step_size,  # for coverage plot
   'd=s'              => \$outDir,
   'bwa_options=s'    => \$bwa_options,
   'bowtie_options=s' => \$bowtie_options,
   'minimap2_options=s'  => \$minimap2_options,
   'snap_options=s'   => \$snap_options,
   'pacbio'           => \$pacbio,
   'debug'            => \$debug,
   'help|?',          sub {Usage()}
);

my $tmp = "$outDir";

## input check ##
unless ( -e $ref_file && $outDir) { &Usage;}
unless ( $paired_files or -e $file_long or -e $singleton) { &Usage; }
if ($paired_files){
  ($file1, $file2) = split /,/,$paired_files;
  unless (-e $file1 && -e $file2) {print "$file1 or $file2 not exists\n";&Usage;}
}

if ($step_size > $window_size) {die "The step_size ($step_size) should be less than window_size ($window_size)\n&Usage";}

## output file variable initialized ##
my $stats_output="$outDir/$prefix.alnstats.txt";
my $vcf_output="$outDir/$prefix.vcf";
my $bcf_output="$outDir/$prefix.raw.bcf";
my $bam_output="$outDir/$prefix.sort.bam";
my $bam_index_output="$outDir/$prefix.sort.bam.bai";
my $pileup_output="$outDir/$prefix.pileup";
my $ref_window_gc="$outDir/$prefix.ref_windows_gc.txt";

if (! -e $outDir){mkdir $outDir;}

my ($bwa_threads)= $bwa_options =~ /-t (\d+)/;
my ($bowtie_threads)= $bowtie_options =~ /-p (\d+)/;
my ($snap_threads)= $snap_options =~ /-t (\d+)/;
my ($minimap2_threads)= $minimap2_options =~ /-t (\d+)/;
my $samtools_threads = $bwa_threads || $bowtie_threads || $snap_threads || $minimap2_threads || "1";
my ($ref_file_name, $ref_file_path, $ref_file_suffix)=fileparse("$ref_file", qr/\.[^.]*/);

# index reference
if ( $aligner =~ /bowtie/i and ! -e "$ref_file.1.bt2"){
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
   $ref_file=&fold($ref_file);
   `bowtie2-build $ref_file $ref_file`;
}
elsif ($aligner =~ /bwa/i and ! -e "$ref_file.bwt"){
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
   $ref_file=&fold($ref_file);
   `bwa index $ref_file`;
}
elsif ($aligner =~ /snap/i){
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
   $ref_file=&fold($ref_file);
   `snap index $ref_file $ref_file.snap `;
}
elsif ($aligner =~ /minimap2/i ){
     # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    $ref_file=&fold($ref_file);
    `minimap2 -d $tmp/$ref_file_name.mmi $ref_file `;
}

## index reference sequence 
`samtools faidx $ref_file`; 

if ($file_long){
   print "Mapping long reads\n";
   if ($aligner =~ /bowtie/i){`bowtie2 -a --local $bowtie_options -x $ref_file -fU $file_long -S $outDir/LongReads$$.sam`;}
   elsif($aligner =~ /bwa/i){
     if ($pacbio){`bwa bwasw -M -H $pacbio_bwa_option -t $bwa_threads $ref_file $file_long -f $outDir/LongReads$$.sam`;}
     else{`bwa bwasw -M -H -t $bwa_threads $ref_file $file_long -f $outDir/LongReads$$.sam`;}
  #my $mapped_Long_reads=`awk '\$3 !~/*/ && \$1 !~/\@SQ/ {print \$1}' /tmp/LongReads$$.sam | uniq - | wc -l`;
  #`echo -e "Mapped_reads_number:\t$mapped_Long_reads" >>$outDir/LongReads_aln_stats.txt`;
   }
   elsif ($aligner =~ /snap/i){`snap single $ref_file.snap $file_long -o $outDir/LongReads$$.sam $snap_options`;}
   elsif ($aligner =~ /minimap2/i){
          `minimap2 -La $minimap2_options  $tmp/$ref_file_name.mmi $file_long > $outDir/LongReads$$.sam`;
   }
   `samtools view -t $ref_file.fai -@ $samtools_threads -uhS $outDir/LongReads$$.sam | samtools sort -@ $samtools_threads -O BAM -T $outDir -o $outDir/LongReads$$.bam - `;
}
if ($paired_files){
   print "Mapping paired end reads\n";
   $offset=fastq_utility::checkQualityFormat($file1);
   if ($offset==64){$bowtie_options=$bowtie_options." --phred64 ";$bwa_options=$bwa_options." -I ";}
   if ($aligner=~/bowtie/i){`bowtie2 $bowtie_options -x $ref_file -1 $file1 -2 $file2 -S $outDir/paired$$.sam`;}
   elsif ($aligner =~ /bwa/i){
      `bwa aln $bwa_options $ref_file $file1 > /tmp/reads_1_$$.sai`;
      `bwa aln $bwa_options $ref_file $file2 > /tmp/reads_2_$$.sai`;
      `bwa sampe -t $bwa_threads -a 100000 $ref_file /tmp/reads_1_$$.sai /tmp/reads_2_$$.sai $file1 $file2 > $outDir/paired$$.sam`;
   }
   elsif ($aligner =~ /snap/i){`snap paired $ref_file.snap $file1 $file2 -o $outDir/paired$$.sam $snap_options`;}
   elsif ($aligner =~ /minimap2/i){
      `minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $file1 $file2 > $outDir/paired$$.sam`;
   }
   `samtools view -t $ref_file.fai -@ $samtools_threads -uhS $outDir/paired$$.sam | samtools sort -@ $samtools_threads -O BAM -T $outDir -o $outDir/paired$$.bam -`;
}

if ($singleton){
   print "Mapping single end reads\n";
   $offset=fastq_utility::checkQualityFormat($singleton);
   if ($offset==64){$bowtie_options=$bowtie_options."--phred64 ";$bwa_options=$bwa_options."-I ";}
   if ($aligner =~ /bowtie/i){`bowtie2 $bowtie_options -x $ref_file -U $singleton -S $outDir/singleton$$.sam`;}
   elsif($aligner =~ /bwa/i){
      `bwa aln $bwa_options $ref_file $singleton > /tmp/singleton$$.sai`;
      `bwa samse -n 50 -t $bwa_threads $ref_file /tmp/singleton$$.sai $singleton > $outDir/singleton$$.sam`;
    }
    elsif($aligner =~ /snap/i){`snap single $ref_file.snap $file_long -o $outDir/singleton$$.sam $snap_options`;}
    elsif ($aligner =~ /minimap2/i){
        `minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $singleton> $outDir/singleton$$.sam`;
    }
    `samtools view -t $ref_file.fai -@ $samtools_threads -uhS $outDir/singleton$$.sam | samtools sort -@ $samtools_threads -O BAM -T $outDir -o $outDir/singleton$$.bam -`;
}

# merge bam files if there are different file type, paired, single end, long..
if ($file_long and $paired_files and $singleton){
   `samtools merge -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam $outDir/LongReads$$.bam`;
}
elsif($file_long and $paired_files){
   `samtools merge -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/LongReads$$.bam`;
}
elsif($paired_files and $singleton){
   `samtools merge -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam`;
}
elsif($singleton and $file_long){
   `samtools merge -@ $samtools_threads $bam_output $outDir/singleton$$.bam $outDir/LongReads$$.bam`;
}
elsif($paired_files){
   `mv $outDir/paired$$.bam $bam_output`;
}
elsif($singleton){
   `mv $outDir/singleton$$.bam $bam_output`;
}
elsif($file_long){
   `mv $outDir/LongReads$$.bam $bam_output`;
}


## index BAM file 
`samtools index $bam_output $bam_index_output`; 

## generate statistical numbers 
print "Generate alignment statistical numbers \n";
`samtools flagstat $bam_output > $stats_output`; 
 
## SNP call
print "SNPs/Indels call...\n";
my $max_depth=100000;
my $min_indel_candidate_depth=3;
my $min_alt_bases=3;
my $min_depth=1;

`bcftools mpileup -d $max_depth -L $max_depth -m $min_indel_candidate_depth -Ov -f $ref_file $bam_output | bcftools call -cO b - > $bcf_output 2>/dev/null`;
`bcftools view -v snps,indels,mnps,ref,bnd,other -Ov $bcf_output | vcfutils.pl varFilter -a$min_alt_bases -d$min_depth -D$max_depth > $vcf_output`;

## derived chimera info 
if ($aligner=~ /bwa/i and $paired_files){ 
   my $proper_paired = `grep "properly paired"  $stats_output | awk '{print \$1}' `;
   my $all_mapped_paired = `grep "with itself and mate mapped"  $stats_output | awk '{print \$1}' `;
   my $chimera = $all_mapped_paired - $proper_paired;
   chomp $chimera;
   `echo -e "Chimera:\t$chimera" >>$stats_output`;
}

## generate genome coverage plots and histograms 
print "Generate genome coverage plots and histograms...\n";
`samtools mpileup -ABQ0 -d10000000 -f  $ref_file $bam_output >$pileup_output`; 

if ($paired_files){
  ## generate proper-paired reads coverage
   `samtools view -@ $samtools_threads -u -h -f 2 $bam_output | samtools mpileup -BQ0 -d10000000 -f $ref_file - | awk '{print \$1"\\t"\$2"\\t"\$4}'  > $outDir/proper_paired$$.coverage`;

  ## generate non-proper-paired reads coverage 2 (properpaired)+4(query unmapped)+8(mate unmapped)
   `samtools view -@ $samtools_threads -u -h -F 14 $bam_output | samtools mpileup -ABQ0 -d10000000 -f $ref_file - | awk '{print \$1"\\t"\$2"\\t"\$4}'  > $outDir/unproper_paired$$.coverage`;
}

# build base coverage hash
open (IN,"$pileup_output");
my %base_hash;
while (<IN>){
   chomp;
   my ($id,$pos,$ref_base,$cov,$seq,$qual)=split /\t/;
   $base_hash{$id}->{$pos}=$cov;
}
close IN;

my %proper_base_hash;
my %unproper_base_hash;
if ($paired_files){
  # build proper_paired mapped reads base coverage hash
   open (IN,"$outDir/proper_paired$$.coverage");
   while (<IN>){
      chomp;
      my ($id,$pos,$cov)=split /\t/ ;
      $proper_base_hash{$id}->{$pos}=$cov;
   }
   close IN;

# build unproper_paired mapped reads base coverage hash
   open (IN,"$outDir/unproper_paired$$.coverage");
   while (<IN>){
      chomp;
      my ($id,$pos,$cov)=split /\t/;
      $unproper_base_hash{$id}->{$pos}=$cov;
   }
   close IN;
}

# get reference informaiton
my $ref_hash=&get_ref_info($ref_file);
$ref_hash=&mapped_reads_per_contigs($bam_output,$ref_hash);

my $stats_print_string="\nRef\tRef_len\tRef_GC%\tMapped_reads\tRef_recovery%\tAvg_fold(x)\tFold_std\tNum_of_Gap\tTotal_Gap_bases";
$stats_print_string .= "\tNum_of_SNPs\tNum_of_INDELs";
$stats_print_string .="\n";
foreach my $ref_name (sort keys %{$ref_hash}){      
   my ($snp_num , $indel_num);
   my $ref_len = $ref_hash->{$ref_name}->{len};
   my $ref_GC = $ref_hash->{$ref_name}->{GC};
   my $mapped_reads = $ref_hash->{$ref_name}->{reads};
   $stats_print_string .= $ref_name."\t".$ref_len."\t".$ref_GC."\t".$mapped_reads."\t";

   # generate coverage file
   my $coverage_output="$outDir/${prefix}_${ref_name}.coverage";
   my $WindowCoverage_output="$outDir/${prefix}_${ref_name}.window_size_coverage";
   my $gap_output="$outDir/${prefix}_${ref_name}.gaps";
   my $coverage_plot="$outDir/${prefix}_${ref_name}_base_coverage.png";
   my $histogram="$outDir/${prefix}_${ref_name}_coverage_histogram.png";
   $stats_print_string .= &window_size_coverage($coverage_output,$WindowCoverage_output,\%base_hash,$gap_output,$ref_name,$ref_len);
  
   my $properpair_coverage_output;
   my $unproperpair_coverage_output;
   my $other_coverage_plot;
   if ($paired_files){
      $properpair_coverage_output="$outDir/${prefix}_${ref_name}.p$$.window_size_coverage";
      $unproperpair_coverage_output="$outDir/${prefix}_${ref_name}.up$$.window_size_coverage";
      $other_coverage_plot="$outDir/${prefix}_${ref_name}_coverage_comparison.png";
      &window_size_coverage("",$properpair_coverage_output,\%proper_base_hash,"",$ref_name,$ref_len);
      &window_size_coverage("",$unproperpair_coverage_output,\%unproper_base_hash,"",$ref_name,$ref_len);
   }

   ($snp_num , $indel_num)= &SNP_INDEL_COUNT("$vcf_output","$ref_name");
   $stats_print_string .= $snp_num ."\t". $indel_num;
   `echo "$stats_print_string" >> $stats_output`;

   $stats_print_string="";  
   
   unless ($debug){
      unlink $WindowCoverage_output;
      unlink $unproperpair_coverage_output;
      unlink $properpair_coverage_output;
      unlink $coverage_output;
   }
}

# clean up
`rm -rf /tmp/*$$*`;
unless ($debug){
   unlink $ref_window_gc;
   `rm $outDir/*$$*`;
};

sub mapped_reads_per_contigs 
{
my $bam_output=shift;
my $ref_hash_r=shift;
open (IN, "samtools idxstats $bam_output |") or die "$!\n";
while (<IN>){
   chomp;
   my ($id,$len,$mapped,$unmapped)=split /\t/,$_;
   next if($id eq '*');
   $id=~ s/\//_/g;
   $ref_hash_r->{$id}->{reads}=$mapped;
}
close IN;
return $ref_hash_r;
}

sub get_ref_info 
{
# Given reference file
# return hash refernece for id as key and len and GC content.
my $file=$_[0];
my %hash;
my $id;
my $seq;
my $seq_len;
my $GC_content;
my $avg_pos=$window_size/2;
open (OUT,">$ref_window_gc") or die "$!\n";
open (IN,$file) or die "$!\n";
while (<IN>){
   chomp;
   if (/>(\S+)/){
      if ($seq){
         $seq_len = length $seq;
         my $GC_num = $seq=~ tr/GCgc/GCgc/;
         $GC_content = sprintf ("%.2f",$GC_num/$seq_len*100);
         $hash{$id}->{len}= $seq_len;
         $hash{$id}->{GC}=$GC_content;
         for (my $i=0; $i<=$seq_len-$window_size;$i=$i+$step_size){
            my $window_seq=substr($seq,$i,$window_size);  
            $GC_num = $window_seq=~ tr/GCgc/GCgc/;  
            $GC_content = $GC_num/$window_size*100;
            if ($i==0){print OUT $id,"\t",$avg_pos,"\t",$GC_content,"\n";}
            else{print OUT $id,"\t",$avg_pos+$i,"\t",$GC_content,"\n";}
         }
      }
      $id=$1;
      $id =~ s/\//_/g;
      $seq="";
   } 
   else{$seq.=$_;}
}
if ($seq){
   $seq_len = length $seq;
   my $GC_num = $seq=~ tr/GCgc/GCgc/;
   $GC_content = sprintf ("%.2f",$GC_num/$seq_len*100);
   $hash{$id}->{len}= $seq_len;
   $hash{$id}->{GC}=$GC_content;
   for (my $i=0; $i<=$seq_len-$window_size;$i=$i+$step_size){
      my $window_seq=substr($seq,$i,$window_size);  
      $GC_num = $window_seq=~ tr/GCgc/GCgc/;  
      $GC_content = $GC_num/$window_size*100;
      if ($i==0){print OUT $id,"\t",$avg_pos,"\t",$GC_content,"\n";}
      else{print OUT $id,"\t",$avg_pos+$i,"\t",$GC_content,"\n";}
   }
}
close IN;
close OUT;
return \%hash;
}

sub window_size_coverage
{
# given output files names and a hash for each ref name and for each base;
# output coverage per base/window_size per ref. output gap per ref. 
# return statistiacl numbers, genome recovery, fold coveage, fold coverage std, gap number, gap total bases.  
my ($coverage_output,$WindowCoverage_output,$base_hash,$gap_output,$ref_name,$ref_len)=@_;
my $pos_cov;
my $cov_sum;
my $step;
my $window_sum =0;
my $step_sum=0;
my @step_sum;
my $avg_cov=0;
my $avg_pos=$window_size/2;
my @step_sum2;
my $step_sum2;
my $gap_length;
my @gap_array;
my $gap_count=0;
my $gap_total_len=0;
my $covered_base_num;
my @cov_array;
my $stats_return;
if ($coverage_output){
   open (OUT, ">$coverage_output") or die "$! $coverage_output\n";
   open (GAP, ">$gap_output" ) or die "$! $gap_output\n";
   print GAP "Start\tEnd\tLength\tRef_ID\n";
}
open (OUT2, ">$WindowCoverage_output") or die "$! $WindowCoverage_output\n";
for (1..$ref_len){
   if ($base_hash->{$ref_name}->{$_}){
      $pos_cov=$base_hash->{$ref_name}->{$_};
      if ($coverage_output){
         print OUT $_,"\t",$pos_cov,"\n";
         if (@gap_array){
            $gap_length = $gap_array[-1] - $gap_array[0]+1;
            print GAP $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t",$ref_name,"\n";
            $gap_count++;
            $gap_total_len += $gap_length;
            @gap_array=();
         }
         $covered_base_num++;
      }
   }
   else{
      $pos_cov=0;
      if ($coverage_output){
         print OUT $_,"\t",$pos_cov,"\n";
         push @gap_array, $_;
      }
   }
   push @cov_array,$pos_cov;
   $cov_sum += $pos_cov;
   $step_sum += $pos_cov;
   if (($_ % $step_size)==0){
      push @step_sum, $step_sum;
      $step_sum=0;
   }
   if ($_ == $window_size){
      $step=1;
      $window_sum = $cov_sum;
      print OUT2 $avg_pos,"\t",$window_sum/$window_size,"\n";
   }
   if ($_ > $window_size){
      $step_sum2 += $pos_cov;
      if (($_-$window_size)%$step_size==0){
         push @step_sum2, $step_sum2;
         $step_sum2=0;
      }
   }
   if ($_ == ($window_size+$step_size*$step)){
      my $previous_step_sum = shift @step_sum;
      my $after_step_sum = shift @step_sum2;
      $window_sum = $window_sum + $after_step_sum - $previous_step_sum;
      $avg_pos = $avg_pos + $step_size; 
      print OUT2 $avg_pos,"\t",$window_sum/$window_size,"\n";
      $step++;
   }
}
if ($coverage_output){
   if (@gap_array){
      $gap_length = $gap_array[-1] - $gap_array[0]+1;
      print GAP $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t",$ref_name,"\n";
      $gap_count++;
      $gap_total_len += $gap_length;
   }
   my ($std_cov,$avg_cov)= &standard_deviation(@cov_array);
   my $percent_genome_coverage = sprintf ("%.2f%%",$covered_base_num/$ref_len*100);
   my $fold = sprintf ("%.2f",$avg_cov);
   my $fold_std = sprintf ("%.2f",$std_cov);
   $stats_return = $percent_genome_coverage."\t".$fold."\t".$fold_std."\t".$gap_count."\t".$gap_total_len."\t";      
   close OUT;
   return ($stats_return);
}
close OUT2;  
}

sub standard_deviation 
{
my(@numbers) = @_;
# Step 1, find the mean of the numbers
my $total1 = 0;
foreach my $num (@numbers){$total1 += $num;}
my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
my $total2 = 0;
foreach my $num (@numbers){$total2 += ($mean1-$num)**2;}
my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
my $std_dev = sqrt($mean2);
return ($std_dev,$mean1);
}

sub fold 
{
# fold and filter reads length by 200 bp.
my $file=$_[0];
my $seq;
my $seq_name;
my $len_cutoff=0;
my $seq_num;
open (IN,$file);
open (OUT,">/tmp/Contig$$.fold");
while(<IN>){
   chomp;
   if(/>/){
      if ($seq and length($seq)>$len_cutoff){
         $seq =~ s/ //g;
         $seq =~ s/(.{100})/$1\n/g;
         chomp $seq;
         print OUT $seq_name,"\n",$seq,"\n";
      }
      $seq_name=$_;
      $seq="";
      $seq_num++;
   }
   else{$seq.=$_;}
}
if ($seq and length($seq)>$len_cutoff){ # last sequence
   $seq =~ s/ //g;
   $seq =~ s/(.{100})/$1\n/g;
   chomp $seq;
   print OUT $seq_name,"\n",$seq,"\n";
}
close IN;
close OUT;
if ($seq_num<1){die "No seqeucne in your reference file\n";}
return ("/tmp/Contig$$.fold");
}

sub SNP_INDEL_COUNT
{
my  $file=shift;
my  $ref=shift;
open (IN,$file) or die "$!";
my $indel_count=0;
my $SNPs_count=0;
$ref =~ s/\|/\\\|/g;
while(<IN>){
   chomp;
   next if (/^#/);
   next if ($_ !~ /$ref/);
   if (/INDEL/){$indel_count++;}
   else{$SNPs_count++;}
}
close IN;
return ($SNPs_count,$indel_count);
}

sub Usage 
{
print <<"END";
Usage: perl $0 
               -p                        leftSequenceFile,rightSequenceFile 
                                         Comma-separated paired-end reads
               -u                        sequenceFile
                                         Provides a file containing single-end reads.
               -long                     long reads file in fasta format.  
                                         --pacbio   <bool> using this flag combined with -long for Pacbio long reads (bwa only) 
               -ref                      reference sequences file in fasta format
               -pre                      output files' prefix (default "ReadsMapping")
               -d                        output directory
               -aligner                  bwa or bowtie or snap (default: bwa)
               -bwa_options <String>     bwa options
                                         type "bwa aln" to see options
                                         default: "-t 4 "
                                         -t        <int> number of threads [4] 
                                         -I        the input is in the Illumina 1.3+ FASTQ-like format
               -bowtie_options <String>  bowtie options
                                         type "bowtie2 -h" to see options
                                         default: "-p 4 -a "  
                                         -p           <int> number of alignment threads to launch [4] 
                                         -a           report all alignments; very slow
                                         --phred64    qualities are Phred+64
               -snap_options             snap options
                                         type "snap paired" to see options
               -window_size              genome coverage plot (default: 1000 bp)
               -step_size                genome coverage plot (default: 200 bp)
               -skip_aln                 <bool> skip the alignment steps, assume bam files were generated 
                                         and with proper prefix,outpurDir.  default: off
               -no_plot                  <bool> default: off
               -no_snp                   <bool> default: off
               -debug                    <bool> default: off 

Synopsis:
      perl $0 -p reads1.fastq,reads2.fastq -u sinlgeton.fastq -long pyroSeq.fasta -ref reference.fasta -pre ReadsMapping -d /outputPath/

END

exit;
}

