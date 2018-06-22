#! /usr/bin/perl
# required: 1. R
#           2. samtools > 1.1
#           3. bwa
#           4. bowtie2
#           5. ContigCoverageFold_plots_from_samPileup.pl
#     input: reads files: forward.fastq and reverse.fastq [or single end files and Long fasta files]
#            reference genome/contigs
#     output: bam file (reads placement from bwa/bowtie + samtools)
#             aln_stats.txt
#             Final_contigs.fasta  (contigs with coverage >=80%)
#             contig_coverage.table (five columns:ID,  Length, GC, Avg_fold, Base_Coverage)
#             plots: coverage vs. contig len, avg depth vs. contig len, GC vs. avg depth
#             
# chienchi@lanl.gov
# 20100824
# 20110923 add filter contigs fasta output
# 20120111 change flag names and add bowtie2
# 20180123 use samtools 1.6

use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(cwd chdir);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use fastq_utility;

$|=1;
$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";

my ($file1, $file2, $paired_files, $ref_file, $outDir,$file_long,$singleton,$pacbio,$offset);
my $numCPU = 4;
my $bwa_options="-t $numCPU ";
my $minimap2_options="-t $numCPU ";
my $bowtie_options="-p $numCPU -a ";
my $cov_cut_off=0;
my $aligner="bwa";
my $pacbio_bwa_option="-b5 -q2 -r1 -z10 ";
my $prefix="ReadsMapping";
my ($file1,$file2);
my $skip_aln;
my $tmp;
GetOptions( 
            'aligner=s' => \$aligner,
            'p=s'       => \$paired_files,
            'ref=s' => \$ref_file, # reference/contigs file
            'pre=s' => \$prefix,
            'long=s' =>  \$file_long,
            'u=s' => \$singleton, # illumina singleton 
            'd=s'   => \$outDir,
            'c=f'   => \$cov_cut_off, 
            'cpu=i' => \$numCPU,
            'bwa_options=s' => \$bwa_options,
            'bowtie_options=s' => \$bowtie_options,
            'minimap2_options=s' => \$minimap2_options,
            'skip_aln'  => \$skip_aln,
            'pacbio'  => \$pacbio,
            'help|?',  sub {Usage()}
);

$tmp = $outDir;

unless ( -e $ref_file && $outDir) { &Usage;}
unless ( $paired_files or -e $file_long or -e $singleton) { &Usage; }
if ($paired_files){
  ($file1, $file2) = split /\s+/,$paired_files;
  unless (-e $file1 && -e $file2) {&Usage;}
}


sub Usage 
{
print <<"END";
Usage: perl $0 
               -p                        'leftSequenceFile rightSequenceFile'
                                         Space separated paired-end reads in quote
               -u                        sequenceFile
                                         Provides a file containing single-end reads.
               -long                     long reads file in fasta format.  
                                         --pacbio   <bool> using this flag combined with -long for Pacbio long reads (bwa only)
               -ref                      reference sequences file in fasta format
               -pre                      output files' prefix (default: "ReadsMapping")
               -d                        output directory
               -aligner                  bwa or bowtie (default: bwa)
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
               -minimap2_options         minimap2_options
                                         type "minimap2" to see options
                                         default: "-t 4 "
               -skip_aln                 <bool> skip alignment step
               -cpu <NUM>                number of CPUs [4]  will overwrite the aligner threads
               -c  <NUM>                 cutoff value of contig coverage for final fasta file. (default:noFilter) 

Synopsis:
      perl $0 -p 'reads1.fastq reads2.fastq' -u sinlgeton.fastq -long pyroSeq.fasta -ref contigs.fasta -pre ReadsToContigs -d /outputPath/

END
exit;
}

if (! -e $outDir)
{
     mkdir $outDir;
}

my $stats_output="$outDir/$prefix.alnstats.txt";
my $bam_output="$outDir/$prefix.sort.bam";
my $bam_index_output="$outDir/$prefix.sort.bam.bai";
my $pileup_output="$outDir/$prefix.pileup";

if ($bwa_options =~ /-t\s+\d+/) { $bwa_options =~ s/-t\s+\d+/-t $numCPU/; } else { $bwa_options .= " -t $numCPU ";}
if ($minimap2_options =~ /-t\s+\d+/){$minimap2_options =~ s/-t\s+\d+/-t $numCPU/;}else{$minimap2_options .= " -t $numCPU ";};
if ($bowtie_options =~ /-p\s+\d+/){$bowtie_options =~ s/-p\s+\d+/-p $numCPU/ ;}else{$bowtie_options .= " -p $numCPU ";}
my $samtools_threads = $numCPU;

my ($ref_file_name, $ref_file_path, $ref_file_suffix)=fileparse("$ref_file", qr/\.[^.]*/);

unlink $bam_output;
unlink $bam_index_output;


unless ($skip_aln){
## fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
#$ref_file=&fold($ref_file);
# index reference
if ( $aligner =~ /bowtie/i and ! -e "$ref_file.1.bt2")
{
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    $ref_file=&fold($ref_file);
    `bowtie2-build $ref_file $ref_file`;
}
elsif ($aligner =~ /bwa/i and ! -e "$ref_file.bwt")
{
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    $ref_file=&fold($ref_file);
    `bwa index $ref_file`;

}
elsif ($aligner =~ /minimap2/i ){
    # fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    $ref_file=&fold($ref_file);
    `minimap2 -d $tmp/$ref_file_name.mmi $ref_file `;
}
 

## index reference sequence 
&executeCommand("samtools faidx $ref_file") if (! -e "$ref_file.fai");

if ($file_long)
{
   print "Mapping Long reads\n";
   if ($aligner =~ /bowtie/i){
     `bowtie2 --local $bowtie_options -x $ref_file -fU $file_long -S $outDir/LongReads$$.sam`;
   }
   elsif ($aligner =~ /bwa/i)
   {
     $bwa_options .= ' -x pacbio ' if ($pacbio);
     #&executeCommand("bwa bwasw -M -H $pacbio_bwa_option -t $bwa_threads $ref_file $file_long -f $outDir/LongReads$$.sam");
     &executeCommand("bwa mem $bwa_options $ref_file $file_long | samtools view -@ $samtools_threads -ubS - |  samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/LongReads$$.bam - ");
  #my $mapped_Long_reads=`awk '\$3 !~/*/ && \$1 !~/\@SQ/ {print \$1}' $tmp/LongReads$$.sam | uniq - | wc -l`;
  #`echo -e "Mapped_reads_number:\t$mapped_Long_reads" >>$outDir/LongReads_aln_stats.txt`;
   }elsif ($aligner =~ /minimap2/i){
     `minimap2 -La $minimap2_options  $tmp/$ref_file_name.mmi $file_long > $outDir/LongReads$$.sam`;
   }

   &executeCommand("samtools view -t $ref_file.fai -uhS $outDir/LongReads$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/LongReads$$.bam - ") if ( -s "$outDir/LongReads$$.sam");
}
if ($paired_files){
   print "Mapping paired end reads\n";
   $offset = fastq_utility::checkQualityFormat($file1);
   my $quality_options="";
   if ($aligner =~ /bowtie/i){
     $quality_options= " --phred64 " if ($offset==64);
     &executeCommand("bowtie2 $bowtie_options $quality_options -x $ref_file -1 $file1 -2 $file2 -S $outDir/paired$$.sam");
   }
   elsif($aligner =~ /bwa_short/i)
   {
     $quality_options= " -I " if ($offset==64);
     &executeCommand("bwa aln $bwa_options $quality_options $ref_file $file1 > $tmp/reads_1_$$.sai");
     &executeCommand("bwa aln $bwa_options $quality_options $ref_file $file2 > $tmp/reads_2_$$.sai");
     &executeCommand("bwa sampe $ref_file $tmp/reads_1_$$.sai $tmp/reads_2_$$.sai $file1 $file2 > $outDir/paired$$.sam");
   }
   elsif($aligner =~ /bwa/i)
   {
     &executeCommand("bwa mem $bwa_options $ref_file $file1 $file2 | samtools view -@ $samtools_threads -ubS - | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/paired$$.bam - ");
   }elsif ($aligner =~ /minimap2/i){
     `minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $file1 $file2 > $outDir/paired$$.sam`;
   }
   &executeCommand("samtools view -t $ref_file.fai -uhS $outDir/paired$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/paired$$.bam - ") if (-s "$outDir/paired$$.sam");
}

if ($singleton)  
{
    print "Mapping single end reads\n";
    $offset = fastq_utility::checkQualityFormat($singleton);
    my $quality_options="";

    if ($aligner =~ /bowtie/i){
       $quality_options= " --phred64 " if ($offset==64);
       &executeCommand("bowtie2 $bowtie_options $quality_options -x $ref_file -U $singleton -S $outDir/singleton$$.sam");
    }
    elsif($aligner =~ /bwa_short/i)
    {
      $quality_options= " -I " if ($offset==64);
      &executeCommand("bwa aln $bwa_options $quality_options $ref_file $singleton > $tmp/singleton$$.sai");
      &executeCommand("bwa samse -n 50 $ref_file $tmp/singleton$$.sai $singleton > $outDir/singleton$$.sam");
    }
    elsif($aligner =~ /bwa/i)
    {
      &executeCommand("bwa mem $bwa_options $ref_file $singleton | samtools view -@ $samtools_threads -ubS -| samtools sort -T $tmp -@ $samtools_threads -O BAM -o  $outDir/singleton$$.bam -");
    } elsif ($aligner =~ /minimap2/i){
      `minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $singleton> $outDir/singleton$$.sam`;
    }
    &executeCommand("samtools view -t $ref_file.fai -uhS $outDir/singleton$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o  $outDir/singleton$$.bam - ") if  (-s "$outDir/singleton$$.sam");
} 

if ($file_long and $paired_files and $singleton){
  &executeCommand("samtools merge -f -@ $samtools_threads -h $outDir/paired$$.bam $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam $outDir/LongReads$$.bam");
}
elsif($file_long and $paired_files)
{
  &executeCommand("samtools merge -f -@ $samtools_threads -h $outDir/paired$$.bam $bam_output $outDir/paired$$.bam $outDir/LongReads$$.bam");
}
elsif($paired_files and $singleton)
{
  &executeCommand("samtools merge -f -@ $samtools_threads -h $outDir/paired$$.bam $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam");
}
elsif($singleton and $file_long)
{
  &executeCommand("samtools merge -f -@ $samtools_threads -h $outDir/singleton$$.bam $bam_output $outDir/singleton$$.bam $outDir/LongReads$$.bam");
}
elsif($paired_files)
{
  `mv $outDir/paired$$.bam $bam_output`;
}
elsif($singleton)
{
  `mv $outDir/singleton$$.bam $bam_output`;
}
elsif($file_long)
{
  `mv $outDir/LongReads$$.bam $bam_output`;
}
} # skip_aln
  # get alignment statistical numbers
  &executeCommand("samtools flagstat $bam_output > $stats_output");
  # generate pileup file for coverage calculation
  &executeCommand("samtools mpileup -ABQ0 -d10000000 -f  $ref_file $bam_output >$pileup_output");


  # pull out un-mapped reads list 
  if ($paired_files)  # paired-end reads
  {
    #`samtools view -X $prefix.sort.bam | awk '\$2 ~ /u/ &&  \$2 ~ /2/ {print \$1"/2"} \$2 ~ /u/ && \$2 ~ /1/ {print \$1"/1"}' > $outDir/un_mapped.list`;
  }
  else
  {
    #`samtools view -X $prefix.sort.bam | awk '\$2 ~ /u/ {print \$1}' > $outDir/un_mapped.list`;
  }

  if (system ("ContigCoverageFold_plots_from_samPileup.pl $ref_file $pileup_output $outDir $cov_cut_off 1>>$stats_output"))
  {
    die "$!\n";
  }

  `samtools index $bam_output $bam_index_output`;
## add mapped reads number to coverage table  
   my $mapped_reads_r = &mapped_reads_per_contigs($bam_output);
   open (IN,"$outDir/${prefix}_coverage.table") or die "$!\n";
   open (OUT,">$outDir/update_table$$") or die "$!\n";
   my $header = <IN>;
   chomp $header;
   print OUT $header."\t"."Mapped_reads\n";
   while (<IN>)
   {
       chomp;
       my @array=split /\t/,$_;
       print OUT $_."\t".$mapped_reads_r->{$array[0]}."\n";
   }
   close IN;
   close OUT;
   `mv $outDir/update_table$$ $outDir/${prefix}_coverage.table`;

  # clean up
  `rm -rf $tmp/*$$*`;
  `rm -rf $outDir/*$$*`;
  unlink "$outDir/$prefix.pileup";


sub mapped_reads_per_contigs {
  my $bam_output = shift;
  my %hash;
  open (IN, "samtools idxstats $bam_output |") or die "$!\n";
  while (<IN>)
  {
      chomp;
      my ($id,$len, $mapped,$unmapped)=split /\t/,$_;
      next if ($id eq '*');
      $hash{$id}=$mapped;
  }
  close IN;
  return \%hash;
}

sub executeCommand 
{
    my $command = shift;
    system($command) == 0
         || die "the command $command failed\n";
}




sub splitContigFile
{
    my $file=$_[0];
    open (IN,$file);
    my $seq_number = `grep -c ">" $file`;
    my $mid_point = int($seq_number/2);
    my $n;
    open (OUT,">$tmp/Contig$$.01");
    open (OUT2,">$tmp/Contig$$.02");
    while(<IN>)
    {
       chomp;
       if (/>/)
       {
          $n++;
          if ($n<$mid_point)
          {
              print OUT $_,"\n";
          }
          else
          {
              print OUT2 $_,"\n";
          }
       }
       else
       {
           if ($n<$mid_point)
          {
              print OUT $_,"\n";
          }
          else
          {
              print OUT2 $_,"\n";
          }
       }
    }
    close IN;
    close OUT;
    close OUT2;
}

sub fold {
    # fold and filter reads length by 200 bp.
    my $file=$_[0];
    my $seq;
    my $seq_name;
    my $len_cutoff=199;
    open (IN,$file);
    open (OUT,">$tmp/Contig$$");
    while(<IN>){
      chomp;
      if(/>/)
      {
         if ($seq and length($seq)>$len_cutoff)
         {
           $seq =~ s/ //g;
           $seq =~ s/(.{100})/$1\n/g;
           chomp $seq;
           print OUT $seq_name,"\n",$seq,"\n";
         }
         $seq_name=$_;
         $seq="";
      }
      else
      {
         $seq.=$_;
      }
    }
    if ($seq and length($seq)>$len_cutoff) # last sequence
    {
         $seq =~ s/ //g;
         $seq =~ s/(.{100})/$1\n/g;
         chomp $seq;
         print OUT $seq_name,"\n",$seq,"\n";
    }
    close IN;
    close OUT;
    return ("$tmp/Contig$$");
}
