#!/usr/bin/env perl
######   Required  ########################################################
#          1. Parallel::ForkManager module from CPAN                      #
#          2. String::Approx  module from CPAN                            #
#          3. R for ploting                                               #
#          4. Jellyfish for kmer counting                                 #
#             (http://www.cbcb.umd.edu/software/jellyfish/)               #
###########################################################################
######   Inputs   #########################################################
#          1. Fastq files either paired-end or unpaired reads or both     #
#              Can input multiple library fastq files but only output     #
#              concatenate trimmed fastq files                            #
#          2. Output directory                                            #
#          3. Other options                                               #
###########################################################################
######   Output   #########################################################
#          1.  Two Paired-ends files if input paired-end reads            #
#          2.  One unpaired reads file                                    #
#          3.  trimming statistical text file                             #
#          4.  quality report pdf file                                    #
###########################################################################
######   Functions     ####################################################
#       1. trim bidirection. 2. add -min_L flag                           #
#              3. phredScore= ord(Q)-$ascii for diffent quality           #
#                 encoding and conversion                                 #
#              4. -n "N" base filter   5. Low complexity filter           #
#              6. multi-threads  (required Parallel::ForkManager)         #
#              7. output ascii conversion 8. stats report                 #
#              9. input paired end reads 10. average read quality filter  #
#             11. Trim artifact 12. replace N                             #
# AUTHOR: CHIEN-CHI LO                                                    #                                                
# Copyright (c) 2013 LANS, LLC All rights reserved                        #
# All right reserved. This program is free software; you can redistribute #
# it and/or modify it under the same terms as Perl itself.                #
#                                                                         #
# LAST REVISED: Aug 2014                                                  # 
###########################################################################
use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Parallel::ForkManager;
use String::Approx;

my $version=1.36;
my $debug=0;

$ENV{PATH}="$Bin/../bin:$ENV{PATH}";
sub Usage {
    my $msg=shift;
    my $short_usage = "perl $0 [options] [-u unpaired.fastq] -p reads1.fastq reads2.fastq -d out_directory";
($msg)?
print "\n $msg\n\n $short_usage\n\n Option -h to see full usage \n\n" :
print <<"END";
     Usage: $short_usage
     Version $version
     Input File: (can use more than once)
            -u            <Files> Unpaired reads
            
            -p            <Files> Paired reads in two files and separate by space
     Trim:
            -mode         "HARD" or "BWA" or "BWA_plus" (default BWA_plus)
                          BWA trim is NOT A HARD cutoff! (see bwa's bwa_trim_read() function in bwaseqio.c)

            -q            <INT> Targets # as quality level (default 5) for trimming
    
            -5end         <INT> Cut # bp from 5 end before quality trimming/filtering 
      
            -3end         <INT> Cut # bp from 3 end before quality trimming/filtering 

            -adapter      <bool> Trim reads with illumina adapter/primers (default: no)
                          -rate   <FLOAT> Mismatch ratio of adapters' length (default: 0.2, allow 20% mismatches)
                          -polyA  <bool>  Trim poly A ( > 15 ) 
                          -keepshort  turn on this will keep short portion of reads instead of keep longer portion of reads
                          -keep5end   keep 5' end (conflict with -keepshort and -keep3end)
                          -keep3end   keep 3' end (conflict with -keepshort and -keep5end)
            					
            -artifactFile  <File>    additional artifact (adapters/primers/contaminations) reference file in fasta format 

     Filters:
            -min_L        <INT> Trimmed read should have to be at least this minimum length (default:50)

            -avg_q        <NUM> Average quality cutoff (default:0, no filtering)
            
            -n            <INT> Trimmed read has more than this number of continuous base "N" will be discarded. 
                          (default: 2, "NN") 

            -lc           <FLOAT> Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence  (default: 0.85)

            -phiX         <bool> Filter phiX reads (slow)
            
     Q_Format:
            -ascii        Encoding type: 33 or 64 or autoCheck (default)
                          Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)

            -out_ascii    Output encoding. (default: 33)
     Output:
            -prefix       <TEXT> Output file prefix. (default: QC)

            -stats        <File> Statistical numbers output file (default: prefix.stats.txt)

            -d            <PATH> Output directory.
     Options:
            -t            <INT > # of CPUs to run the script (default:2 )

            -split_size   <INT> Split the input file into several sub files by sequence number (default: 1000000) 

            -qc_only      <bool> no Filters, no Trimming, report numbers.

            -kmer_rarefaction     <bool>   
                          Turn on the kmer calculation. Turn on will slow down ~10 times. (default:Calculation is off.)
                          (meaningless if -subset is too small)
                          -m  <INT>     kmer for rarefaction curve (range:[2,31], default 31)

            -subset       <INT>   Use this nubmer x split_size for qc_only and kmer_rarefaction  
                                  (default: 10,  10x1000000 SE reads, 20x1000000 PE reads)

            -discard      <bool> Output discarded reads to prefix.discard.fastq (default: 0, not output)
 
            -substitute   <bool> Replace "N" in the trimmed reads with random base A,T,C ,or G (default: 0, off)
 
            -trim_only    <bool> No quality report. Output trimmed reads only.

            -replace_to_N_q  <INT>  For NextSeq data, to replace base G to N when below this quality score (default:0, off)

            -debug        <bool> keep intermediate files
END
exit(1);
}

# magic number of quality score
my $highest_illumina_score=41;
my $lowest_illumina_score=0;
# Options Variable initialization
my $thread=2;
my $opt_q=5;
my $opt_min_L=50;
my $opt_avg_cutoff=0;
my $trim_5_end=0;
my $trim_3_end=0;
my $ascii;
my $mode="BWA_plus";
my $N_num_cutoff=2;
my $replace_N;
my $replace_to_N_q=0;
my $is_NextSeq=0;
my $out_offset=33;
my $low_complexity_cutoff_ratio=0.85;
my $subfile_size=1000000;
#my $subfile_size=100000;
my $kmer_rarefaction_on=0; 
my $kmer=31; 
my $prefix="QC";
my $plots_file;
my $stats_output;
my $trimmed_reads1_fastq_file;
my $trimmed_reads2_fastq_file;
my $trimmed_unpaired_fastq_file;
my $trimmed_discard_fastq_file;
my @paired_files;
my @unpaired_files;
my $outDir;
my $output_discard;
my $qc_only=0;
my $trim_only=0;
my $stringent_cutoff=0;
my $filter_adapter=0;
my $keep_short_after_adapter_trim=0;
my $keep_5end_after_adapter_trim=0;
my $keep_3end_after_adapter_trim=0;
my $trim_polyA;
my $filter_phiX=0;
my $filterAdapterMismatchRate=0.2;
my $artifactFile;
my $subsample_num=10;

## not used yet for shihai's algorithm
my $minicut=4;
my $lowperc=0;
my $upperc=0.2;
my $discperc=0.35;
my $cutlimit1=0.3;
my $cutlimit2=0.3;

# Options
GetOptions("q=i"          => \$opt_q,
           "min_L=i"      => \$opt_min_L,
           "avg_q=f"      => \$opt_avg_cutoff,
           "5end=i"       => \$trim_5_end,
           "3end=i"       => \$trim_3_end,
           "mode=s"       => \$mode,
           "p=s{,}"       => \@paired_files,
           "u=s{,}"       => \@unpaired_files,
           "ascii=i"      => \$ascii,
           "n=i"          => \$N_num_cutoff,
           'stringent_q=i'=> \$stringent_cutoff,
           "lc=f"         => \$low_complexity_cutoff_ratio,
           "out_ascii=i"  => \$out_offset,
           "t|threads=i"  => \$thread,
           'kmer_rarefaction' => \$kmer_rarefaction_on,
           'm=i'          => \$kmer,
           'split_size=i' => \$subfile_size,
           'prefix=s'     => \$prefix,
           'd=s'          => \$outDir,
           'stats=s'      => \$stats_output,
           'discard'      => \$output_discard,
           'substitute'   => \$replace_N,
           'qc_only'      => \$qc_only,
           'trim_only'    => \$trim_only,
           'replace_to_N_q=i' => \$replace_to_N_q,
           'subset=i'     => \$subsample_num,
           'debug'        => \$debug,
           'adapter'      => \$filter_adapter,
           'keepshort'    => \$keep_short_after_adapter_trim,
           'keep5end'     => \$keep_5end_after_adapter_trim,
           'keep3end'     => \$keep_3end_after_adapter_trim,
           'polyA'        => \$trim_polyA,
           'phiX'         => \$filter_phiX,
           'rate=f'       => \$filterAdapterMismatchRate,
           'artifactFile=s'  => \$artifactFile,
           'R1=s'         => \$trimmed_reads1_fastq_file,   # for galaxy impelementation
           'R2=s'         => \$trimmed_reads2_fastq_file,   # for galaxy impelementation
           'Ru=s'         => \$trimmed_unpaired_fastq_file, # for galaxy impelementation
           'Rd=s'         => \$trimmed_discard_fastq_file,  # for galaxy impelementation
           'QRpdf=s'      => \$plots_file,                  # for galaxy impelementation
           "version"      => sub{print "Version: $version\n";exit;},
           "help|?"       => sub{Usage()} );


####   Input check  ####
Usage("Missing input files.") unless @unpaired_files or @paired_files;
my @make_paired_paired_files;
my %file;
if (@paired_files)
{
    if (scalar(@paired_files) % 2) { Usage("Please check paired data input are even file numbers\n") ;}
    map { if(is_file_empty($_)){ Usage("Please check paired data input at flag -p.\n $_ doesn't not exist or empty."); $file{basename($_)}=1;} } @paired_files;
    #make pair in a new array 'read1_1 read1_2', 'read2_1 read2_2' ...
    for(my$i=0;$i<=$#paired_files;$i=$i+2)
    {            
        if (&is_paired($paired_files[$i], $paired_files[$i+1]))
        {
            push @make_paired_paired_files, "$paired_files[$i] $paired_files[$i+1]";
        }
        else
        {
            print ("The sequence names of the paired end reads in $paired_files[$i],$paired_files[$i+1] are not matching.\nWill use them as single end reads\n");
            push @unpaired_files, $paired_files[$i],$paired_files[$i+1];
            delete $file{basename($paired_files[$i])};
            delete $file{basename($paired_files[$i+1])};
        }
    }
}

if (@unpaired_files)
{
    map { if(is_file_empty($_))
          { 
              Usage("Please check unpaired data input at flag -u.\n $_ doesn't not exist or empty.");
          } 
          if ($file{basename($_)}) 
          {
              Usage("The single end file, $_,has been used in the paired end data.")
          }
        } @unpaired_files;
}
Usage("Missing output directory at flag -d  ") unless $outDir;

#########################


if ($mode =~ /hard/i)
{
    print "Hard trimming is used. \n" if (!$qc_only);
    $mode="hard";
}
elsif ($mode =~ /BWA_plus/i)
{
    print "Bwa extension trimming algorithm is used. \n" if (!$qc_only);
    $mode="BWA_plus";
}
elsif ($mode =~ /BWA/i)
{
    print "Bwa trimming algorithm is used. \n" if (!$qc_only);
    $mode="BWA";
}
else   # default
{
    print "Not recognized mode $mode. Bwa extension trimming algorithm is used. \n" if (!$qc_only);
    $mode="BWA_plus";
}

###### Output file initialization #####
# temp files for plotting
my $quality_matrix="$outDir/$prefix.quality.matrix";
my $avg_quality_histogram="$outDir/$prefix.for_qual_histogram.txt";
my $base_matrix="$outDir/$prefix.base.matrix";
my $nuc_composition_file="$outDir/$prefix.base_content.txt";
my $length_histogram="$outDir/$prefix.length_count.txt";
my $qa_quality_matrix="$outDir/qa.$prefix.quality.matrix";
my $qa_avg_quality_histogram="$outDir/qa.$prefix.for_qual_histogram.txt";
my $qa_base_matrix="$outDir/qa.$prefix.base.matrix";
my $qa_nuc_composition_file="$outDir/qa.$prefix.base_content.txt";
my $qa_length_histogram="$outDir/qa.$prefix.length_count.txt";
my $fastq_count="$outDir/fastqCount.txt";

# output files
$trimmed_discard_fastq_file="$outDir/$prefix.discard.fastq" if (!$trimmed_discard_fastq_file);
$plots_file="$outDir/${prefix}_qc_report.pdf" if (!$plots_file);
$trimmed_unpaired_fastq_file="$outDir/$prefix.unpaired.trimmed.fastq" if (!$trimmed_unpaired_fastq_file);
$trimmed_reads1_fastq_file="$outDir/$prefix.1.trimmed.fastq" if (!$trimmed_reads1_fastq_file);
$trimmed_reads2_fastq_file="$outDir/$prefix.2.trimmed.fastq" if (!$trimmed_reads2_fastq_file);
$stats_output="$outDir/$prefix.stats.txt" if (!$stats_output);
   
# temp files for kmer counting plot   
my $kmer_rarefaction_file="$outDir/$prefix.Kmercount.txt"; #
my $kmer_files="$outDir/$prefix.KmerFiles.txt"; #
my $kmer_histogram_file="$outDir/$prefix.kmerH.txt"; #
#######################################

######   Output check  ################
if (! -e $outDir)
{
  mkdir $outDir;
}

if (-e $trimmed_reads1_fastq_file)
{
   print "The output $trimmed_reads1_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads1_fastq_file");
}
if (-e $trimmed_reads2_fastq_file)
{
   print "The output $trimmed_reads2_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_reads2_fastq_file");
}
if (-e $trimmed_unpaired_fastq_file)
{
   print "The output $trimmed_unpaired_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_unpaired_fastq_file");
}
if (-e $plots_file)
{
   print "The output $plots_file file exists and will be overwritten.\n";
   system ("rm $plots_file");
}
if (-e $trimmed_discard_fastq_file)
{
   print "The output $trimmed_discard_fastq_file file exists and will be overwritten.\n";
   system ("rm $trimmed_discard_fastq_file");
}
system ("rm $kmer_files") if (-e $kmer_files);
#######################################

if ($qc_only)
{
    # set all Filters to pass
   # $opt_q=0;
   # $opt_min_L=0;
   # $opt_avg_cutoff=0;
    #$N_num_cutoff=1000;
    #$low_complexity_cutoff_ratio=1;
}

my $orig_opt_q = $opt_q;
my ( $total_count,$total_count_1, $total_count_2, $total_len_1,$total_len_2, $total_len);
my ( $total_num, $trimmed_num,$total_raw_seq_len, $total_trimmed_seq_len);
my ( $trim_seq_len_std, $trim_seq_len_avg, $max, $min, $median);
my ( $paired_seq_num, $total_paired_bases )=(0,0);
my ( @split_files, @split_files_2);
my ( $readsFilterByLen, $basesFilterByLen ) = (0,0);
my ( $readsTrimByQual, $basesTrimByQual ) = (0,0);
my ( $readsFilterByNN, $basesFilterByNN ) = (0,0);
my ( $readsFilterByPhiX, $basesFilterByPhiX ) = (0,0);
my ( $readsTrimByAdapter, $basesTrimByAdapter ) = (0,0);
my ( $readsFilterByAvgQ, $basesFilterByAvgQ ) = (0,0);
my ( $readsFilterByLowComplexity, $basesFilterByLowComplexity ) = (0,0);
$N_num_cutoff=1 if (!$N_num_cutoff);
my (%position, %qa_position);
my (%AverageQ, %qa_AverageQ);
my (%base_position,%qa_base_position);
my (%base_content, %qa_base_content);
my (%len_hash, %qa_len_hash);
my (%filter_stats, %qa_filter_stats);
my %EachAdapter;
my %EachReplaceN;
my ( $i_file_name, $i_path, $i_suffix );
  
my ($phiX_id,$phiX_seq) = &read_phiX174 if ($filter_phiX);
$filter_adapter=1 if ($artifactFile);
$filter_adapter=1 if ($trim_polyA);

open(my $fastqCount_fh, ">$fastq_count") or die "Cannot write $fastq_count\n";
  foreach my $input (@unpaired_files,@make_paired_paired_files){
     print "Processing $input file\n";
     #print $STATS_fh "Processing $input file\n";
     my ($reads1_file,$reads2_file) = split /\s+/,$input;
  
     # check file 
     if(&is_file_empty($reads1_file)<0) { die "The file $reads1_file doesn't exist or empty.\n";}
     if(&is_file_empty($reads2_file)<0 and $reads2_file) { die "The file $reads2_file doesn't exist or empty.\n";}

     # check quality offset
     if (! $ascii){$ascii = &checkQualityFormat($reads1_file)}

     # check NextSeq platform
     $is_NextSeq = ( &is_NextSeq($reads1_file) )?1:0;
     if( $is_NextSeq ){
        if ($orig_opt_q < 20){
	  $opt_q = 20;
	  warn "The input looks like NextSeq data and the quality level (-q) is adjusted to $opt_q for trimming.\n";
        }
     }else{ $opt_q = $orig_opt_q;}

    #split
    ($total_count_1,$total_len_1,@split_files) = &split_fastq($reads1_file,$outDir,$subfile_size);
    printf $fastqCount_fh ("%s\t%d\t%d\t%.2f\n",basename($reads1_file),$total_count_1,$total_len_1,$total_len_1/$total_count_1);
    if ($reads2_file)
    {
        ($total_count_2,$total_len_2,@split_files_2) = &split_fastq($reads2_file,$outDir,$subfile_size);
        printf $fastqCount_fh ("%s\t%d\t%d\t%.2f\n",basename($reads2_file),$total_count_2,$total_len_2,$total_len_2/$total_count_2);
    }
     $total_count += $total_count_1 + $total_count_2;
     $total_len += $total_len_1 + $total_len_2;
     my $random_num_ref = &random_subsample(scalar(@split_files),$subsample_num);
     $subsample_num -= scalar(@split_files);
     $subsample_num -= scalar(@split_files) if ($reads2_file);

    my $pm = new Parallel::ForkManager($thread);
 # data structure retrieval and handling from multi-threading 
   # $pm->run_on_start(
   #   sub { my ($pid,$ident)=@_;
   #    #print "** $ident started, pid: $pid\n";
   #   }
   # );

    $pm -> run_on_finish ( # called BEFORE the first call to start()
      sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $nums_ref) = @_;
        
        print qq|child process $ident, exit on singal $exit_signal\n| if ($exit_signal);
        print qq|core dump from child process $pid! on $ident\n| if ($core_dump);
        
        # retrieve data structure from child
        if (defined($nums_ref)) {  # children are not forced to send anything
          my $random_sub_raw_ref;
          $random_sub_raw_ref=$nums_ref->{random_sub_raw} if ( $random_num_ref->{$ident} and !$qc_only);
          my ($total_score, $avg_score);
          $total_num += $nums_ref->{raw_seq_num};
          $trimmed_num += $nums_ref->{trim_seq_num};
          $total_raw_seq_len += $nums_ref->{total_raw_seq_len};
          $total_trimmed_seq_len += $nums_ref->{total_trim_seq_len};
          my $processed_num= $nums_ref->{raw_seq_num};
          $trim_seq_len_std = $nums_ref->{trim_seq_len_std};
          $trim_seq_len_avg = $nums_ref->{trim_seq_len_avg};
          $max = $nums_ref->{max};
          $min = $nums_ref->{min};
          $median = $nums_ref->{median};
          $paired_seq_num += $nums_ref->{paired_seq_num};
          $total_paired_bases +=  $nums_ref->{total_paired_bases};
          my %filter;
          if ($nums_ref->{filter})
          {
            %filter=%{$nums_ref->{filter}};
            map {$filter_stats{$_}->{basesNum} += $filter{$_}->{basesNum};
                 $filter_stats{$_}->{readsNum} += $filter{$_}->{readsNum}; 
              } keys %filter; 
          }
          if ($replace_N)
          {
              my %tmp_ReplaceN;
              %tmp_ReplaceN  = %{$nums_ref->{replaceN}} if ($nums_ref->{replaceN});
              map {
                    $EachReplaceN{$_} += $tmp_ReplaceN{$_};
                  } keys %tmp_ReplaceN;
          }
          if ($filter_adapter)
          {
              my %tmp_Adapters;
              %tmp_Adapters= %{$nums_ref->{adapter}} if ($nums_ref->{adapter});
              map {
                    $EachAdapter{$_}->{readsNum} += $tmp_Adapters{$_}->{readsNum};
                    $EachAdapter{$_}->{basesNum} += $tmp_Adapters{$_}->{basesNum};
                  } keys %tmp_Adapters;
          }
          if ($nums_ref->{ReadAvgQ})
          {
              my %temp_avgQ = %{$nums_ref->{ReadAvgQ}};
              map {$AverageQ{$_}->{bases} += $temp_avgQ{$_}->{basesNum};
                   $AverageQ{$_}->{reads} += $temp_avgQ{$_}->{readsNum}; 
              } keys %temp_avgQ; 
          }
          my %temp_position= %{$nums_ref->{qual}} if ($nums_ref->{qual});
          my %temp_base_position= %{$nums_ref->{base}} if ($nums_ref->{base});
          my %qa_temp_position= %{$random_sub_raw_ref->{qual}} if ($random_sub_raw_ref);
          my %qa_temp_base_position= %{$random_sub_raw_ref->{base}} if ($random_sub_raw_ref);

          for my $score ($lowest_illumina_score..$highest_illumina_score)
          {   
              if ($nums_ref->{qual})
              {
                  foreach my $pos (keys (%temp_position))
                  {
                      $position{$pos}->{$score} += $temp_position{$pos}->{$score};
                      $total_score +=  $score * $temp_position{$pos}->{$score};
                  }
              }
              if ($random_sub_raw_ref)
              {
                  foreach my $pos (keys (%qa_temp_position))
                  {
                      $qa_position{$pos}->{$score} += $qa_temp_position{$pos}->{$score};
                  }
              }
          }
          for my $nuc ("A","T","C","G","N")
          {
              if ($nums_ref->{base})
              {
                  foreach my $pos (keys (%temp_base_position))
                  {
                      $base_position{$pos}->{$nuc} += $temp_base_position{$pos}->{$nuc};
                  }
              }
              if ($random_sub_raw_ref)
              {
                  foreach my $pos (keys (%qa_temp_base_position))
                  {
                      $qa_base_position{$pos}->{$nuc} += $qa_temp_base_position{$pos}->{$nuc};
                  }
              }
          }
          
          my %tmp_base_content = %{$nums_ref->{Base_content}} if ($nums_ref->{Base_content});
          my %qa_tmp_base_content = %{$random_sub_raw_ref->{Base_content}} if ($random_sub_raw_ref);
          for my $nuc ("A","T","C","G","N","GC")
          {
              while (my ($key, $value)= each %{$tmp_base_content{$nuc}})
              {
                 $base_content{$nuc}->{$key} += $value if ($nums_ref->{Base_content});
              }
              if ($random_sub_raw_ref)
              {
                 while (my ($key, $value)= each %{$qa_tmp_base_content{$nuc}})
                 {
                    $qa_base_content{$nuc}->{$key} += $value;
                 }
              }
          } 
          if ($nums_ref->{ReadLen})
          {
              while (my ($key, $value)= each %{$nums_ref->{ReadLen}} )
              {
		         $len_hash{$key} += $value;   
              }
          }
          if ( $random_num_ref->{$ident}){
              open (KMEROUT, ">>$kmer_files") or die "$!\n";
              if ($qc_only)
              {
                  print KMEROUT $nums_ref->{file_name_1},"\t",$nums_ref->{raw_seq_num_1},"\n";
                  print KMEROUT $nums_ref->{file_name_2},"\t" if ($nums_ref->{file_name_2});
                  print KMEROUT $nums_ref->{raw_seq_num_2},"\n" if ($nums_ref->{file_name_2});
              }
              else
              {
                  print KMEROUT $nums_ref->{trim_file_name_1},"\t",$nums_ref->{trim_seq_num_1},"\n";
                  print KMEROUT $nums_ref->{trim_file_name_2},"\t" if ($nums_ref->{trim_file_name_2});
                  print KMEROUT $nums_ref->{trim_seq_num_2},"\n" if ($nums_ref->{trim_file_name_2});
              }
              close KMEROUT; 
              if ($random_sub_raw_ref)
              {
                  while (my ($key, $value)= each %{$random_sub_raw_ref->{ReadLen}} )
                  {
		         $qa_len_hash{$key} += $value;   
                  }
                  my %qa_temp_avgQ = %{$random_sub_raw_ref->{ReadAvgQ}};
                  map {$qa_AverageQ{$_}->{bases} += $qa_temp_avgQ{$_}->{basesNum};
                       $qa_AverageQ{$_}->{reads} += $qa_temp_avgQ{$_}->{readsNum}; 
                  } keys %qa_temp_avgQ; 
 
                  if ($random_sub_raw_ref->{filter})
                  {
                      %filter=%{$random_sub_raw_ref->{filter}};
                      map {$qa_filter_stats{$_}->{basesNum} += $filter{$_}->{basesNum};
                           $qa_filter_stats{$_}->{readsNum} += $filter{$_}->{readsNum}; 
                          } keys %filter; 
                  }
              }
              # print Dumper($random_sub_raw_ref),"\n";
          }

          #print $STATS_fh " Processed $total_num/$total_count\n";
          print "Processed $total_num/$total_count\n";
          if ( $nums_ref->{total_trim_seq_len} )
          {
              if (!$trim_only)
              {
                  printf (" Post Trimming Length(Mean, Std, Median, Max, Min) of %d reads with Overall quality %.2f\n",$nums_ref->{trim_seq_num}, $total_score/$nums_ref->{total_trim_seq_len});
              }
              else
              {
                  printf (" Post Trimming Length(Mean, Std, Median, Max, Min) of %d reads\n",$nums_ref->{trim_seq_num});
              }
              printf (" (%.2f, %.2f, %.1f, %d, %d)\n",$trim_seq_len_avg,$trim_seq_len_std,$median,$max,$min);
          #unlink $split_files[$ident];
          #unlink $split_files_2[$ident] if ($split_files_2[$ident]);
          }
          else
          {
              print "All reads are trimmed/filtered\n";
          } 
        } else {  # problems occuring during storage or retrieval will throw a warning
          print qq|No message received from child process $pid! on $ident\n|;
        }
      }
    );

  foreach my $i(0..$#split_files)
  {
    next if ($qc_only and ! $random_num_ref->{$i});
    $pm->start($i) and next;
    my $hash_ref = &qc_process($split_files[$i],$split_files_2[$i],$random_num_ref->{$i});
    &run_kmercount($hash_ref,$kmer,1) if ( $random_num_ref->{$i} and $kmer_rarefaction_on );
    &run_kmercount($hash_ref,$kmer,2) if ( $random_num_ref->{$i} and $kmer_rarefaction_on and $split_files_2[$i]);
    $pm->finish(0, $hash_ref);
   
  }
  $pm->wait_all_children;
  # clean up
  foreach my $i(0..$#split_files)
  {
    unlink $split_files[$i];
    unlink $split_files_2[$i] if ($split_files_2[$i]);
  }
} #end foreach $input

close $fastqCount_fh;

# concatenate each thread's trimmed reads and files clean up.
if (! $qc_only)
{
    if (@unpaired_files) 
    {
      foreach my $input (@unpaired_files){
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
    
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if ($output_discard){
            if (system("cat $outDir/${i_file_name}_?????_trim_discard.fastq >> $trimmed_discard_fastq_file")) { die "cat failed: $!" }
            `rm $outDir/${i_file_name}_?????_trim_discard.fastq`;
        }
      }
    }
    
    if (@make_paired_paired_files)
    {
      foreach my $input (@make_paired_paired_files){
        my ($reads1_file,$reads2_file) = split /\s+/,$input;
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads1_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads1_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if ($output_discard){ 
           if (system("cat $outDir/${i_file_name}_?????_trim_discard.fastq >> $trimmed_discard_fastq_file")) { die "cat failed: $!" }
           `rm $outDir/${i_file_name}_?????_trim_discard.fastq`;
        }
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$reads2_file", qr/\.[^.]*/ );
        ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$i_file_name", qr/\.[^.]*/ ) if ($i_suffix =~ /gz$/);
        if (system("cat $outDir/${i_file_name}_?????_trim.fastq >> $trimmed_reads2_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim.fastq`;
        if (system("cat $outDir/${i_file_name}_?????_trim_unpaired.fastq >> $trimmed_unpaired_fastq_file")) { die "cat failed: $!" }
        `rm $outDir/${i_file_name}_?????_trim_unpaired.fastq`;
      }
    }  

&print_quality_report_files($qa_quality_matrix,$qa_avg_quality_histogram,$qa_base_matrix,$qa_nuc_composition_file,$qa_length_histogram,\%qa_position,\%qa_AverageQ,\%qa_base_position,\%qa_base_content,\%qa_len_hash) if (!$trim_only);
}



&Kmer_rarefaction() if ($kmer_rarefaction_on); 
&print_quality_report_files($quality_matrix,$avg_quality_histogram,$base_matrix,$nuc_composition_file,$length_histogram,\%position,\%AverageQ,\%base_position,\%base_content,\%len_hash) if (!$trim_only);
&print_final_stats();
&plot_by_R() if (!$trim_only);

unless ($debug){
unlink $nuc_composition_file;
unlink $quality_matrix;
unlink $base_matrix;
unlink $avg_quality_histogram;
unlink $length_histogram;
unlink $qa_nuc_composition_file;
unlink $qa_quality_matrix;
unlink $qa_base_matrix;
unlink $qa_avg_quality_histogram;
unlink $qa_length_histogram;
unlink $kmer_files;
    if ($kmer_rarefaction_on)
    {
#        unlink $kmer_rarefaction_file;
#        unlink $kmer_histogram_file;
    }
}

# END MAIN
exit(0);

sub build_initial_quality_matrix  # not used yet
{
    my $fastqFile=shift;
    open (my $fh, "$fastqFile") or die "$fastqFile $!";
    my %basequal;
    my $total_reads=0;
    while(<$fh>)
    {
        last if ($total_reads>2000000);
        my $name=$_;
        my $seq=<$fh>;
        $seq =~ s/\n//g;
        while ($seq !~ /\+/)
        {
           $seq .= <$fh>;
           $seq =~ s/\n//g;
        }
        my $q_id_pos=index($seq,"+");
        $seq = substr($seq, 0, $q_id_pos);
        my $seq_len = length $seq;
        my $qual_seq=<$fh>;
        $qual_seq =~ s/\n//g;
        my $qual_seq_len = length $qual_seq;
        while ( $qual_seq_len < $seq_len )
        {
           last if ( $qual_seq_len == $seq_len);
           $qual_seq .= <$fh>;
           $qual_seq =~ s/\n//g;
           $qual_seq_len = length $qual_seq;
        }
        my @qual_seq=split //, $qual_seq;

        if (rand() <=0.2) {
           $total_reads++;
           for my $pos(0..$#qual_seq)
           {
              push @{$basequal{$pos}}, ord($qual_seq[$pos])-$ascii;
           }
        }
    }
    close $fh;
    @{$basequal{$_}}= sort {$a <=> $b} @{$basequal{$_}} foreach (keys %basequal);

    return \%basequal;
}

sub run_kmercount 
{
    my $hash_ref =shift;
    my $kmer =shift;
    my $mate = shift;
    my $file= ($qc_only)?$hash_ref->{"file_name_$mate"}:$hash_ref->{"trim_file_name_$mate"};
    my $cmd = "jellyfish count -C -s 512M -c 3 -m $kmer -o $file$prefix $file";
    if (system($cmd)){
      die "$!\n";
    }
    
}  


sub Kmer_rarefaction {                                                               
    my $old_merge;                                                                     
    my $count=0;                                                                       
    my $new_merge="$outDir/${prefix}tmpKmerMerge${$}$count";                           
    my (@kmercount,$distinct_kmer,$total_kmer);                          
    #opendir (DIR, $outDir);                                                            
    #my @files=grep { /${prefix}_0/ } readdir(DIR);                                     
    #closedir DIR;  
    open (KMER, ">$kmer_rarefaction_file") or die "$!\n";                              
    open (FILES, "$kmer_files") or die "$!\n";
    while (<FILES>)
    #for my $i(0..$#files)                                                              
    {                                          
      chomp;                                        
      #my $subset_kmer_file="$outDir/$files[$i]";                         
      my ($subset_kmer_file, $seq_num) = split /\t/, $_; 
      $subset_kmer_file = $subset_kmer_file.$prefix."_0";
      #my ($subset_file_name) = $subset_kmer_file =~ /(\S+)${prefix}_0/; 
      #$sequence_num += $subfile_size_hash{$subset_file_name};                          
      if ($count==0)                                                                   
      {                                                                                
         @kmercount=`jellyfish stats $subset_kmer_file | grep -A 1 'Distinct' | awk '{print \$2}'`;                                                                         
         ($distinct_kmer) = $kmercount[0] =~ /(\d+)/;                                  
         ($total_kmer) = $kmercount[1] =~ /(\d+)/;                                     
         print KMER $seq_num . "\t". $distinct_kmer. "\t". $total_kmer."\n";      
         $old_merge=$subset_kmer_file;                                                 
         $count++;                                                                     
      }                                                                                
      else                                                                             
      {                                                                                
         `jellyfish merge -o $new_merge $old_merge $subset_kmer_file `;                
         @kmercount=`jellyfish stats $new_merge | grep -A 1 'Distinct' | awk '{print \$2}'`;                                                                                
         ($distinct_kmer) = $kmercount[0] =~ /(\d+)/;                                  
         ($total_kmer) = $kmercount[1] =~ /(\d+)/;                                     
         unlink $old_merge;                                                            
         unlink $subset_kmer_file;                                                     
         $old_merge=$new_merge;                                                        
         $new_merge="$outDir/${prefix}tmpKmerMerge${$}$count";                         
         print KMER $seq_num . "\t". $distinct_kmer. "\t". $total_kmer."\n";      
         $count++;                                                                     
      }                                                                                
    }                                     
    close FILES;                                             
    close KMER;                                                                        
    `jellyfish histo -o $kmer_histogram_file $old_merge`;                              
    unlink $old_merge;                                                                 
} 

sub print_final_stats{
    open (my $fh, ">$stats_output") or die "$!\t$stats_output\n";
    $readsFilterByLen = $filter_stats{len}->{readsNum};
    $basesFilterByLen = $filter_stats{len}->{basesNum};
    $readsTrimByQual = $filter_stats{qualTrim}->{readsNum};
    $basesTrimByQual = $filter_stats{qualTrim}->{basesNum};
    $readsFilterByNN = $filter_stats{NN}->{readsNum};
    $basesFilterByNN = $filter_stats{NN}->{basesNum};
    $readsFilterByPhiX = $filter_stats{phiX}->{readsNum};
    $basesFilterByPhiX = $filter_stats{phiX}->{basesNum};
    $readsTrimByAdapter = $filter_stats{adapter}->{readsNum};
    $basesTrimByAdapter = $filter_stats{adapter}->{basesNum};
    $readsFilterByAvgQ = $filter_stats{AvgQ}->{readsNum};
    $basesFilterByAvgQ = $filter_stats{AvgQ}->{basesNum};
    $readsFilterByLowComplexity = $filter_stats{lowComplexity}->{readsNum};
    $basesFilterByLowComplexity = $filter_stats{lowComplexity}->{basesNum};
    # stats
    if ($qc_only)
    {
      print $fh "\n";
      print $fh "Reads #: $total_count\n";
      print $fh "Total bases: $total_len\n";
      printf $fh ("Reads Length: %.2f\n",$total_len/$total_count);
      print $fh "Processed $total_num reads for quality check only\n";
      
      printf $fh ("  Reads length < %d bp: %d (%.2f %%)\n", $opt_min_L, $readsFilterByLen , ($readsFilterByLen)/$total_num*100);
      printf $fh ("  Reads have %d continuous base \"N\": %d (%.2f %%)\n", $N_num_cutoff, $readsFilterByNN , ($readsFilterByNN)/$total_num*100);
      printf $fh ("  Low complexity Reads  (>%.2f%% mono/di-nucleotides): %d (%.2f %%)\n", $low_complexity_cutoff_ratio*100, $readsFilterByLowComplexity , ($readsFilterByLowComplexity)/$total_num*100);
      printf $fh ("  Reads < average quality %.1f: %d (%.2f %%)\n", $opt_avg_cutoff, $readsFilterByAvgQ , ($readsFilterByAvgQ)/$total_num*100);
      printf $fh ("  Reads hits to phiX sequence: %d (%.2f %%)\n", $readsFilterByPhiX , ($readsFilterByPhiX)/$total_num*100) if ($filter_phiX);
      if ($filter_adapter)
      {
        printf $fh ("  Reads with Adapters/Primers: %d (%.2f %%)\n", $readsTrimByAdapter , ($readsTrimByAdapter)/$total_num*100);
        foreach my $adapter_id (keys %EachAdapter)
        {
            my $affect_reads = $EachAdapter{$adapter_id}->{readsNum};
            my $affect_bases = $EachAdapter{$adapter_id}->{basesNum};
            printf $fh ("    %s %d reads (%.2f %%) %d bases (%.2f %%)\n", $adapter_id, $affect_reads, $affect_reads/$total_num*100,$affect_bases, $affect_bases/$total_raw_seq_len*100);
        }
      }
    } 
    else
    {
 #     print $fh "\nQC process\n";
      print $fh "Before Trimming\n";
      print $fh "Reads #: $total_num\n";
      print $fh "Total bases: $total_raw_seq_len\n";
      printf $fh ("Reads Length: %.2f\n",$total_raw_seq_len/$total_num);
    
      print $fh "\nAfter Trimming\n";
      printf $fh ("Reads #: %d (%.2f %%)\n",$trimmed_num, $trimmed_num/$total_num*100);
      printf $fh ("Total bases: %d (%.2f %%)\n",$total_trimmed_seq_len,$total_trimmed_seq_len/$total_raw_seq_len*100);
      if ($trimmed_num)
      {
        printf $fh ("Mean Reads Length: %.2f\n",$total_trimmed_seq_len/$trimmed_num); 
      }
      else
      {
        printf $fh "Mean Reads Length: 0\n";
      }
    
      if (@make_paired_paired_files){
        printf $fh ("  Paired Reads #: %d (%.2f %%)\n",$paired_seq_num, $paired_seq_num/$trimmed_num*100);
        printf $fh ("  Paired total bases: %d (%.2f %%)\n",$total_paired_bases,$total_paired_bases/$total_trimmed_seq_len*100);
        printf $fh ("  Unpaired Reads #: %d (%.2f %%)\n", $trimmed_num - $paired_seq_num, ($trimmed_num - $paired_seq_num)/$trimmed_num*100);
        printf $fh ("  Unpaired total bases: %d (%.2f %%)\n", $total_trimmed_seq_len - $total_paired_bases , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100);
      }
    
      printf $fh ("\nDiscarded reads #: %d (%.2f %%)\n", $total_num - $trimmed_num , ($total_num - $trimmed_num)/$total_num*100);
      printf $fh ("Trimmed bases: %d (%.2f %%)\n", $total_raw_seq_len - $total_trimmed_seq_len, ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100);
      printf $fh ("  Reads Filtered by length cutoff (%d bp): %d (%.2f %%)\n", $opt_min_L, $readsFilterByLen , ($readsFilterByLen)/$total_num*100);
      printf $fh ("  Bases Filtered by length cutoff: %d (%.2f %%)\n", $basesFilterByLen , ($basesFilterByLen)/$total_raw_seq_len*100);
      printf $fh ("  Reads Filtered by continuous base \"N\" (%d): %d (%.2f %%)\n", $N_num_cutoff, $readsFilterByNN , ($readsFilterByNN)/$total_num*100);
      printf $fh ("  Bases Filtered by continuous base \"N\": %d (%.2f %%)\n", $basesFilterByNN , ($basesFilterByNN)/$total_raw_seq_len*100);
      printf $fh ("  Reads Filtered by low complexity ratio (%.1f): %d (%.2f %%)\n", $low_complexity_cutoff_ratio, $readsFilterByLowComplexity , ($readsFilterByLowComplexity)/$total_num*100);
      printf $fh ("  Bases Filtered by low complexity ratio: %d (%.2f %%)\n", $basesFilterByLowComplexity , ($basesFilterByLowComplexity)/$total_raw_seq_len*100);
      if ($opt_avg_cutoff>0)
      {
        printf $fh ("  Reads Filtered by avg quality (%.1f): %d (%.2f %%)\n", $opt_avg_cutoff, $readsFilterByAvgQ , ($readsFilterByAvgQ)/$total_num*100);
        printf $fh ("  Bases Filtered by avg quality: %d (%.2f %%)\n", $basesFilterByAvgQ , ($basesFilterByAvgQ)/$total_raw_seq_len*100);
      }
      if ($filter_phiX)
      {
        printf $fh ("  Reads Filtered by phiX sequence: %d (%.2f %%)\n", $readsFilterByPhiX , ($readsFilterByPhiX)/$total_num*100);
        printf $fh ("  Bases Filtered by phiX sequence: %d (%.2f %%)\n", $basesFilterByPhiX , ($basesFilterByPhiX)/$total_raw_seq_len*100);
      }
      printf $fh ("  Reads Trimmed by quality (%.1f): %d (%.2f %%)\n", $opt_q, $readsTrimByQual , ($readsTrimByQual)/$total_num*100);
      printf $fh ("  Bases Trimmed by quality: %d (%.2f %%)\n", $basesTrimByQual , ($basesTrimByQual)/$total_raw_seq_len*100);
      if ($trim_5_end)
      {
        printf $fh ("  Reads Trimmed with %d bp from 5' end\n", $trim_5_end);
      }
      if ($trim_3_end)
      {
        printf $fh ("  Reads Trimmed with %d bp from 3' end\n", $trim_3_end);
      }
      if ($filter_adapter){
        printf $fh ("  Reads Trimmed with Adapters/Primers: %d (%.2f %%)\n", $readsTrimByAdapter , ($readsTrimByAdapter)/$total_num*100);
        printf $fh ("  Bases Trimmed with Adapters/Primers: %d (%.2f %%)\n", $basesTrimByAdapter , ($basesTrimByAdapter)/$total_raw_seq_len*100);
        foreach my $adapter_id (keys %EachAdapter)
        {
            my $affect_reads = $EachAdapter{$adapter_id}->{readsNum};
            my $affect_bases = $EachAdapter{$adapter_id}->{basesNum};
            printf $fh ("    %s %d reads (%.2f %%) %d bases (%.2f %%)\n", $adapter_id, $affect_reads, $affect_reads/$total_num*100,$affect_bases, $affect_bases/$total_raw_seq_len*100);
        }
      }
      if ($replace_N)
      {
          printf $fh ("\nN base random substitution: A %d, T %d, C %d, G %d\n",$EachReplaceN{A},$EachReplaceN{T},$EachReplaceN{C},$EachReplaceN{G});
      }
    } # end qc only
    

    return (0);
}


sub plot_by_R
{

    open (R,">$outDir/tmp$$.R");
    print R <<RSCRIPT; 
if(file.exists(\"$qa_length_histogram\")){
  pdf(file = \"$plots_file\",width=15,height=7)
}else{
  pdf(file = \"$plots_file\",width=10,height=8)
}
def.par <- par(no.readonly = TRUE) # get default parameters

#Summary
par(family="mono")
SummaryStats<-readLines("$stats_output")
plot(0:1,0:1,type=\'n\',xlab=\"\",ylab=\"\",xaxt=\'n\',yaxt=\'n\',bty=\'n\')
if ($qc_only){
  for (i in 1:length(SummaryStats)){
     text(0.05,1-0.04*(i-1),SummaryStats[i],adj=0,font=2,cex=1)
  }
}else{
  if ($paired_seq_num > 0) {
      adjust <-14
      abline(h=0.73,lty=2)
  }else{
      adjust<-11
      abline(h=0.85,lty=2)
  }
  for (i in 1:length(SummaryStats)){
     if (i>5 && i<adjust){
       text(0.45,1-0.035*(i-6),SummaryStats[i],adj=0,font=2,cex=0.9)
     }else if(i >=adjust){
       text(0.05,1-0.035*(i-6),SummaryStats[i],adj=0,font=2,cex=0.9)
     }else{
       text(0.05,1-0.035*(i-1),SummaryStats[i],adj=0,font=2,cex=0.9)
     }
  }
}
#title(paste(\"$prefix\",\"QC report\"),sub = 'DOE Joint Genome Institute/Los Alamos National Laboratory', adj = 0.5, col.sub='darkblue',font.sub=2,cex.sub=0.8)
title("QC stats")
par(def.par)#- reset to default

#lenght histogram
length_histogram <- function(length_count_file, xlab,ylab){
  lengthfile<-read.table(file=length_count_file)
  lengthList<-as.numeric(lengthfile\$V1)
  lengthCount<-as.numeric(lengthfile\$V2)
  lenAvg<-sum(lengthList * lengthCount)/sum(lengthCount)
  lenStd<-sqrt(sum(((lengthList - lenAvg)**2)*lengthCount)/sum(lengthCount))
  lenMax<-max(lengthList[lengthCount>0])
  lenMin<-min(lengthList[lengthCount>0])
  totalReads<-sum(lengthCount)
  barplot(lengthCount/1000000,names.arg=lengthList,xlab=xlab,ylab=ylab,cex.names=0.8)
  legend.txt<-c(paste("Mean",sprintf ("%.2f",lenAvg),"±",sprintf ("%.2f",lenStd)),paste("Max",lenMax),paste("Min",lenMin))
  legend('topleft',legend.txt,bty='n')
  if (totalReads< $trimmed_num){
      mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
  }
  return(totalReads)
}
if(file.exists(\"$qa_length_histogram\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    qa.readsCount<-length_histogram(\"$qa_length_histogram\","Input Length","Count (millions)")
    readsCount<-length_histogram(\"$length_histogram\","Trimmed Length","")
}else{
    readsCount<-length_histogram(\"$length_histogram\","Length","Count (millions)")
}
par(def.par)#- reset to default
title("Reads Length Histogram")

#readGC plot

readGC_plot <- function(base_content_file, totalReads, fig_x_start,fig_x_end,xlab,ylab,new){
	baseP<-read.table(file=base_content_file)
	Apercent<-baseP\$V2[which(baseP\$V1=="A")]
	ApercentCount<-baseP\$V3[which(baseP\$V1=="A")]
	Tpercent<-baseP\$V2[which(baseP\$V1=="T")]
	TpercentCount<-baseP\$V3[which(baseP\$V1=="T")]
	Cpercent<-baseP\$V2[which(baseP\$V1=="C")]
	CpercentCount<-baseP\$V3[which(baseP\$V1=="C")]
	Gpercent<-baseP\$V2[which(baseP\$V1=="G")]
	GpercentCount<-baseP\$V3[which(baseP\$V1=="G")]
	#Npercent<-baseP\$V2[which(baseP\$V1=="N")]
	#NpercentCount<-baseP\$V3[which(baseP\$V1=="N")]
	GCpercent<-baseP\$V2[which(baseP\$V1=="GC")]
	GCpercentCount<-baseP\$V3[which(baseP\$V1=="GC")]
	aAvg<-sum(Apercent * ApercentCount)/sum(ApercentCount)
	aStd<-sqrt(sum(((Apercent - aAvg)**2)*ApercentCount)/sum(ApercentCount))
	tAvg<-sum(Tpercent * TpercentCount)/sum(TpercentCount)
	tStd<-sqrt(sum(((Tpercent - tAvg)**2)*TpercentCount)/sum(TpercentCount))
	cAvg<-sum(Cpercent * CpercentCount)/sum(CpercentCount)
	cStd<-sqrt(sum(((Cpercent - cAvg)**2)*CpercentCount)/sum(CpercentCount))
	gAvg<-sum(Gpercent * GpercentCount)/sum(GpercentCount)
	gStd<-sqrt(sum(((Gpercent - gAvg)**2)*GpercentCount)/sum(GpercentCount))
	#nAvg<-sum(Npercent * NpercentCount)/sum(NpercentCount)
	#nStd<-sqrt(sum(((Npercent - nAvg)**2)*NpercentCount)/sum(NpercentCount))
	gcAvg<-sum(GCpercent * GCpercentCount)/sum(GCpercentCount)
	gcStd<-sqrt(sum(((GCpercent - gcAvg)**2)*GCpercentCount)/sum(GCpercentCount))
	GCaggregate<-tapply(GCpercentCount,list(cut(GCpercent,breaks=c(seq(0,100,1)))),FUN=sum)
	Aaggregate<-tapply(ApercentCount,list(cut(Apercent,breaks=c(seq(0,100,1)))),FUN=sum)
	Taggregate<-tapply(TpercentCount,list(cut(Tpercent,breaks=c(seq(0,100,1)))),FUN=sum)
	Caggregate<-tapply(CpercentCount,list(cut(Cpercent,breaks=c(seq(0,100,1)))),FUN=sum)
	Gaggregate<-tapply(GpercentCount,list(cut(Gpercent,breaks=c(seq(0,100,1)))),FUN=sum)

	par(fig=c(fig_x_start,(fig_x_end-fig_x_start)*0.75+fig_x_start,0,1),mar=c(5,6,4,2),xpd=FALSE,cex.main=1.2,new=new)
	plot(GCaggregate/1000000,xlim=c(0,100),type="h",lwd=4,xlab=paste(xlab,"GC (%)"),ylab=ylab,lend=2)
	legend.txt<-c(paste("GC",sprintf ("%.2f%%",gcAvg),"±",sprintf ("%.2f",gcStd)))
	legend('topright',legend.txt,bty='n')
	if (totalReads< $trimmed_num){
            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
        }

	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.75,1), mar=c(3, 2, 2, 2),new=TRUE,cex.main=1)
	legend.txt<-c(paste("A",sprintf ("%.2f%%",aAvg),"±",sprintf ("%.2f",aStd)))
	plot(Aaggregate/1000000,xlim=c(0,50),type="h",lwd=2,main=legend.txt,xlab="",ylab="",)

	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.5,0.75),mar=c(3, 2, 2, 2),new=TRUE)
	legend.txt<-c(paste("T",sprintf ("%.2f%%",tAvg),"±",sprintf ("%.2f",tStd)))
	plot(Taggregate/1000000,xlim=c(0,50),type="h",lwd=2,main=legend.txt,xlab="",ylab="")

	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0.25,0.5),mar=c(3, 2, 2, 2),new=TRUE)
	legend.txt<-c(paste("C",sprintf ("%.2f%%",cAvg),"±",sprintf ("%.2f",cStd)))
	plot(Caggregate/1000000,xlim=c(0,50),type="h",lwd=2,main=legend.txt,xlab="",ylab="")

	par(fig=c((fig_x_end-fig_x_start)*0.75+fig_x_start,fig_x_end,0,0.25),mar=c(3, 2, 2, 2),new=TRUE)
	legend.txt<-c(paste("G",sprintf ("%.2f%%",gAvg),"±",sprintf ("%.2f",gStd)))
	plot(Gaggregate/1000000,xlim=c(0,50),type="h",lwd=2,main=legend.txt,xlab="",ylab="")
}
if(file.exists(\"$qa_nuc_composition_file\")){
    readGC_plot(\"$qa_nuc_composition_file\",qa.readsCount,0,0.5,"Input Reads","Number of reads (millions)",FALSE)
    readGC_plot(\"$nuc_composition_file\",readsCount,0.5,1,"Trimmed Reads","",TRUE)
    #abline(v=0.5,lty=2,xpd=TRUE)
} else {
    readGC_plot(\"$nuc_composition_file\",readsCount,0,1,"","Number of reads (millions)",FALSE)
}
par(def.par)#- reset to default
title("Reads GC content",adj=0)


#ATCG composition per base ATCG plot
baseM<-read.table(file=\"$base_matrix\")
ATCG_composition_plot <- function(base_matrix_file,totalReads,xlab,ylab,xlab_adj) {
	baseM<-read.table(file=base_matrix_file)
	aBase<-baseM\$V1
	tBase<-baseM\$V2
	cBase<-baseM\$V3
	gBase<-baseM\$V4
	nBase<-baseM\$V5

	aPer<-(aBase/rowSums(baseM))*100
	tPer<-(tBase/rowSums(baseM))*100
	cPer<-(cBase/rowSums(baseM))*100
	gPer<-(gBase/rowSums(baseM))*100

        ymax<-floor(max(aPer,tPer,cPer,gPer))
        ymin<-floor(min(aPer,tPer,cPer,gPer))
        if((ymin-5)>0){ymin <- ymin-5}else{ymin<-0}
	xpos<-seq(1,length(aBase),1)
	plot(xpos,aPer,col='green3',type='l',xaxt='n',xlab=xlab,ylab=ylab ,ylim=c(ymin,ymax+5))
	lines(xpos,tPer,col='red')
	lines(xpos,cPer,col='blue')
	lines(xpos,gPer,col='black')
	axis(1,at=xpos,labels=xpos+xlab_adj)
	legend('topright',c('A','T','C','G'),col=c('green3','red','blue','black'),box.col=0,lwd=1)
	if (totalReads< $trimmed_num){
            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
        }
	return(nBase)
}
if(file.exists(\"$qa_base_matrix\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    qa.nBase<-ATCG_composition_plot(\"$qa_base_matrix\",qa.readsCount,"Input Reads Base",'Base content (%)',0)
    nBase<-ATCG_composition_plot(\"$base_matrix\",readsCount,"Trimmed Reads Base","",$trim_5_end)
}else{
    nBase<-ATCG_composition_plot(\"$base_matrix\",readsCount,"",'Base content (%)',0)
}
par(def.par)#- reset to default
title("Nucleotide Content Per Cycle")


#N composition per Base plot
N_composition_plot<-function(BaseArray,totalReads,xlab,ylab,xlab_adj){
  xpos<-seq(1,length(BaseArray),1)
  plot(xpos,BaseArray/totalReads*1000000,col='red',type='l',xaxt='n',xlab=xlab,ylab=ylab,ylim=c(0,max(BaseArray/totalReads*1000000)))
  axis(1,at=xpos,labels=xpos+xlab_adj)
  legend('topright',paste("Total bases: ",sum(BaseArray)),bty='n')
  if (totalReads< $trimmed_num){
      mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
  }
}
if (sum(nBase) >0){
  if(file.exists(\"$qa_base_matrix\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    N_composition_plot(qa.nBase,qa.readsCount,"Input reads Position","N Base count per million reads",0)
    N_composition_plot(nBase,readsCount,"Trimmed reads Position","",$trim_5_end)
  }else{
    N_composition_plot(nBase,readsCount,"Position","N Base count per million reads",0)
  } 
}
par(def.par)#- reset to default
title("N Nucleotide Content Per Cycle")

if(file.exists(\"$kmer_rarefaction_file\")){
par(mar=c(5,6,4,2))
kmerfile<-read.table(file=\"$kmer_rarefaction_file\")
sampling_size<-sum(kmerfile\$V1)
sampling<-""
if(sampling_size< $trimmed_num){
    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")
}
cumSeqNum<-cumsum(kmerfile\$V1);
plot(cumSeqNum/1000000,kmerfile\$V3/1000000,xlab=paste(\"Number of Sequence (million)\",sampling), ylab=\"Number of Distinct K-mer (million,k=$kmer)\",type=\'l\',lty=2)
lines(cumSeqNum/1000000,kmerfile\$V2/1000000,col='blue',lwd=2)
title(\"Kmer Rarefaction Curve\")
y<-kmerfile\$V2/1000000
x<-cumSeqNum/1000000
lres<-lm(y~x)
# y=ax+b
a<-format(coef(lres)[[2]], digits = 2)
b<-format(coef(lres)[[1]], digits = 2)
par(def.par)
}

if(file.exists(\"$kmer_histogram_file\")){
par(mar=c(5,6,4,2))
sampling<-""
if(sampling_size<$trimmed_num){
    sampling<-paste(\"(Sampling\",format(sampling_size/1000000,digit=3),\"M Reads)\")
}
kmerHfile<-read.table(file=\"$kmer_histogram_file\")
barplot(kmerHfile\$V2[3:length(kmerHfile\$V1)],names.arg=kmerHfile\$V1[3:length(kmerHfile\$V1)],xlab=\"K-mer Count (k=$kmer)\",log=\'y\',ylab=\"Number of K-mer\",main=paste(\"K-mer Frequency Histogram\",sampling),col=\'darkcyan\',border=\'darkcyan\')
total_kmer<-sum(kmerHfile\$V2)
legend('topleft',paste(\"Total: \",total_kmer),bty='n')
par(fig=c(0.5,0.9,0.5,1), new=TRUE)
barplot(kmerHfile\$V2[3:length(kmerHfile\$V1)],xlim=c(0,100),names.arg=kmerHfile\$V1[3:length(kmerHfile\$V1)],log=\'y\',col=\'darkcyan\',border=\'darkcyan\')
par(fig=c(0,1,0,1))
par(def.par)#- reset to default

}

# read avg quality count barplot 
quality_histogram<-function(qual_histogram_file,totalReads,xlab,ylab){
	Qhist_file<-read.table(file=qual_histogram_file,header=TRUE)
	cumulate<-cumsum(Qhist_file\$readsNum)
	par(mar=c(5,6,5,4))
	if (missing(ylab)){ylab2<-""} else {ylab2<-ylab}
	plot(Qhist_file\$Score,Qhist_file\$readsNum/1000000,type='h',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),xlab=xlab, ylab=ylab2,lwd=12,lend=2)
	par(new=TRUE)
	plot(Qhist_file\$Score,cumulate/sum(Qhist_file\$readsNum)*100,type='l',xlim=c(max(Qhist_file\$Score),min(Qhist_file\$Score)),yaxt='n',xaxt='n',ylab="",xlab="",col='blue',lwd=3)
	axis(4,col='blue',col.ticks='blue',col.axis='blue')
	if (missing(ylab)){
	  mtext(side=4,'Cumulative Percentage',line=2,col='blue')
	}
	Qover20Reads<-sum(as.numeric(Qhist_file\$readsNum[Qhist_file\$Score>=20]))
	Qover20ReadsPer<-sprintf("%.2f%%",Qover20Reads/sum(Qhist_file\$readsNum)*100)
	Qover20Bases<-sum(as.numeric(Qhist_file\$readsBases[Qhist_file\$Score>=20]))
	Qover20AvgLen<-sprintf("%.2f",Qover20Bases/Qover20Reads)
	TotalBases<-sum(as.numeric(Qhist_file\$readsBases))
	mtext(side=3,paste("Number of Q>=20 reads:",formatC(Qover20Reads,format="d",big.mark=","),"(",Qover20ReadsPer,")",", mean Length:",Qover20AvgLen),adj=0,cex=0.8,line=0.3)
	if (totalReads< $trimmed_num){
           mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),line=1,adj=0,cex=0.8)
        }
    return(TotalBases)
}
if(file.exists(\"$qa_avg_quality_histogram\"))
{
    par(mfrow=c(1,2))
    qa.totalBases<-quality_histogram(\"$qa_avg_quality_histogram\",qa.readsCount,"Input Reads Avg Score","Reads Number (millions)")
    totalBases<-quality_histogram(\"$avg_quality_histogram\",readsCount,"Trimmed Reads Avg Score")
}else{
    totalBases<-quality_histogram(\"$avg_quality_histogram\",readsCount,"Avg Score","Reads Number (millions)")
}
par(def.par)#- reset to default
title('Reads Average Quality Histogram')

# read in matrix file for the following three plots
quality_boxplot<-function(quality_matrix_file,totalReads,totalBases,xlab,ylab,xlab_adj){
	z<-as.matrix(read.table(file=quality_matrix_file))
	x<-1:nrow(z)
	y<-1:ncol(z)
	y<-y-1

	#quality boxplot per base
	is.wholenumber <- function(x, tol = .Machine\$double.eps^0.5)  abs(x - round(x)) < tol
	plot(1:length(x),x,type='n',xlab=xlab,ylab=ylab, ylim=c(0,max(y)+1),xaxt='n')
	axis(1,at=x,labels=x+xlab_adj)

	for (i in 1:length(x)) {
	  total<-sum(z[i,])
	  qAvg<-sum(y*z[i,])/total
	  if (is.wholenumber(total/2))
	  {
		 med<-( min(y[cumsum((z[i,]))>=total/2]) + min(y[cumsum((z[i,]))>=total/2+1]) )/2
	  }
	  else
	  {
		 med<-min(y[cumsum((z[i,]))>=ceiling(total/2)])
	  }

	  if (is.wholenumber(total/4))
	  {
		 Q1<-( min(y[cumsum((z[i,]))>=total/4]) + min(y[cumsum((z[i,]))>=total/4+1]) )/2
	  }
	  else
	  {
		 Q1<-min(y[cumsum((z[i,]))>=round(total/4)])
	  }

	  if (is.wholenumber(total/4*3))
	  {
		 Q3<-( min(y[cumsum((z[i,]))>=total/4*3]) + min(y[cumsum((z[i,]))>=total/4*3+1]) )/2
	  }
	  else
	  {
		 Q3<-min(y[cumsum((z[i,]))>=round(total/4*3)])
	  }
	  maxi<-max(y[z[i,]>0])
	  mini<-min(y[z[i,]>0])
	  #if (Q1 == 'Inf') {Q1 = maxi}
	  if (Q3 == 'Inf') {Q3 = maxi}
	  IntQ<-Q3-Q1
	  mini<-max(mini,Q1-1.5*IntQ)
	  maxi<-min(maxi,Q3+1.5*IntQ)
	  rect(i-0.4,Q1,i+0.4,Q3,col='bisque')
	  lines(c(i,i),c(Q3,maxi),lty=2)
	  lines(c(i,i),c(mini,Q1),lty=2)
	  lines(c(i-0.4,i+0.4),c(mini,mini))
	  lines(c(i-0.4,i+0.4),c(maxi,maxi))
	  lines(c(i-0.4,i+0.4),c(med,med))
	  #points(i,qAvg,col='red')
	  reads_num<-formatC(totalReads,format="d",big.mark=",")
	  reads_base<-formatC(totalBases,format="d",big.mark=",")
	  abline(h=20, col = "gray60")
	  legend("bottomleft",c(paste("# Reads: ",reads_num),paste("# Bases:",reads_base)))
	  if (totalReads< $trimmed_num){
              mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
          }
	## for outliers
	#points()
	}
}
if (file.exists(\"$qa_quality_matrix\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    quality_boxplot(\"$qa_quality_matrix\",qa.readsCount,qa.totalBases,"Input Reads Position","Quality score",0)
    quality_boxplot(\"$quality_matrix\",readsCount,totalBases,"Trimmed Reads Position","",$trim_5_end)
}else{
    quality_boxplot(\"$quality_matrix\",readsCount,totalBases,"Position","Quality score",0)
}
par(def.par)#- reset to default
title("Quality Boxplot Per Cycle")

#quality 3D plot
quality_3d_plot<-function(quality_matrix_file,totalReads,xlab,ylab){
	z<-as.matrix(read.table(file=quality_matrix_file))
	x<-1:nrow(z)
	y<-1:ncol(z)
	y<-y-1
    persp(x,y,z/1000000,theta = 50, phi = 30, expand = 0.7, col = rev(terrain.colors(length(z),alpha=0.8)),border=NA,ntick=10,ticktype="detailed",xlab=xlab,ylab=ylab,zlab="",r=6,shade=0.75)
    mtext(side=2, "Frequency (millions)",line=2)
	if (totalReads< $trimmed_num){
            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),adj=0)
        }
}
if (file.exists(\"$qa_quality_matrix\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    quality_3d_plot(\"$qa_quality_matrix\",qa.readsCount,"Input Reads Position","Q Score")
    quality_3d_plot(\"$quality_matrix\",readsCount,"Trimmed Reads Position","Q Score")
}else{
    quality_3d_plot(\"$quality_matrix\",readsCount,"Position","Q Score")
}
par(def.par)#- reset to default
title("Quality 3D plot. (Position vs. Score vs. Frequency)")

#Quality count bar plot
upper_limit<-41
quality_count_histogram<-function(quality_matrix_file,totalReads,highestScore,xlab,ylab){
    z<-as.matrix(read.table(file=quality_matrix_file));
    col<-colSums(z)
    less30columnNum<-length(col)-highestScore+30-1
    atleast30columnNum<-highestScore-30+1
    color<-c(rep('blue',less30columnNum),rep('darkgreen',atleast30columnNum))
    over30per<-sprintf("%.2f%%",sum(col[(less30columnNum+1):length(col)])/sum(col)*100)
    countInM<-col/1000000
    avgQ<-sprintf("%.2f",sum(seq(0,41,1)*col)/sum(col))
    plot(seq(0,highestScore,1),countInM,col=color,type='h',ylab=ylab,xlab=xlab,lwd=12,lend=2,bty='n')
    abline(v=29.5,col='darkgreen')
    text(30,(max(countInM)-min(countInM))*0.9,labels=">=Q30",cex=0.8,adj=0,col='darkgreen')
    text(30,(max(countInM)-min(countInM))*0.85,labels=over30per,cex=0.8,adj=0,col='darkgreen')
    mtext(side=3,paste("Average: ",avgQ),adj=0)
	if (totalReads< $trimmed_num){
            mtext(side=3,paste(\"(Sampling\",formatC(totalReads/1000000,digits=3),\"M Reads)\"),line=1 ,adj=0)
        }
}
if (file.exists(\"$qa_quality_matrix\")){
    par(mfrow=c(1,2),mar=c(5,6,4,2))
    quality_count_histogram(\"$qa_quality_matrix\",qa.readsCount,upper_limit,"Input Reads Q score","Total (million)")
    quality_count_histogram(\"$quality_matrix\",readsCount,upper_limit,"Trimmed Reads Q score","")
}else{
    quality_count_histogram(\"$quality_matrix\",readsCount,upper_limit,"Q score","Total (million)")
}
par(def.par)#- reset to default
title("Quality report")

tmp<-dev.off()

quit()

RSCRIPT

    close R;
    system ("R --vanilla --silent --quiet < $outDir/tmp$$.R 1>/dev/null");
    unless ($debug){system ("rm $outDir/tmp$$.*");}
    system ("rm $outDir/Rplots.pdf") if (-e "$outDir/Rplots.pdf");
    system ("rm Rplots.pdf") if (-e "Rplots.pdf");
}

sub print_quality_report_files{
    my $quality_matrix=shift;
    my $avg_quality_histogram=shift;
    my $base_matrix=shift;
    my $nuc_composition_file=shift;
    my $length_histogram=shift;
    my $position_r=shift;
    my $AverageQ_r=shift;
    my $base_position_r=shift;
    my $base_content_r=shift;
    my $len_hash_r=shift;
    open (OUT, ">$quality_matrix");
    open (BASE, ">$base_matrix");
    foreach my $pos2(sort {$a<=>$b} keys %{$position_r}){
        my $q_string;
        my $n_string;
        for ($lowest_illumina_score..$highest_illumina_score)
        {
            if ($position_r->{$pos2}->{$_}){
              $q_string .=  $position_r->{$pos2}->{$_}."\t";
            }
            else
            {
              $q_string .= "0\t";
            }
        }
        for my $base ("A","T","C","G","N")
        { 
            if ($base_position_r->{$pos2}->{$base}){
              $n_string .= $base_position_r->{$pos2}->{$base}."\t";
            }
            else
            {
              $n_string .= "0\t";
            }
        }
        $q_string =~ s/\t$/\n/;
        $n_string =~ s/\t$/\n/;
        print OUT $q_string;
        print BASE $n_string;
    }
    close OUT;
    close BASE;

    open (OUT2, ">$avg_quality_histogram");
    my $Key=$highest_illumina_score;
    my $h_print_string;
    while ($Key >= 0)
    {
        if(defined $AverageQ_r->{$Key}->{reads}){
           $h_print_string .= "$Key\t".
           $AverageQ_r->{$Key}->{reads}. "\t".$AverageQ_r->{$Key}->{bases}."\n";
        }else{
           $h_print_string .= "$Key\t".
              "0\t".
              "0\n";
        }
        --$Key;
    }
    print OUT2 "Score\treadsNum\treadsBases\n";
    print OUT2 $h_print_string; 
    close OUT2;

    open (OUT3, ">$nuc_composition_file");
    my %base_content=%{$base_content_r};
    for my $nuc ("A","T","C","G","N","GC")
    {
      foreach my $key ( sort {$a<=>$b} keys %{$base_content{$nuc}})
      {
             print OUT3 "$nuc\t$key\t${$base_content{$nuc}}{$key}\n";
      }
    } 
    close OUT3;
 
    open (OUT4,">$length_histogram");
	my @len_list= sort {$a<=>$b} keys %{$len_hash_r};
    for my $key (1..$len_list[-1])
    {
	if ($len_hash_r->{$key}) 
        {
            print OUT4 "$key\t$len_hash_r->{$key}\n";
	}
	else
	{
	    print OUT4 "$key\t0\n";
	}
    }
    close OUT4;
}

sub is_file_empty 
{
    #check file exist and non zero size
    my $file=shift;
    my $empty=1;
    if (-e $file) {$empty=0};
    if (-z $file) {$empty=1};
    return $empty;
}

sub qc_process {
  my ($input, $input2,$random_select_file_flag) = @_;
  my ($h1,$s,$s_trimmed,$h2,$q, $q_trimmed); my ($len, $len_q, $trim_len)=(0,0,0);
  my ($r2_h1,$r2_s,$r2_s_trimmed,$r2_h2,$r2_q,$r2_q_trimmed); my ($r2_len, $r2_len_q, $r2_trim_len)=(0,0,0);
  my %stats;
  my %random_sub_stats;
  $stats{filter}->{adapter}->{readsNum}=0;
  my ($seq_r,$random_sub_seq_r);;
  my $avg_q;
  my ($pos5,$pos3);
  my ($raw_seq_num_1,$raw_seq_num_2,$total_raw_seq_len);
  my ($trim_seq_num_1,$trim_seq_num_2,$trim_seq_len, $total_trim_seq_len)=(0,0,0,0);
  my @trim_seq_len=();
  my ($paired_seq_num,$total_paired_bases); 
  my (%tmp1,%tmp2);
  my ($drop_1,$drop_2)=(0,0);
  my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input", qr/\.[^.]*/ );
  my $trim_output_1="$outDir/${i_file_name}_trim.fastq";
  my $trim_output_2;
  my $trim_output_unpaired;
  my $trim_output_discard="$outDir/${i_file_name}_trim_discard.fastq";
  open(IN,"$input") or die "$input\t$!";
  if (! $qc_only)
  {
     open(OUT,"> $trim_output_1");
     open(DISCARD,">$trim_output_discard") if ($output_discard);
  }
  if ($input2) # paired mate
  {
     open(IN2,"$input2") or die "$input2\t$!";
     my ( $i_file_name, $i_path, $i_suffix ) = fileparse( "$input2", qr/\.[^.]*/ );
     $trim_output_2="$outDir/${i_file_name}_trim.fastq";
     $trim_output_unpaired="$outDir/${i_file_name}_trim_unpaired.fastq";
     if (! $qc_only)
     {
         open(OUT2,"> $trim_output_2");
         open(UNPAIR,"> $trim_output_unpaired");
     }
  }

  while ($h1 = <IN>) {  # read first header
        $drop_1=0;
        $raw_seq_num_1++;
        $s = <IN>;  # read sequence
        chomp $s;
	$len = length($s);
        $h2 = <IN>;  # read second header
        $q = <IN>;  # read quality scores
        chomp $q;
        $len_q = length($q);
        $total_raw_seq_len += $len;
        $s_trimmed=$s;
        $q_trimmed=$q;
        $pos5=0;
        $pos3=$len+1;
        $trim_len=$len;
        if ($filter_adapter)
        {
            (my $adapterFlag, $s_trimmed,$pos5,$pos3, my $adapterID)=&filter_adapter($s,$filterAdapterMismatchRate);
            if ($adapterFlag) 
            {
                $trim_len=length($s_trimmed);
                $q_trimmed=($trim_len)?substr($q,$pos5,$pos3-$pos5-1):"";
                $stats{adapter}->{$adapterID}->{readsNum}++;
                $stats{adapter}->{$adapterID}->{basesNum} += ($len - $trim_len);
                $stats{filter}->{adapter}->{readsNum}++;
                $stats{filter}->{adapter}->{basesNum} += ($len - $trim_len);
            }
        }
        if ( $len != $len_q){
            $stats{filter}->{len_ne}->{readsNum}++;
            $stats{filter}->{len_ne}->{basesNum}+=$trim_len;
            warn "$h1 sequence length is no equal to quality string length. It will be filtered.\n";
            $drop_1=1;
        }
        if ($trim_5_end && !$qc_only && ($trim_len>$trim_5_end))
        {
            $s_trimmed=substr($s_trimmed,$trim_5_end);
            $q_trimmed=substr($q_trimmed,$trim_5_end);
            $trim_len = length ($s_trimmed);
            $pos5=$pos5 + $trim_5_end;
            
        }
        if ($trim_3_end && !$qc_only && ($trim_len>$trim_3_end))
        {
            $s_trimmed=substr($s_trimmed,0,$trim_len-$trim_3_end);
            $q_trimmed=substr($q_trimmed,0,$trim_len-$trim_3_end);
            $trim_len = length ($s_trimmed);
            $pos3=$pos3-$trim_3_end;
        }
        #apply length filter
        if ($trim_len < $opt_min_L || $trim_len == 0 || $trim_len<=($trim_5_end+$trim_3_end) )
        {
            $stats{filter}->{len}->{readsNum}++;
            $stats{filter}->{len}->{basesNum}+=$trim_len;
            $drop_1=1;
        }
        if ($qc_only)
        {
            $drop_1 =0;
            $s_trimmed=$s;
            $q_trimmed=$q;
            $pos5=0;
            $pos3=$len+1;
            $trim_len=$len;
        } 
        else # do quality trim
        {
            if ($drop_1==0){
               my $before_trim_len=$trim_len;
               if ($mode =~ /hard/i)
               {
                   ($s_trimmed,$q_trimmed,$pos5,$pos3)= &hard_trim ($trim_len,$s,$q,$pos5,$pos3);
               }
               elsif ($mode =~ /BWA_plus/i)
               {
                   ($s_trimmed,$q_trimmed,$pos5,$pos3)= &bwa_trim_plus ($trim_len,$s,$q,$pos5,$pos3);
               }
               elsif ($mode =~ /BWA/i)
               {
                   ($s_trimmed,$q_trimmed,$pos5,$pos3)= &bwa_trim ($trim_len,$s,$q,$pos5,$pos3);
               }
               $trim_len=length($q_trimmed);
               $stats{filter}->{qualTrim}->{basesNum}+= ($before_trim_len - $trim_len);
               $stats{filter}->{qualTrim}->{readsNum}++ if ( ($before_trim_len - $trim_len)>0);
               # apply length filter, trimming, "N" filter
               if ($trim_len < $opt_min_L || $trim_len == 0 )
               {
                   $stats{filter}->{len}->{readsNum}++;
                   $stats{filter}->{len}->{basesNum}+=$trim_len;
                   $drop_1=1;
               }
            }
        }
        # "N" filter
        if ($s_trimmed =~ /N{$N_num_cutoff,}/i and $drop_1==0)
        {
            $stats{filter}->{NN}->{readsNum}++;
            $stats{filter}->{NN}->{basesNum}+=$trim_len;
            $drop_1=1;
            $drop_1=0 if ($qc_only);
        }
        if ($filter_phiX and $drop_1==0)
        {
           my $phiX_Flag=&filter_phiX($s_trimmed,$filterAdapterMismatchRate);
           if ($phiX_Flag)
           {
               $stats{filter}->{phiX}->{readsNum}++; 
               $stats{filter}->{phiX}->{basesNum}+=$trim_len; 
               $drop_1=1;
           }
           $drop_1=0 if ($qc_only);
        }
        if ($random_select_file_flag and !$qc_only and !$trim_only)
        {
            my ($tmp1,$tmp2);
            ($random_sub_seq_r,$tmp1,$tmp2)=&get_base_and_quality_info($s,$q,$len,0,$len+1,\%random_sub_stats,1);
            %random_sub_stats=%{$random_sub_seq_r};
        }
        if ($drop_1==0 and ! $trim_only)
        {
          #     print "drop1\t",$drop_1,"\t";
            ($seq_r,$drop_1,$s_trimmed)=&get_base_and_quality_info($s_trimmed,$q_trimmed,$trim_len,$pos5,$pos3,\%stats,$qc_only);
            %stats=%{$seq_r};
         #   print $drop_1,"\n";
        }
        if ($drop_1==0){  # pass all filters...
            $q_trimmed=&quality_encoding_coversion($q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
            $trim_seq_num_1++;
            push @trim_seq_len, $trim_len;
            $total_trim_seq_len += $trim_len;
        }
        
        if ($input2) {
           $drop_2=0;
           $r2_h1 = <IN2>;
           $raw_seq_num_2++;
           $r2_s = <IN2>;  # mate read sequence
           chomp $r2_s;
           $r2_len = length($r2_s);
           $r2_h2 = <IN2>;  # mate read second header
           $r2_q = <IN2>;  # mate read quality scores
           chomp $r2_q;
           $r2_len_q = length($r2_q);
           $total_raw_seq_len += $r2_len;
           $r2_s_trimmed=$r2_s;
           $r2_q_trimmed=$r2_q;
           $pos5=0;
           $pos3=$r2_len+1;
           $r2_trim_len=$r2_len;

           if ($filter_adapter)
           {
               $drop_2=&filter_adapter($r2_s,$filterAdapterMismatchRate);
               (my $adapterFlag, $r2_s_trimmed,$pos5,$pos3, my $adapterID)=&filter_adapter($r2_s,$filterAdapterMismatchRate);
               if ($adapterFlag) 
               {
                   $r2_trim_len=length($r2_s_trimmed);
                   $r2_q_trimmed=($r2_trim_len)?substr($r2_q,$pos5,$pos3-$pos5-1):"";
                   $stats{adapter}->{$adapterID}->{readsNum}++;
                   $stats{adapter}->{$adapterID}->{basesNum} += ($r2_len - $r2_trim_len);
                   $stats{filter}->{adapter}->{readsNum}++;
                   $stats{filter}->{adapter}->{basesNum} += ($r2_len - $r2_trim_len);
               }
           }
           if ( $r2_len != $r2_len_q)
           {
               $stats{filter}->{len_ne}->{readsNum}++;
               $stats{filter}->{len_ne}->{basesNum}+=$r2_trim_len;
               warn "$r2_h1 sequence length is no equal to quality string length. It will be filtered.\n";
               $drop_2=1;
           }
           if ($trim_5_end && !$qc_only && ($trim_5_end< $r2_trim_len))
           {
               $r2_s_trimmed=substr($r2_s_trimmed,$trim_5_end);
               $r2_q_trimmed=substr($r2_q_trimmed,$trim_5_end);
               $r2_trim_len = length ($r2_s_trimmed);
               $pos5=$pos5+$trim_5_end;
           }
           if ($trim_3_end && !$qc_only && ($trim_3_end< $r2_trim_len))
           {
               $r2_s_trimmed=substr($r2_s_trimmed,0,$r2_trim_len-$trim_3_end);
               $r2_q_trimmed=substr($r2_q_trimmed,0,$r2_trim_len-$trim_3_end);
               $r2_trim_len = length ($r2_s_trimmed);
               $pos3=$pos3 - $trim_3_end;
           }
           if ($r2_trim_len < $opt_min_L || $r2_trim_len == 0 || $r2_trim_len <= ($trim_5_end + $trim_3_end) )
           {
               $stats{filter}->{len}->{readsNum}++;
               $stats{filter}->{len}->{basesNum}+=$r2_trim_len;
               $drop_2=1;
           }
           if ($qc_only)
           {
               $drop_2 =0;
               $r2_s_trimmed=$r2_s;
               $r2_q_trimmed=$r2_q;
               $pos5=0;
               $pos3=$r2_len+1;
               $r2_trim_len=$r2_len;
           } 
           else  # do quality trim
           {
               if ($drop_2==0){
                   my $before_trim_len=$r2_trim_len;
                   if ($mode =~ /hard/i)
                   {
                       ($r2_s_trimmed, $r2_q_trimmed, $pos5, $pos3)=&hard_trim ($r2_trim_len,$r2_s,$r2_q,$pos5,$pos3);
                   }
                   elsif ($mode =~ /BWA_plus/i)
                   {
                       ($r2_s_trimmed, $r2_q_trimmed, $pos5, $pos3)=&bwa_trim_plus ($r2_trim_len,$r2_s,$r2_q,$pos5,$pos3);
                   }
                   elsif ($mode =~ /BWA/i)
                   {
                       ($r2_s_trimmed, $r2_q_trimmed, $pos5, $pos3)=&bwa_trim ($r2_trim_len,$r2_s,$r2_q,$pos5,$pos3);
                   }
                   $r2_trim_len=length($r2_q_trimmed);
                   $stats{filter}->{qualTrim}->{basesNum}+= ($before_trim_len - $r2_trim_len);
                   $stats{filter}->{qualTrim}->{readsNum}++ if ( ($before_trim_len - $r2_trim_len)>0);
                   if ($r2_trim_len < $opt_min_L || $r2_trim_len == 0 )
                   {
                       $stats{filter}->{len}->{readsNum}++;
                       $stats{filter}->{len}->{basesNum} += $r2_trim_len;
                       $drop_2=1;
                   }
               }
           }
           if ($r2_s_trimmed =~ /N{$N_num_cutoff,}/i and $drop_2==0)
           {
               $stats{filter}->{NN}->{readsNum}++;
	       $stats{filter}->{NN}->{basesNum}+=$r2_trim_len;
               $drop_2=1;
               $drop_2=0 if ($qc_only);
           }
           if ($filter_phiX and $drop_2==0)
           {
               my $phiX_Flag=&filter_phiX($r2_s_trimmed,$filterAdapterMismatchRate);
               if ($phiX_Flag)
               {
                   $stats{filter}->{phiX}->{readsNum}++; 
                   $stats{filter}->{phiX}->{basesNum}+=$r2_trim_len; 
                   $drop_2=1;
               }
               $drop_2=0 if ($qc_only);
           }
           if ($random_select_file_flag and !$qc_only and !$trim_only)
           {
               my ($tmp1,$tmp2);
               ($random_sub_seq_r,$tmp1,$tmp2)=&get_base_and_quality_info($r2_s,$r2_q,$r2_len,0,$r2_len+1,\%random_sub_stats,1);
               %random_sub_stats=%{$random_sub_seq_r};
           }
           if ($drop_2==0 and !$trim_only)
           {
           #    print "drop2\t",$drop_2,"\t";
               ($seq_r,$drop_2,$r2_s_trimmed)=&get_base_and_quality_info($r2_s_trimmed,$r2_q_trimmed,$r2_trim_len,$pos5,$pos3,\%stats,$qc_only);
               %stats=%{$seq_r};
            #   print $drop_2,"\n";
           }
          
           if ($drop_2==0){ # pass all filters
               $r2_q_trimmed=&quality_encoding_coversion($r2_q_trimmed,$ascii,$out_offset) if ($ascii != $out_offset);
               $trim_seq_num_2++;
               push @trim_seq_len, $r2_trim_len;
               $total_trim_seq_len += $r2_trim_len;
           }
        }
        if ($drop_1==0 and $drop_2==0)
        {
            $paired_seq_num = $paired_seq_num + 2 if ($input2);
            $total_paired_bases += $trim_len + $r2_trim_len if ($input2);
        }
 
        # output trimmed files
        if (! $qc_only)
        {
            if ($drop_1==0 and $drop_2==0)
            {
	        print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                print OUT2 $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==1 and $drop_2==0)
            {
                print DISCARD $h1.$s."\n".$h2.$q."\n" if ($output_discard);
                print UNPAIR $r2_h1.$r2_s_trimmed."\n".$r2_h2.$r2_q_trimmed."\n";
            }
            elsif ($drop_1==0 and $drop_2==1)
            {
                if ($input2)
                {
                   print UNPAIR $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n" if ($output_discard);
                }
                else
                {
                   print OUT $h1.$s_trimmed."\n".$h2.$q_trimmed."\n";
                }
            }
            elsif ($output_discard)
            {   
                print DISCARD $h1.$s."\n".$h2.$q."\n";
                if ($input2)
                {
                   print DISCARD $r2_h1.$r2_s."\n".$r2_h2.$r2_q."\n";
                }
            }
        }
  } # end while

  my ($trim_seq_len_std,$trim_seq_len_avg,$max,$min,$median)=&standard_deviation(@trim_seq_len);
  
  close(IN);
  close(OUT);
  close IN2 if ($input2);
  close OUT2 if ($input2);
  close UNPAIR if ($input2);
  close DISCARD;
  $stats{raw_seq_num}=$raw_seq_num_1+$raw_seq_num_2;
  $stats{raw_seq_num_1}=$raw_seq_num_1;
  $stats{raw_seq_num_2}=$raw_seq_num_2;
  $stats{trim_seq_num}=$trim_seq_num_1+$trim_seq_num_2;
  $stats{trim_seq_num_1}=$trim_seq_num_1;
  $stats{trim_seq_num_2}=$trim_seq_num_2;
  $stats{total_raw_seq_len}=$total_raw_seq_len;
  $stats{total_trim_seq_len}=$total_trim_seq_len;
  $stats{trim_seq_len_std}=$trim_seq_len_std;
  $stats{trim_seq_len_avg}=$trim_seq_len_avg;
  $stats{max}=$max;
  $stats{min}=$min;
  $stats{median}=$median;
  $stats{paired_seq_num}=$paired_seq_num;
  $stats{total_paired_bases}=$total_paired_bases;
  $stats{file_name_1}=$input;
  $stats{file_name_2}=$input2;
  $stats{trim_file_name_1}=$trim_output_1;
  $stats{trim_file_name_2}=$trim_output_2;
  $stats{random_sub_raw}=$random_sub_seq_r;
  return \%stats;
}

sub filters  # not used
{
    my $s = shift;
    my $q = shift;
    my $drop = 0;
    # "N" base filter
    my $N_num = $s=~ tr/Nn/Nn/;
    if ( $N_num > $N_num_cutoff and $drop==0) {$drop=1;}
    #low complexity filter
    if (&low_complexity_filter($s,$low_complexity_cutoff_ratio) and $drop==0) {$drop=1;}
    
    return $drop;
} 
     
sub bwa_trim
{
        # original bwa trim implementation  at 3' end
        my ($len,$s,$q,$final_pos_5,$final_pos_3) = @_;
        # trim 3' end
        my $pos_3 = $final_pos_3-1;
 #       my $final_pos_3 =$pos_3+1;  
        my $area=0;  
        my $maxArea=0;
        
	while ($pos_3>0 and $area>=0) {
		$area += $opt_q - (ord(substr($q,$pos_3-1,1))-$ascii);
		if ($area > $maxArea) {
			$maxArea = $area;
			$final_pos_3 = $pos_3;
		}
		$pos_3--;
	}
        # trim 5' end
  #      my $pos_5=1; 
  #      my $final_pos_5=0;
  #      $maxArea=0;
  #      $area=0;
#	while ($pos_5<$final_pos_3 and $area>=0) {
#		$area += $opt_q - (ord(substr($q,$pos_5-1,1))-$ascii);
#		if ($area > $maxArea) {
#			$maxArea = $area;
#			$final_pos_5 = $pos_5;
#		}
#		$pos_5++;
#	}
        
        $s = substr($s,$final_pos_5,$final_pos_3-$final_pos_5-1);
        $q = substr($q,$final_pos_5,$final_pos_3-$final_pos_5-1);
        return ($s,$q,$final_pos_5,$final_pos_3);
}

sub hard_trim 
{
    my ($len,$s,$q, $final_pos_5, $final_pos_3) =@_;
    # trim 3' end
    my $pos_3 = $final_pos_3-1;   
    while($pos_3 >0)
    {
        if ($opt_q < (ord(substr($q,$pos_3-1,1))-$ascii) )
        {
            $final_pos_3=$pos_3+1;
            last;
        }
        $pos_3--;
    }
    # trim 5' end
    my $pos_5=$final_pos_5; 
    while($pos_5<$pos_3)
    {
        if ($opt_q < (ord(substr($q,$pos_5,1))-$ascii) )
	{
	    $final_pos_5=$pos_5;
	    last;
	}
        $pos_5++;
    }
    $s = substr($s,$final_pos_5,$final_pos_3-$final_pos_5-1);
    $q = substr($q,$final_pos_5,$final_pos_3-$final_pos_5-1);
    return ($s,$q,$final_pos_5,$final_pos_3);
}

sub bwa_trim_plus
{
    # bwa trim implementation from both 5' and 3' end
    # at least scan 5 bases from both end and 2 more bases after the negative area
    my ($len,$s,$q,$final_pos_5,$final_pos_3) = @_;
    my $at_least_scan=5;
    my $num_after_neg=2;
    $at_least_scan = $len if $len <5;
    $num_after_neg = $len if $len <2;
    # trim 3' end
    my $pos_3 = $final_pos_3-1;
 #   my $final_pos_3 =$pos_3+1;  
    my $area=0;  
    my $maxArea=0;
    while ($at_least_scan) 
    {
        $at_least_scan--;    
        if ($pos_3>$num_after_neg and $area>=0) {$at_least_scan=$num_after_neg;}
	    $area += $opt_q - (ord(substr($q,$pos_3-1,1))-$ascii);
	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_3 = $pos_3;
	    }
	    $pos_3--;
    }
    # trim 5' end
    my $pos_5=$final_pos_5+1; 
  #  my $final_pos_5=0;
    $maxArea=0;
    $area=0;
    $at_least_scan=5;
    $at_least_scan = $len if $len <5;
    while ($at_least_scan) 
    {
        $at_least_scan--; 
        if ($pos_5<($final_pos_3-$num_after_neg) and $area>=0) {$at_least_scan=$num_after_neg;}
	    $area += $opt_q - (ord(substr($q,$pos_5-1,1))-$ascii);

	    if ($area > $maxArea) {
		    $maxArea = $area;
		    $final_pos_5 = $pos_5;
	    }
	    $pos_5++;
    }
    

    if ($final_pos_3<=$final_pos_5)
    {
        # since at least scan 5 bases, we need to avoid negative length on substr
        # return empty and it will be filtered by read length;
        $s = substr($s,0,0);
        $q = substr($q,0,0);
    }
    else
    {
        $s = substr($s,$final_pos_5,$final_pos_3-$final_pos_5-1);
        $q = substr($q,$final_pos_5,$final_pos_3-$final_pos_5-1);
    }
    return ($s,$q,$final_pos_5,$final_pos_3);
}


sub get_base_and_quality_info
{
     my $s=shift;
     my $q=shift;
     my $len=shift;
     my $start_pos=shift;
     my $end_pos=shift;
     my $stats=shift;
     my $qc_only=shift;
     my $total_q=0;
     my $avg_q=0;
     my $drop=0;
     my $low_complexity_flag=0;
     my $new_s=$s;
     my ($a_Base,$t_Base,$c_Base,$g_Base,$n_Base)=(0,0,0,0,0);
     my ($a_Base_percent,$t_Base_percent,$c_Base_percent,$g_Base_percent,$n_Base_percent);
     my $base_percent_string;
     my ($previous_base,$dinucleotide,%di_nucleotide);
     my @N_pos;
     my %seq=%{$stats};
     $start_pos++;
     $end_pos--;
     
     for my $pos($start_pos..$end_pos)
     {
         my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
         #$drop=1 if ($q_digit < $stringent_cutoff);
         $seq{qual}->{$pos}->{$q_digit}++;
         $total_q += $q_digit; 
         my $base=uc(substr($s,$pos-$start_pos,1));
         if ($q_digit < $replace_to_N_q and $is_NextSeq and $base eq "G"){
             substr($new_s,$pos-$start_pos,1,"N");
         }
         $a_Base++ if ($base =~ /A/);  
         $t_Base++ if ($base =~ /T/);  
         $c_Base++ if ($base =~ /C/);  
         $g_Base++ if ($base =~ /G/);  
         if ($base =~ /N/)
         {
             $n_Base++;
             if ($replace_N and !$qc_only)
             {
                my $rand_base = &rand_base;
                substr($new_s,$pos-$start_pos,1,$rand_base);
                $seq{replaceN}->{$rand_base}++;
             }
         }
         if ($previous_base and ($previous_base ne $base))
         {
            $dinucleotide = $previous_base.$base;
            $di_nucleotide{$dinucleotide}++;
         }
         $previous_base=$base;
         $seq{base}->{$pos}->{$base}++;
     }
     
     $avg_q = ($len>0)? int ($total_q/$len):0;
   #  if (!$qc_only)
   #  {
         # apply filters:, average quality, low complexity
         #$drop=1 if ($n_Base>$N_num_cutoff);
         if ($avg_q < $opt_avg_cutoff)
         {
             $seq{filter}->{AvgQ}->{readsNum}++;
             $seq{filter}->{AvgQ}->{basesNum}+=$len;
             $drop=1;
         }
         $low_complexity_flag=1 if ( $len >0 and (($a_Base/$len) > $low_complexity_cutoff_ratio) and $drop==0);
         $low_complexity_flag=1 if ( $len >0 and (($t_Base/$len) > $low_complexity_cutoff_ratio) and $drop==0);
         $low_complexity_flag=1 if ( $len >0 and (($c_Base/$len) > $low_complexity_cutoff_ratio) and $drop==0);
         $low_complexity_flag=1 if ( $len >0 and (($g_Base/$len) > $low_complexity_cutoff_ratio) and $drop==0);
         if ($drop == 0)
         {
           foreach my $di_nuc (keys %di_nucleotide)
           {
             if (($di_nucleotide{$di_nuc}*2/$len) > $low_complexity_cutoff_ratio)
             {
              $low_complexity_flag=1;
              delete $di_nucleotide{$di_nuc};
              last; #end for loop
             }
             delete $di_nucleotide{$di_nuc};
           }
         }
         if ($low_complexity_flag)
         {
             $seq{filter}->{lowComplexity}->{readsNum}++;
             $seq{filter}->{lowComplexity}->{basesNum}+=$len;
             $drop=1;
         }
  #   }
     $drop=0 if ($qc_only);
     if ($drop == 1)  #substract the position matrix by one because the dropped read
     { 
         for my $pos($start_pos..$end_pos)
         {
            my $q_digit=ord(substr($q,$pos-$start_pos,1))-$ascii;
            $seq{qual}->{$pos}->{$q_digit}--;
            my $base=uc(substr($s,$pos-$start_pos,1));
            $seq{base}->{$pos}->{$base}--;
         }
         return (\%seq,$drop,$new_s);
     }
     else
     {
         $a_Base_percent = ($len>0)? sprintf("%.2f",$a_Base/$len*100):0;
         $t_Base_percent = ($len>0)? sprintf("%.2f",$t_Base/$len*100):0;
         $c_Base_percent = ($len>0)? sprintf("%.2f",$c_Base/$len*100):0;
         $g_Base_percent = ($len>0)? sprintf("%.2f",$g_Base/$len*100):0;
         $n_Base_percent = ($len>0)? sprintf("%.2f",$n_Base/$len*100):0;
         $seq{ReadAvgQ}->{$avg_q}->{readsNum}++;
         $seq{ReadAvgQ}->{$avg_q}->{basesNum}+=$len;
         $seq{ReadLen}->{$len}++;
         $seq{Base_content}->{A}->{$a_Base_percent}++;
         $seq{Base_content}->{T}->{$t_Base_percent}++;
         $seq{Base_content}->{C}->{$c_Base_percent}++;
         $seq{Base_content}->{G}->{$g_Base_percent}++;
         $seq{Base_content}->{N}->{$n_Base_percent}++;
         $seq{Base_content}->{GC}->{$c_Base_percent+$g_Base_percent}++;
 
         return (\%seq,$drop,$new_s);
     }
} 

sub rand_base
{
    my %hash = ( 0 => "A", 1 => "T", 2 => "C", 3 => "G");
    return $hash{int(rand(4))};
}

sub standard_deviation {
  my(@numbers) = @_;
  #Prevent division by 0 error in case you get junk data
  return undef unless(scalar(@numbers));
  @numbers= sort {$b<=>$a} @numbers;
  my $max=$numbers[0];
  my $min=$numbers[-1];
  my $median;
  if (scalar(@numbers) % 2) {
      $median = $numbers[int(scalar(@numbers)/2)];
  } else {
      $median = ($numbers[scalar(@numbers)/2] + $numbers[scalar(@numbers)/2 - 1]) / 2;
  }

  # Step 1, find the mean of the numbers
  my $total1 = 0;
  foreach my $num (@numbers) {
  $total1 += $num;
  }
  my $mean1 = $total1 / (scalar @numbers);

  # Step 2, find the mean of the squares of the differences
  # between each number and the mean
  my $total2 = 0;
  foreach my $num (@numbers) {
  $total2 += ($mean1-$num)**2;
  }
  my $mean2 = $total2 / (scalar @numbers);

  # Step 3, standard deviation is the square root of the
  # above mean
  my $std_dev = sqrt($mean2);
  return ($std_dev,$mean1,$max,$min,$median);
}


sub low_complexity_filter {  # not used 
	my ($seq,$cut_off_ratio) = @_;
	my $seq_len = length ($seq);
    my @low_complex_array=("A","T","C","G","AT","AC","AG","TA","TC","TG","CA","CT","CG","GA","GT","GC");
    my $filter=0;
    for my $item (@low_complex_array)
    {
      my $item_len = length ($item);
      my $num_low_complex = $seq =~ s/$item/$item/g;
      if (($num_low_complex*$item_len/$seq_len) >= $cut_off_ratio)
      {
        $filter=1;
        last; #end for loop
      }
    }
    return ($filter);
}

sub random_subsample
{
  my $total_num=shift;
  my $output_num=shift;
  my %get;
  if ($output_num <= 0)
  {
     $get{1}=0;
  }
  else
  { 
    while (1)
    {
      $get{int(rand($total_num))}=1;
      my $num_of_random=scalar (keys %get);
      last if ($num_of_random == $output_num or $num_of_random == $total_num);
    }
  }
  return (\%get);
}
sub quality_encoding_coversion 
{  
    # given quality acsii string, input offset and output offset
    my $q_string=shift;
    my $input_offset=shift;
    my $out_offset=shift;
    $q_string =~ s/[;<=>?]/@/g if ($input_offset == 64); # deal with score < 0
    $q_string=~ s/(\S)/chr(ord($1)-$input_offset+$out_offset)/eg;
    return($q_string);
}

sub open_file
{
    my ($file) = @_;
    my $fh;
    my $pid;
    if ( $file=~/\.gz$/i ) { $pid=open($fh, "gunzip -c $file 2>&1 |") or die ("gunzip -c $file: $!"); }
    else { $pid=open($fh,'<',$file) or die("$file: $!"); }
    return ($fh,$pid);
}

sub is_NextSeq{
    $SIG{'PIPE'}=sub{};
    my $file=shift;
    my ($fh,$pid)= open_file($file);
    my $head=<$fh>;
    close $fh;
    kill 9, $pid; # avoid gunzip broken pipe
    
    $SIG{'PIPE'} = 'DEFAULT';
    ($head =~/^\@NS\d+/)?
        return 1:
        return 0;
}

sub checkQualityFormat {
    # $offset_value=&checkQualityFormat($fastq_file)
    # $offset_value = -1 means not proper fastq format.
    $SIG{'PIPE'} = sub{}; # avoid gunzip broken pipe
    my $fastq_file=shift;
    # open the files
    my ($fh,$pid) = open_file($fastq_file);
    # initiate
    my @line;
    my $l;
    my $number;
    my $offset=33;
    my $line_num=0;
    # go thorugh the file
    my $first_line=<$fh>;
    if ($first_line !~ /^@/) {$offset=-1; return $offset;}
    OUTER:while(<$fh>){
      $line_num++;
      # if it is the line before the quality line
      if($_ =~ /^\+/){
    
       $l = <$fh>; # get the quality line
       chomp $l;
       @line = split(//,$l); # divide in chars
       for(my $i = 0; $i <= $#line; $i++){ # for each char
        $number = ord($line[$i]); # get the number represented by the ascii char
      
        # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
        if($number > 104){
	  $offset=33;  # pacbio CCS reads
          last OUTER;
        }
        elsif($number > 74){ # if solexa/illumina
          $offset=64;
          #die "This file is solexa/illumina format\n"; # print result to terminal and die
          # read a few more quality line to make sure the offset 
          last OUTER if ($line_num > 100);
        }elsif($number < 59){ # if sanger
          $offset=33;
          #die "This file is sanger format\n"; # print result to terminal and die
          last OUTER;
        }
       }
      }
    }
   close $fh;
   kill 9 ,$pid;
   $SIG{'PIPE'} = 'DEFAULT';
   return $offset;
}

sub split_fastq {
   my ($input,$outDir,$subfile_size)=@_;

   my $seqThisFile = 0;
   my $fileCount = 0;
   my $name="";
   my $seq="";
   my $q_name="";
   my $qual_seq="";
   my @subfiles;
   my $seq_num=0;
   my $total_seq_length=0;
   my ($file_name, $i_path, $i_suffix_1)=fileparse("$input", qr/\.[^.]*/);
   ($file_name, $i_path, $i_suffix)=fileparse("$file_name", qr/\.[^.]*/) if ($i_suffix_1 =~ /gz$/);
   my $file_ext = sprintf("%05d", $fileCount);
   my $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
   push @subfiles, $out_file_name;
   open (OUTFILE, ">" . $out_file_name ) or die ("Cannot open file for output: $!");
   my ($fh,$pid) = open_file($input);
   while (<$fh>)
   {
          die "$_\n" if (/invalid/);
          last if (eof);
          next if (/^$/);
          $name = $_;
          $name = '@'. "seq_$seq_num\n" if ($name =~ /No name/);
          $seq=<$fh>;
          $seq =~ s/\n//g;
          die "$seq\n" if ($seq=~/invalid/);
          while ($seq !~ /\+/)
          {
             die "$seq\nFormat ERROR" if ($seq=~/invalid/ or eof);
             $seq .= <$fh>;
             $seq =~ s/\n//g;
          }
          my $q_name_pos=index($seq,"+");
          $q_name = substr($seq,$q_name_pos);
          $seq = substr($seq, 0, $q_name_pos);
          my $seq_len = length $seq;
          $qual_seq=<$fh>;
          $qual_seq =~ s/\n//g;
          my $qual_seq_len = length $qual_seq;
          die "$qual_seq\n" if ($qual_seq =~ /invalid/);
          while ( $qual_seq_len < $seq_len )
          {
              die "$qual_seq\nFormat ERROR\n" if ($qual_seq =~ /invalid/ or eof);
              last if ( $qual_seq_len == $seq_len);
              $qual_seq .= <$fh>;
              $qual_seq =~ s/\n//g;
              $qual_seq_len = length $qual_seq;
          }
          $total_seq_length += $seq_len;
          print (OUTFILE "$name$seq\n$q_name\n$qual_seq\n");
          $seq_num++;
          $seqThisFile++;

          if ($seqThisFile == $subfile_size)
          {
              $fileCount++;
	      $seqThisFile = 0;
	      close (OUTFILE) or die( "Cannot close file : $!");
	      $file_ext = sprintf("%05d", $fileCount);
	      $out_file_name= $outDir . "/". $file_name . "_" . $file_ext . ".fastq";
              open (OUTFILE, ">" .  $out_file_name ) or die ("Cannot open file for output: $!") if (! eof);
              push @subfiles, $out_file_name if (! eof);
          }

   }
   my $average_len = $total_seq_length/$seq_num;
   if ( $average_len < $opt_min_L) { print "ERROR: The input ($file_name) average length $average_len < minimum cutoff length(opt_min_L) $opt_min_L.\n"; }
   close ($fh)  or die( "Cannot close file : $!");
   close (OUTFILE) or die( "Cannot close file : $!") if (! eof OUTFILE);
   return ($seq_num,$total_seq_length,@subfiles);
}
   
sub filter_adapter
{
    my ($s,$mismatchRate)=@_;
    my @match;
    $mismatchRate = $mismatchRate*100;
    my $adapter;
    my $adapter_count=0;
    my $adapter_name;
    my $s_len=length($s);
    my $pos5=0;
    my $pos3=$s_len+1;
    my %adapterSeqs=(
              "cre-loxp-forward" => "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG",
              "cre-loxp-reverse" => "AGCATATTGAAGCATATTACATACGATATGCTTCAATAATGC",
              "TruSeq-adapter-1" => "GGGGTAGTGTGGATCCTCCTCTAGGCAGTTGGGTTATTCTAGAAGCAGATGTGTTGGCTGTTTCTGAAACTCTGGAAAA",
              "TruSeq-adapter-3" => "CAACAGCCGGTCAAAACATCTGGAGGGTAAGCCATAAACACCTCAACAGAAAA",
              "PCR-primer-1"     => "CGATAACTTCGTATAATGTATGCTATACGAAGTTATTACG",
              "PCR-primer-2"     => "GCATAACTTCGTATAGCATACATTATACGAAGTTATACGA",
              "Nextera-primer-adapter-1" => "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
              "Nextera-primer-adapter-2" => "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
              "Nextera-junction-adapter-1" => "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"
              #"Nextera-junction-adapter-2" => "CTGTCTCTTATACACATCT",
              #"Nextera-junction-adapter-3" => "AGATGTGTATAAGAGACAG"
	);
    $adapterSeqs{"polyA"} = "AAAAAAAAAAAAAAAAAAAA" if ($trim_polyA);
    &read_artifactFile($artifactFile,\%adapterSeqs) if ($artifactFile);
    foreach my $key (keys %adapterSeqs)
    {
        $adapter = $adapterSeqs{$key};
        @match = String::Approx::aslice($adapter, ["i", "S ${mismatchRate}% I 0 D 0"], $s);
        if (defined $match[0][0])
        {
            $adapter_count++;
            if ($adapter_count >1){
                last;
            }
            $adapter_name=$key;
            my $index=$match[0][0];
            my $match_len=$match[0][1];
            if ($keep_5end_after_adapter_trim){
                    substr($s,$index,$s_len-$index,"");
                    $pos3=length($s)+1;
            }elsif($keep_3end_after_adapter_trim){
                    substr($s,0,$index+$match_len,"");
                    $pos5=$index+$match_len;
            }else{
                if ( int($s_len/2)-$index < ($match_len/2) )  # longer left
                {
                    if ($keep_short_after_adapter_trim){
                        substr($s,0,$index+$match_len,"");
                        $pos5=$index+$match_len;
                    }else{
                        substr($s,$index,$s_len-$index,"");
                        $pos3=length($s)+1;
                    }
                }
                else  #longer right
                {
                    if ($keep_short_after_adapter_trim){
                        substr($s,$index,$s_len-$index,"");
                        $pos3=length($s)+1;
                    }else{
                        substr($s,0,$index+$match_len,"");
                        $pos5=$index+$match_len;
                    }
                }
            }
            # same adapter sencond match
            my $match = String::Approx::amatch($adapter, ["i", "S ${mismatchRate}% I 0 D 0"], $s);
            if ($match){
                $adapter_count++;
                last;
            }
        }
    }

    if ($adapter_count>1){
        ## filter read if there are more than one adapter (or same adapter match twice) in a read
        return (1,"",$pos5,$pos3,$adapter_name);
    }elsif($adapter_count){
        return (1,$s,$pos5,$pos3,$adapter_name);
    }else{
        return (0,$s,$pos5,$pos3,"");
    }
}

sub filter_phiX
{
    my ($s,$mismatchRate)=@_;
    my $match_flag;
    $mismatchRate = $mismatchRate*100;
    $match_flag = String::Approx::amatch($s, ["i", "S ${mismatchRate}% I 0 D 0"], $phiX_seq);
    #print $mismatchRate,"\t",$match[0][0],"\n";

    if ($match_flag)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

sub read_artifactFile 
{
     my $file =shift;
     my $adapter_ref=shift;
     open (ARTIFACT,$file);
     my $id;
     my $seq;
     while(<ARTIFACT>)
     {
        chomp;
        if (/^>(\S+)/)
        {
           if ($seq)
           { 
               $adapter_ref->{$id}=$seq; 
           }
           $seq="";
           $id=$1;
        }
        else
        {
           $seq .= $_;
        }
     }
     if ($seq)
     { 
        $adapter_ref->{$id}=$seq; 
     }
     close ARTIFACT;
     return $adapter_ref;
}

sub read_phiX174
{
  my $id= "PhiX174_NC_001422";
  my $seq="GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAA
AAATTATCTTGATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGAC
TGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTT
GCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAG
TGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTAGATATGAGTCACATTTTGTT
CATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATCTGAGTCCGAT
GCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTT
TCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCG
AAGATGATTTCGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTG
CTCGTCGCTGCGTTGAGGCTTGCGTTTATGGTACGCTGGACTTTGTGGGATACCCTCGCT
TTCCTGCTCCTGTTGAGTTTATTGCTGCCGTCATTGCTTATTATGTTCATCCCGTCAACA
TTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTACGGAAAACATTATTAATGGCG
TCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTACGCGCAGGAA
ACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCGGAAGGAG
TGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTT
GCGAGGTACTAAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATT
TTAATTGCAGGGGCTTCGGCCCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGG
CGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCG
TCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGT
TGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTACTGTAGACAT
TTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAA
GGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGC
CGCTTTTCTTGGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGG
TTATTTGAATATCTATAACAACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGC
TAACCCTAATGAGCTTAATCAAGATGATGCTCGTTATGGTTTCCGTTGCTGCCATCTCAA
AAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTTTCTCGCCAAATGACGACTTC
TACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGCATACTGACCA
AGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAAC
CTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGG
CTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACA
GACCTATAAACATTCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCT
TGCGCTTGTTCGTTTTCCGCCTACTGCGACTAAAGAGATTCAGTACCTTAACGCTAAAGG
TGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTGTATGGCAACTTGCCGCCGCG
TGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGTTTAAGATTGC
TGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGA
AGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCG
CCACCATGATTATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGT
TAAATTTAATGTGACCGTTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTC
GTGATAAAAGATTGAGTGTGAGGTTATAACGCCGAAGCGGTAAAAATTTTAATTTTTGCC
GCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGCTTAGGAGTTTAATCATGTTT
CAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGTTCTCACTTCT
GTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTA
TATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATG
GATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGAT
GCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACC
CTCCCGACTGCCTATGATGTTTATCCTTTGAATGGTCGCCATGATGGTGGTTATTATACC
GTCAAGGACTGTGTGACTATTGACGTCCTTCCCCGTACGCCGGGCAATAACGTTTATGTT
GGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAAT
CAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTG
CTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAG
GCGGTCAAAAAGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATA
CTGTAGGCATGGGTGATGCTGGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACC
CTGATGAGGCCGCCCCTAGTTTTGTTTCTGGTGCTATGGCTAAAGCTGGTAAAGGACTTC
TTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGATAAGTTGCTTGATTTGGTTG
GACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTATCTTGCTGCTG
CATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGG
TTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAG
AGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTT
CACGCCAGAATACGAAAGACCAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGA
AGGAGTCTACTGCTCGCGTTGCGTCTATTATGGAAAACACCAATCTTTCCAAGCAACAGC
AGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCAAACGGCTGGTCAGTATTTTA
CCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGACTTAGTTCATC
AGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTT
CTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAG
CTGTTGCCGATACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTA
ATTTGTCTAGGAAATAACCGTCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTC
CAAATCTTGGAGGCTTTTTTATGGTTCGTTCTTATTACCCTTCTGAATGTCACGCTGATT
ATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTATTGAGGCTTGTGGCATTTCTA
CTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGCATCAAGCTCT
TGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATG
TTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGT
TAATGGATGAATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTA
TAGACCACCGCCCCGAAGGGGACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGC
AGTTTTGCCGCAAGCTGGCTGCTGAACGCCCTCTTAAGGATATTCGCGATGAGTATAATT
ACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATTGCTGGAGGCCTCCACTATGA
AATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAGGCTCATGCTG
ATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTT
ATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTG
CCGAGGGTCGCAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTG
AGTATGGTACAGCTAATGGCCGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTC
CTACAGGTAGCGTTGACCCTAATTTTGGTCGTCGGGTACGCAATCGCCGCCAGTTAAATA
GCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGG
ACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAGCCGCTTAAAG
CTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATA
TGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGC
TGTCGCTACTTCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAA
TGCTCACAATGACAAATCTGTCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACG
ACGCGACGCCGTTCAACCAGATATTGAAGCAGAACGCAAAAAGAGAGATGAGATTGAGGC
TGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACCTGTGACGACAAATCTGCTCA
AATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA";
   $seq=~ s/\n//g;
   $seq .=  ReverseComplement($seq);
   return ($id,$seq);
}

sub ReverseComplement{
        my $dna = $_[0];
    my $ReverseCompSeq = reverse ($dna);
        $ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
        return($ReverseCompSeq);
}

sub is_paired
{
    $SIG{'PIPE'} = sub{}; # avoid gunzip broken pipe
    my $paired_1=shift;
    my $paired_2=shift;
    my ($fh1,$pid1)=open_file($paired_1);
    my ($fh2,$pid2)=open_file($paired_2);
    my $count=0;
    my $check_num=1000;
    my $is_paired = 1;
    ## check top 1000 sequences paired by matching names.
    for ($count..$check_num)
    {
        my $id1=<$fh1>;
        my $seq1=<$fh1>;
        my $q_id1=<$fh1>;
        my $q_seq1=<$fh1>;
        my $id2=<$fh2>;
        my $seq2=<$fh2>;
        my $q_id2=<$fh2>;
        my $q_seq2=<$fh2>;
        my ($name1) = $id1 =~ /(\S+)/; 
        $name1 =~ s/\.\d$//;
        $name1 =~ s/\/\d$//;
        my ($name2) = $id2 =~ /(\S+)/; 
        $name2 =~ s/\.\d$//;
        $name2 =~ s/\/\d$//;
        if ($name1 ne $name2)
        {
            $is_paired=0;
        }
    }
    close $fh1;
    close $fh2;
    kill 9, $pid1;
    kill 9, $pid2;
    $SIG{'PIPE'} = 'DEFAULT';
    return $is_paired;
}

1;


