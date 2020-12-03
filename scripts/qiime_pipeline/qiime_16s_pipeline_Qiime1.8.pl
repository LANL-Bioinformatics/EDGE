#!/usr/bin/env perl
#  Based on Qiime v1.8.0
###########################################################################
# AUTHOR: CHIEN-CHI LO                                                    #                                                
# Copyright (c) 2014 LANS, LLC All rights reserved                        #
# All right reserved. This program is free software; you can redistribute #
# it and/or modify it under the same terms as Perl itself.                #
#                                                                         #
# LAST REVISED: June 2014                                                 #
###########################################################################

use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use POSIX qw(strftime);
use FindBin qw($RealBin);
use lib "$RealBin/../lib/perl5";
use lib "$RealBin/../lib/";
use File::Tee qw(tee);
#use Data::Dumper;

my $version=0.2;
my $Qiime_version="v1.8.0";
my $debug=0;

# set up environments
$ENV{PATH}= "/users/218819/scratch/tools/epd_free-7.3-2-rh5-x86_64/bin/:$ENV{PATH}";


sub Usage {
print <<"END";
Usage: perl $0 [options] -p 'reads1.fastq reads2.fastq' -b barcode.fastq -m mapping.txt -o out_directory 
     Version $version
     Input File: (can use more than once)
            -p or -u      -u     <File> Unpaired reads fastq
                          -p     <File> Paired reads in two fastq files and separate by space in quote
            
            -m            <File> MAPPING file.  #SampleID BarcodeSequence LinkerPrimerSequence sampleType...
     Output:
            -o            <PATH> Output directory.
     Options:
            -c            <INT> # of CPUs to run the script (default:4 )

            -b            <File> Barcodes fastq
            
            -barcode_len  Length of Barcode [default: 6]
 
            -q            PHRED_QUALITY_THRESHOLD
                          the maximum unacceptable Phred quality score (e.g.,
                          for Q20 and better, specify -q 19) [default: 3]
	
            -n            SEQUENCE_MAX_N
                          maximum number of N characters allowed in a sequence
                          to retain it -- this is applied after quality
                          trimming, and is total over combined paired end reads
                          if applicable [default: 1]

	    -min_per_read_length_fraction  MIN_PER_READ_LENGTH_FRACTION
                          min number of consecutive high quality base calls to
                          include a read (per single end read) as a fraction of
                          the input read length [default: 0.5]
  
            -min_otu_size the minimum otu size (in number of sequences) to
                          retain the otu [default: 2]

            -similarity   sequence similarity threshold (for blast, cdhit,
                          uclust, uclust_ref, usearch, usearch_ref, usearch61,
                          or usearch61_ref) [default: 0.94]
            
            -t            <STRING>  Project title

            -e            SAMPLING_DEPTH, --sampling_depth=SAMPLING_DEPTH
                          Sequencing depth to use for even sub-sampling and
                          maximum rarefaction depth. You should review the
                          output of the 'biom summarize-table' command to decide
                          on this value.[default: 1000]

            -debug        <bool> keep intermediate files
END

}


my $CPUs=4;
my $barcode_len=6;
my $quality_cutoff=3;
my $min_reads_q_fraction=0.5;
my $min_otu_size=2;
my $seq_max_n=1;
my $similarity=0.94;
my @paired_files;
my @unpaired_files;
my @barcode_files;
my $outDir;
my @mapping_files;
my $title;
my $sampling_depth_cutoff=1000;

my $db_path="/mnt/lustre/refdb/usrdb/Qiime";
my $reference_seqs="$db_path/gg_13_8_otus/rep_set/94_otus.fasta";
my $reference_tree="$db_path/gg_13_8_otus/trees/94_otus.tree";
my $reference_tax="$db_path/gg_13_8_otus/taxonomy/94_otu_taxonomy.txt";

# Options
GetOptions(  
           "p=s"          => \@paired_files,
           "u=s"          => \@unpaired_files,
           "b=s"          => \@barcode_files,
           "m=s"          => \@mapping_files,
           "o=s"          => \$outDir,
           "t=s"          => \$title,
           "n=i"	  => \$seq_max_n,
           "similarity=f"  => \$similarity,
           "min_otu_size=i" => \$min_otu_size,
           "min_per_read_length_fraction=f" => \$min_reads_q_fraction,
           "barcode_type=s"  => \$barcode_len,
           "q=i"          => \$quality_cutoff,
           "e=i"          => \$sampling_depth_cutoff,
           "c|cpu=i"      => \$CPUs,
           "debug"        => \$debug,
           "version"      => sub{print "Version: $version\n";exit;},
           "help|?"       => sub{Usage()} );


@paired_files = split(/,/,join(',',@paired_files));
@unpaired_files = split(/,/,join(',',@unpaired_files));
@mapping_files = split(/,/,join(',',@mapping_files));
@barcode_files = split(/,/,join(',',@barcode_files));

die "Missing input files.\n", &Usage() unless @unpaired_files or @paired_files or @barcode_files;
die "Missing output directory.\n",&Usage() unless $outDir;
die "Missing mapping file.\n",&Usage() unless @mapping_files;

`mkdir -p $outDir`;
$outDir=Cwd::abs_path("$outDir");

my $error_log_file="$outDir/errorLog.txt";
my $process_log_file = "$outDir/processLog.txt";
my ($error_log_fh, $process_log_fh)=&_logProcess($error_log_file,$process_log_file);

my $parameter_file = "$outDir/parameter.txt";

`echo \"pick_otus:similarity $similarity\" >> $parameter_file`;

$parameter_file = " -p $parameter_file ";
#&check_mapping_files(\@mapping_files);

my $before_demultiplex_fastq;
my $before_demultiplex_barcode;
my $before_demultiplex_mapping;
my $merge_mapping_file= &merge_mapping_file;

my @new_fastq;
if (!@barcode_files)
{
        my @old_fastq = (@paired_files)?@paired_files:@unpaired_files;
        my ($fastq_files_r,$barcode_files_r)=&get_barcode_internal_fastq(\@old_fastq,$barcode_len);
        if (@paired_files) 
        {
            @paired_files=@{$fastq_files_r};
            $barcode_len = $barcode_len * 2;
        }
        else
        {
            @unpaired_files=@{$fastq_files_r};
        }
        @barcode_files=@{$barcode_files_r};
}

if (@paired_files)
{
    ($before_demultiplex_fastq,$before_demultiplex_barcode)=&join_paired_ends;
    $before_demultiplex_mapping .= "$merge_mapping_file," foreach @paired_files;
    $before_demultiplex_mapping =~ s/,$//;
}
else
{
    $before_demultiplex_fastq = join(',',@unpaired_files);
    $before_demultiplex_barcode = join(',',@barcode_files);
    $before_demultiplex_mapping = join(',',@mapping_files);    
}

&split_libraries($before_demultiplex_fastq,$before_demultiplex_barcode,$before_demultiplex_mapping);

my ($biom,$rep_set_fna,$rep_set_tre)=&pick_otus;

my ($sampling_depth,$num_sample,$otu_tables)=&determine_sampling_depth_and_otu_table($biom);

$sampling_depth = $sampling_depth_cutoff if ($sampling_depth_cutoff);

my $analysis_outputDir="$outDir/analysis";

system("rm -rf $analysis_outputDir") if (-e $analysis_outputDir);
`mkdir -p $analysis_outputDir`;
#system("gzip -fc $rep_set_fna > $analysis_outputDir/rep_set.fna.gz");
#system("gzip -fc $rep_set_tre > $analysis_outputDir/rep_set.tre.gz");
system("cp $rep_set_fna $analysis_outputDir/rep_set.fna");
system("cp $rep_set_tre $analysis_outputDir/rep_set.tre");
system("gzip -fc $biom > $analysis_outputDir/table.biom.gz");
#system("gzip -fc $otu_tables > $analysis_outputDir/OTUs.table.txt.gz");
system("cp $otu_tables $analysis_outputDir/OTUs.table.txt");
if ($num_sample>1)
{
    &core_diversity_analysis($biom,$sampling_depth,$analysis_outputDir);
    &heatmap($biom,$analysis_outputDir);
    &update_index_html;
}
else
{
    &heatmap($biom,$analysis_outputDir);
    &summarize_taxa($biom,$analysis_outputDir);
 #   &alpha_diversity_analysis;
}


close $error_log_fh;
close $process_log_fh;

system("rm -rf jobs") if (-e "jobs");

print "\nTotal" ;
&printRunTime($^T);


### END MAIN ###
### Below are subroutins ###
sub summarize_taxa
{
    my $biom = shift;
    my $dir = shift;
    my $outputDir="$dir/taxa_summary";
    my $cmd="summarize_taxa.py -i $biom -o $outputDir"; 
    if (-e "$outputDir/taxa.finish")
    {
        print "\nTaxanomy summarized\n";
        return 0;
    }
    &process_cmd($cmd,"Generating taxanomy summary");
    system("touch $outputDir/taxa.finish");
}

sub heatmap
{
    my $biom = shift;
    my $dir = shift;
    my $outputDir="$dir/heatmap";
    my $cmd="make_otu_heatmap_html.py -i $biom -o $outputDir"; 
    if (-e "$outputDir/heatmap.finish")
    {
        print "\nHeatMap generated\n";
        return 0;
    }
    &process_cmd($cmd,"Generating Heatmap");
    system("touch $outputDir/heatmap.finish");
}

sub get_barcode_internal_fastq
{
    my $files_ref=shift;
    my $barcode_len=shift;
    my $outputDir="$outDir/newFastq";
    `mkdir -p $outputDir`;
    my @return_fastq_files;
    my @return_barcodes_files;
    my $fastq_file_list = "$outputDir/fastq.list";
    my $barcode_file_list = "$outputDir/barcode.list";
    if (-e "$outputDir/getBarcode.finish")
    {
        @return_fastq_files=&readFileList($fastq_file_list);
        @return_barcodes_files=&readFileList($barcode_file_list);
        return (\@return_fastq_files,\@return_barcodes_files);
    }
    print "\nGet barcode internal fastq\n";
    foreach my $file (@{$files_ref})
    {
        if (@paired_files)
        {
            my ($R1,$R2)=split /\s+/,$file;
            my ($R1_no_barcodes_fastq,$R1_barcodes)=&split_barcode($R1,$barcode_len,$outputDir);
            my ($R2_no_barcodes_fastq,$R2_barcodes)=&split_barcode($R2,$barcode_len,$outputDir);
            my $barcodes=&catBarcodes($R1_barcodes,$R2_barcodes,$outputDir);
            push @return_fastq_files, "$R2_no_barcodes_fastq $R1_no_barcodes_fastq";
           # push @return_fastq_files, "$R1_no_barcodes_fastq $R2_no_barcodes_fastq";
            push @return_barcodes_files, $barcodes;
        }
        else
        {
            my ($file_no_barcodes_fastq,$barcodes)=&split_barcode($file,$barcode_len,$outputDir);
            push @return_fastq_files, $file_no_barcodes_fastq;
            push @return_barcodes_files, $barcodes;
        }
    }
    open (my $fh, ">$fastq_file_list") or die "Cannot write $fastq_file_list\n";
    open (my $fh2, ">$barcode_file_list") or die "Cannot write $barcode_file_list\n";
    for (0..$#return_barcodes_files)
    {
        print $fh $return_barcodes_files[$_],"\n";
        print $fh2 $return_barcodes_files[$_],"\n";
    }
    close $fh;
    close $fh2;
    system("touch $outputDir/getBarcode.finish");
    return (\@return_fastq_files,\@return_barcodes_files);
}

sub open_file
{
    my ($file) = @_;
    my $fh;
    if ( $file=~/\.gz$/i ) { open($fh, "gunzip -c $file |") or die ("gunzip -c $file: $!"); }
    else { open($fh,'<',$file) or die("$file: $!"); }
    return $fh;
}


sub readFileList
{
    my $file = shift;
    my @files;
    open(my $fh,$file) or die "Cannot open $file\n";
    while(<$fh>)
    {
        chomp;
        $_ =~ s/ //g;
        push @files, $_;
    }
    close $fh;
    return(@files);
}

sub catBarcodes 
{
    my $b1=shift;
    my $b2=shift;
    my $outputDir=shift;
    my ($name,$path,$suffix) = fileparse($b1,qr/\.barcodes\.[^.]*/);
    my $barcode_file = "$outputDir/$name.catbarcodes.fastq";
    open(my $ofh,">$barcode_file") or die "Cannot write $barcode_file\n";
    my $ib1= open_file($b1);
    my $ib2= open_file($b2);
    while(<$ib1>)
    {
        my $id=$_;
        my $seq=<$ib1>;
        chomp $seq;
        my $q_id=<$ib1>;
        my $q_seq=<$ib1>;
        chomp $q_seq;
        my $b2_id=<$ib2>;
        my $b2_seq=<$ib2>;
        chomp $b2_seq;
        my $b2_q_id=<$ib2>;
        my $b2_q_seq=<$ib2>;
        chomp $b2_q_seq;
        #print $ofh "$id$seq$b2_seq\n$q_id$q_seq$b2_q_seq\n";
        print $ofh "$b2_id$seq$b2_seq\n$b2_q_id$q_seq$b2_q_seq\n";
    }
    close $ib1;
    close $ib2;
    close $ofh;
    return $barcode_file;
}
sub split_barcode
{
    my $file=shift;
    my $barcode_len=shift;
    my $outputDir=shift;
    my $fh = open_file($file);
    my ($name,$path,$suffix) = fileparse($file,qr/\.[^.]*/);
    my $barcode_file="$outputDir/$name.barcodes.fastq";
    my $new_fastq= "$outputDir/$name.nobarcodes.fastq";
    open (my $oBarcode,">$barcode_file") or die "Cannot write $barcode_file\n";
    open (my $oFastq,">$new_fastq") or die "Cannot write $new_fastq\n";
    while(<$fh>)
    {
        my $id=$_;
        my $seq=<$fh>;
        chomp $seq;
        my $q_id=<$fh>;
        my $q_seq=<$fh>;
        chomp $q_seq;
        my $barcode_seq = substr($seq,0,$barcode_len,"");   
        my $barcode_q_seq = substr($q_seq,0,$barcode_len,"");
        print $oBarcode "$id$barcode_seq\n$q_id$barcode_q_seq\n";
        print $oFastq "$id$seq\n$q_id$q_seq\n";
    } 
    close $fh;
    close $oBarcode;
    close $oFastq;
    return ($new_fastq,$barcode_file);
}

sub check_mapping_files
{
    my $mapping_files = shift;
    my $outputDir="$outDir/checkMappingFile";
 #   if (-e "$outputDir/check.finished")
  #  {
   #     print "\nChecking Mapping Files Done\n";
    #    return 0;
    #}
    system("rm -rf $outputDir");
    `mkdir -p $outputDir`;
    foreach my $file ( @{$mapping_files})
    {
        my ($name,$path,$suffix) = fileparse($file,qr/\.[^.]*/);
        my $log = "$outputDir/$name.log";
        my $cmd="validate_mapping_file.py -m $file -o $outputDir";
        &process_cmd($cmd,"Checking $file");
        my $ret=`grep Error $log`;
        if ($ret) 
        {
            exit;
        }
    }
   # system("touch $outputDir/check.finished");
}

sub join_paired_ends 
{
    my $outputDir="$outDir/fastq-join_joined";
    my $tmp_dir="$outputDir/pair";
    my $joined_barcode;
    my $joined_file;
    my @join_files;
    my @join_barcode_files;
    `mkdir -p $outputDir`;
    if (-e "$outputDir/join.finished")
    {
        print "\nJoining Paired-end Files Done\n";
        open (my $fh, "$outputDir/join.finished") or die "Cannot open $outputDir/join.finished\n";
        $joined_file = <$fh>;   chomp $joined_file;
        $joined_barcode = <$fh>;  chomp $joined_barcode;
        return($joined_file,$joined_barcode);
    }
    foreach my $file_i (0..$#paired_files)
    {
        my ($R1,$R2)=split /\s+/,$paired_files[$file_i];   
        $R1=Cwd::abs_path("$R1");
        $R2=Cwd::abs_path("$R2"); 
        my $cmd="join_paired_ends.py -f $R1 -r $R2 -b $barcode_files[$file_i] -o $tmp_dir$file_i";
        &process_cmd($cmd,"Joining $paired_files[$file_i]");
        system ("echo \"$R1 $R2\" > $tmp_dir$file_i/file.txt");
        push @join_files, "$tmp_dir$file_i/fastqjoin.join.fastq"; 
        push @join_barcode_files, "$tmp_dir$file_i/fastqjoin.join_barcodes.fastq";
    }
    #if (scalar(@paired_files)>1)
    #{
    #    &process_cmd("cat ${tmp_dir}*/*join.fastq > $joined_file", "concatenate joined fastq files");
    #    &process_cmd("cat ${tmp_dir}*/*join_barcodes.fastq > $joined_barcode", "concatenate barcode update fastq files");
    #}
    #else
    #{
    #    system("mv $outDir/tmp0/fastqjoin.join.fastq  $joined_file");
    #    system("mv $outDir/tmp0/fastqjoin.join_barcodes.fastq $joined_barcode");
    #}
    #system("rm -rf ${tmp_dir}*");
    $joined_file = join(",",@join_files);
    $joined_barcode = join(",",@join_barcode_files);
    system("echo $joined_file > $outputDir/join.finished");
    system("echo $joined_barcode >> $outputDir/join.finished");
    
    return ($joined_file,$joined_barcode);
}

sub merge_mapping_file
{
    my $merged_mapping_file="$outDir/combined_mapping.txt";
    my @tmp = $merged_mapping_file;
    if (-s $merged_mapping_file)
    {
        print "\nMerging Mapping Files Done\n";
        &check_mapping_files(\@tmp);
        return $merged_mapping_file;
    }
    my $cmd="merge_mapping_files.py -o $merged_mapping_file -m ". join( ',' ,@mapping_files);
    if (scalar(@mapping_files)>1)
    {
        &process_cmd($cmd,"Merge mapping files");
    }
    else
    {
        system("cp $mapping_files[0] $merged_mapping_file");
    }
    &check_mapping_files(\@tmp);
    return $merged_mapping_file;
}

sub split_libraries 
{
     my $outputDir="$outDir/slout";
     my $tmpDir = "$outputDir/sloutTmp";
     my $fastq_file=shift;
     my $barcode_file=shift;
     my $mapping_file=shift;
     `mkdir -p $outputDir`;
     if (-e "$outputDir/split.finished")
     {
        print "\nDe-multiplexing Done\n";
        return 0;
     }
     my $cmd;
     $cmd="split_libraries_fastq.py -o $outputDir -i $fastq_file -b $barcode_file -m $mapping_file --barcode_type $barcode_len -q $quality_cutoff -p $min_reads_q_fraction -n $seq_max_n ";
     &process_cmd($cmd,"De-multiplexing the samples");
     
     system("touch $outputDir/split.finished");
}

sub pick_otus
{
    my $outputDir="$outDir/otus";
    my $biom="$outputDir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom";
    my $rep_set_fna="$outputDir/rep_set.fna";
    my $rep_set_tre="$outputDir/rep_set.tre";
    my $seq = "$outDir/slout/seqs.fna";
    if (-e "$outputDir/pickOTU.finished")
    {
        print "\nPicking OTUs Done\n";
        return ($biom,$rep_set_fna,$rep_set_tre);
    }
    my $cmd="pick_open_reference_otus.py -f -i $seq -r $reference_seqs -o $outputDir -s 0.01 $parameter_file --min_otu_size $min_otu_size -aO $CPUs";
    &process_cmd($cmd,"Pick OTUs");
    system("touch $outputDir/pickOTU.finished");
    return($biom,$rep_set_fna,$rep_set_tre);
}

sub determine_sampling_depth_and_otu_table
{
    my $biom=shift;;
    my $biom_summary="$outDir/otus/biom_table_summary.txt";
    my $otu_tables="$outDir/otus/otu_table_tabseparated.txt";
    my $sampling_depth;
    my $num_sample;
    my $cmd="biom summarize-table -i $biom -o $biom_summary --suppress-md5";
    my $otu_table_cmd="biom convert -i $biom -o $otu_tables -b --header-key=\"taxonomy\"";
    if ( ! -e $otu_tables)
    {
        &process_cmd($otu_table_cmd,"Generating otus table");
    }
    if ( ! -e $biom_summary)
    {
        &process_cmd($cmd,"Generating biom table summmary");
    }
    ($sampling_depth,$num_sample)=&get_depth_cutoff_from_biom_summary_table($biom_summary);
    return ($sampling_depth, $num_sample,$otu_tables);
}

sub get_depth_cutoff_from_biom_summary_table {

    my $biom_summary = shift;
    my $sampling_depth;
    my $num_sample;
    open (my $fh, $biom_summary) or die "Cannot open $biom_summary\n";
    while(<$fh>)
    {
        ($num_sample) = $_ =~ /Num samples:\s+(\d+)/ if (!$num_sample);
        last if (/detail/);
    }
    my @sample_counts;
    my $pass_cutoff=0;
    while(<$fh>)
    {
        chomp;
        $_ =~ s/^\s+//;
        my ($id,$value) = split /\s+/,$_;
        next if( $value<100); # too few to be analyzed;
	$pass_cutoff++ if( $value > $sampling_depth_cutoff); 
        push @sample_counts, int($value);
    }
    close $fh;
    if ($pass_cutoff == 0 ){
        die "No sample size  > $sampling_depth_cutoff. Stop process. See file $biom_summary";
    }
    my ($mean_depth,$std,$median_depth,$q1,$q3) = &statistical_calculation_on_array(@sample_counts);
    $sampling_depth = $mean_depth - 2*$std;
    if ($mean_depth < 100)
    {
	die "Average sample size < 100. Too small to proceed. See file $biom_summary\n";
    }
    if ($sampling_depth<100)
    {
        my $outlier_cutoff = $q1 - 1.5*($q3-$q1);
        $sampling_depth = ($outlier_cutoff>100)? $outlier_cutoff:$q1-100;
    }
    return ($sampling_depth, $num_sample);
}

sub statistical_calculation_on_array {
    my (@values) = @_;
    my $total_elements = scalar(@values);
    #Prevent division by 0 error in case you get junk data
    return undef unless($total_elements);

    my @numbers = sort {$a<=>$b} @values;
    # Step 1, find the mean of the numbers
    my $total1 = 0;
    foreach my $num (@numbers) {
        $total1 += $num;
    }
    my $mean1 = $total1 / $total_elements;

    # Step 2, find the mean of the squares of the differences
    # between each number and the mean
    my $total2 = 0;
    foreach my $num (@numbers) {
        $total2 += ($mean1-$num)**2;
    }
    my $mean2 = $total2 / $total_elements;

    # Step 3, standard deviation is the square root of the
    # above mean
    my $std_dev = sqrt($mean2);

    my $median;
    if ($total_elements % 2) {
        $median = $numbers[int($total_elements/2)];
    } else {
        $median = ($numbers[$total_elements/2] + $numbers[$total_elements/2 - 1]) / 2;
    }
    my $q1=$numbers[int($total_elements/4)];  # 1st quartile
    my $q3=$numbers[int((3*$total_elements)/4)]; # 3rd quartile

    return ($mean1,$std_dev,$median);
}

sub core_diversity_analysis 
{
    my $biom=shift;;
    my $sampling_depth=shift;
    my $outputDir=shift;
    my $DieCatch=1;
    if (-e "$outputDir/diversity.finished")
    {
        print "\nDiversity Analysis Done\n";
        return 0;
    }
    my $cmd = "core_diversity_analyses.py --recover_from_failure -i $biom -o $outputDir -m $merge_mapping_file -e $sampling_depth -t $reference_tree $parameter_file -aO $CPUs";
    my $return=&process_cmd($cmd,"Perform Diversity analysis and taxanomy summary",$DieCatch);
    if ($return)
    {
	$cmd .= " -w > core_analysis_commands.sh";
	&process_cmd($cmd,"Perform Diversity analysis and taxanomy summary");
	&process_cmd("sh core_analysis_commands.sh", "skip error step and run remaining analysis");
	unlink "core_analysis_commands.sh";
    }
    system("touch $outputDir/diversity.finished");
    system("cp $outputDir/index.html $outputDir/index.html.org");
}

sub update_index_html
{
    my $index = "$outDir/analysis/index.html";
    system ("rm -rf $outDir/analysis/checkMappingFile");
    system ("mv $outDir/checkMappingFile $outDir/analysis/checkMappingFile");
    open (my $fh, "$index.org") or die "Cannot open $index.org\n";
    open (my $ofh, ">$index") or die "Cannot write $index\n";
    while(<$fh>)
    {
        if (/^<table/)
        {
            print $ofh "<h1>Project: $title</h1>\n";
            print $ofh $_;
        }
        elsif(/Master/)
        {
            print $ofh $_;
            print $ofh "<tr><td>Mapping File Check</td><td> <a href=\"./checkMappingFile/combined_mapping.html\" target=\"_blank\">combined_mapping.html</a></td></tr>\n";
            print $ofh "<tr><td>Representative Sequences</td><td> <a href=\"./rep_set.fna\" target=\"_blank\">rep_set.fna</a></td></tr>\n";
            print $ofh "<tr><td>Representative Sequences Tree</td><td> <a href=\"./rep_set.tre\" target=\"_blank\">rep_set.tre</a></td></tr>\n";
            print $ofh "<tr><td>OTUs Table</td><td> <a href=\"./OTUs.table.txt\" target=\"_blank\">OTUs.table.txt</a></td></tr>\n";
            print $ofh "<tr><td>Biom Table</td><td> <a href=\"./table.biom.gz\" target=\"_blank\">table.biom.gz</a></td></tr>\n";
            print $ofh "<tr><td>OTUs HeatMap</td><td> <a href=\"./heatmap/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.html\" target=\"_blank\">heatmap.html</a></td></tr>\n";
        }
        else
        {
            print $ofh $_;
        }
    }
    close $fh;
    close $ofh;
}


sub process_cmd {
    my ($cmd, $msg, $dieCatch) = @_;


    if ($msg) {
        print "\n\n";
        print "###########################################################################\n";
        print &getTimeString."  $msg\n";
        print "###########################################################################\n";
    }
    
    print "CMD: $cmd\n";
    if ($debug) {
        print STDERR "\n\n-WAITING, PRESS RETURN TO CONTINUE ...";
        my $wait = <STDIN>;
        print STDERR "executing cmd.\n\n";
        
    }
    

    my $time_start = time();
    
    my $ret = system($cmd);
    my $time_end = time();

    if ($ret) {
	if ($dieCatch)
        {
	    print $ret;
        }else
	{
	    die "Error, CMD: $cmd died with ret $ret";
        }
    }

    my $number_minutes = sprintf("%.1f", ($time_end - $time_start) / 60);
    
 #   print "TIME: $number_minutes min. for $cmd\n";
    print "TIME: $number_minutes min.\n";
    

    return $ret;
}

sub printRunTime {
  my $time=shift;
  my $runTime = time() - $time;
  printf(" Running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime % 3600) / 60), 
  int($runTime % 60));
}

sub getTimeString
{
    my $now_string = strftime "[%Y %b %e %H:%M:%S]", localtime;
    #$now_string =~ s/\s//g;
    return $now_string;
}

sub _logProcess
{
    my $error_log_file = shift;
    my $process_log_file =shift;
    # capture error log
    open (my $error_log_fh, ">> $error_log_file") or die "Failed to open $error_log_file\n$!";
    $SIG{__WARN__} = sub {print $error_log_fh @_;print STDOUT @_};
    $SIG{__DIE__} = sub {print $error_log_fh @_;print STDOUT @_;exit 1};
    #Caputre log
    open (my $process_log_fh, ">> $process_log_file") or die "Failed to open $process_log_file\n$!";
    tee(STDOUT, '>>', "$process_log_file");
    tee(STDERR, '>>', "$error_log_file");
    # print whole running command and config file to process log
    print $process_log_fh "\n". &getTimeString." Job Start\n";
    print $process_log_fh qx/ps -o args $$/;
    print $process_log_fh "\nPipeline Version: $version based on Qiime $Qiime_version\n";
    print $process_log_fh "\nProject Dir: $outDir\n";
    return ($error_log_fh, $process_log_fh);
}
