#!/usr/bin/env perl
#  Based on Qiime v1.9.1
###########################################################################
# AUTHOR: CHIEN-CHI LO                                                    #                                                
# Copyright (c) 2014-2015 LANS, LLC All rights reserved                   #
# All right reserved. This program is free software; you can redistribute #
# it and/or modify it under the same terms as Perl itself.                #
#                                                                         #
# LAST REVISED: Sep 2015                                                  #
###########################################################################

use strict;
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use File::Basename;
use POSIX qw(strftime);
use FindBin qw($RealBin);
use lib "$RealBin/../lib/perl5";
use lib "$RealBin/../lib/";
use File::Copy;
#use File::Tee qw(tee);
use Data::Dumper;

my $version=0.6;
my $Qiime_version="v1.9.1";
my $debug=0;

my $working_dir = getcwd;
my $EDGE_HOME =($ENV{EDGE_HOME})? $ENV{EDGE_HOME}: abs_path("$RealBin/../../");
# set up environments
$ENV{PATH}= "$EDGE_HOME/thirdParty/Anaconda2/bin:/scratch-218819/apps/anaconda2/bin:$RealBin/scripts/:$ENV{PATH}";

sub Usage {
	my $msg=shift;
	my $short_usage = "perl $0 [options] -p reads1.fastq reads2.fastq -m mapping.txt -o out_directory ";
($msg)?
print "\n $msg\n\n $short_usage\n\n Option -h to see full usage \n\n":
print <<"END";
Usage: $short_usage
     Version $version
     Based on Qiime $Qiime_version
     Input File: (can use more than once)
            -p or -u      -u     <File> Unpaired reads fastq
                          -p     <File> Paired reads in two fastq files and separate by space
               or -d      -d     <PATH/directory> Contains Demultiplexed fastq files. To use this option, 
                                 the mapping file need a 'Files' column with filenames for each sampleID.
            
            -m            <File> MAPPING file.  #SampleID BarcodeSequence LinkerPrimerSequence sampleType ... Description
			  comma separated for multiple mapping file
     Output:
            -o            <PATH> Output directory.
     Options:
            -target       Greengenes, SILVA, or ITS. [default: Greengenes]
                          Greengenes is for 16s.
                          SILVA is for 16s/18s.
                          ITS is from https://unite.ut.ee/ and for fungal rDNA ITS sequences.

            -b            <File> Barcodes fastq
            
            -barcode_len  Length of Barcode [default: 6]

            -pe_orientation  Paired Reads Orientation; FR (default) or RF
                          FR: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) 
                              is in the antisense strand (reverse)
                          RF: first read (/1) of fragment pair is sequenced as anti-sense (reverse(R)), and 
                              second read (/2) is in the sense strand (forward(F));

            -demultiplex_fa  DE-MULTIPLEX FASTA file 
                          provide multiple previous demultiplex fasta by comma
			  separated files and will bypass the split fastq step. 
			  The fasta file will be concatenated as input for OTU clustering.
			  e.g: /path/run1/seqs.fna,/path/run2/seqs.fna

	    -UPARSE	  <bool> use UPARSE pipeline clusters NGS amplicon reads into OTUs
			  Edgar, R.C. (2013) Nature Methods Pubmed:23955772
 
            -q            PHRED_QUALITY_THRESHOLD
                          the maximum unacceptable Phred quality score (e.g.,
                          for Q20 and better, specify -q 19) [default: 3]

            -phred_offset The ascii offset to use when decoding phred scores (33 or 64)
                          [default: determined automatically]
	
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

 	    -filter_taxa           
			  Comma-separated list of taxa to discard from OTU table
			  e.g: p__Bacteroidetes,p__Firmicutes

            -substract_NTC   A LIST OF SAMPLE ID separated by comma
		          substarct observation count from No Template Control (NTC)
			  and remove NTC samples from OTU table

            -e            Sequencing_depth_min_cutoff
                          Filter sample less this amount of sequences.
                          The minimium of sequenceing depth of samples after
                          this filter will be  use for even sub-sampling and
                          maximum rarefaction depth. You should review the
                          output of the 'biom summarize-table' command to decide
                          on this value.[default: 1000]

            -chart_type   This is the type of chart to plot (i.e. pie, bar or
                          area). The user has the ability to plot multiple
                          types, by using a comma-separated list (e.g. area,pie)
                          [default: area,bar]

            -c            <INT> # of CPUs to run the script (default:4 )

            -t            <STRING>  Project title

            -debug        <bool> keep intermediate files
END
exit(1);
}


my $CPUs=4;
my $barcode_len=6;
my $pe_orientation="fr";
my $quality_cutoff=3;
my $phred_offset;
my $min_reads_q_fraction=0.5;
my $min_otu_size=2;
my $seq_max_n=1;
my $similarity=0.94;
my @paired_files;
my @unpaired_files;
my @barcode_files;
my $demultiplex_dir;
my @demultiplex_fa_files;
my $outDir;
my @mapping_files;
my $title;
my $sampling_depth_cutoff=1000;
my $negative_taxa;
my $ntc_list;
my $UPARSE_opt;
my $target_ref;
my $chartType="area,bar";
#my $db_path="/mnt/lustre/refdb/usrdb/Qiime";
#my $reference_seqs="$db_path/gg_13_8_otus/rep_set/94_otus.fasta";
#my $reference_tree="$db_path/gg_13_8_otus/trees/94_otus.tree";
#my $reference_tax="$db_path/gg_13_8_otus/taxonomy/94_otu_taxonomy.txt";

# Options
GetOptions(  
           "p=s{,}"          => \@paired_files,
           "u=s"          => \@unpaired_files,
           "b=s"          => \@barcode_files,
           "m=s"          => \@mapping_files,
           "d|demultiplex_dir=s"   => \$demultiplex_dir,
           "target=s"     => \$target_ref,
	   "demultiplex_fa=s"  =>\@demultiplex_fa_files,
	   "filter_taxa=s" => \$negative_taxa,
	   "substract_NTC=s"	=> \$ntc_list,
           "pe_orientation=s"	=> \$pe_orientation,
	   "UPARSE"	  => \$UPARSE_opt,
           "o=s"          => \$outDir,
           "t|title=s"          => \$title,
           "n=i"	  => \$seq_max_n,
           "similarity=f"  => \$similarity,
           "min_otu_size=i" => \$min_otu_size,
           "min_per_read_length_fraction=f" => \$min_reads_q_fraction,
           "barcode_type=s"  => \$barcode_len,
           "q=i"          => \$quality_cutoff,
           "phred_offset=i" => \$phred_offset,
           "chart_type=s" => \$chartType,
           "e=i"          => \$sampling_depth_cutoff,
           "c|cpu=i"      => \$CPUs,
           "debug"        => \$debug,
           "version"      => sub{print "Version: $version\n";exit;},
           "help|?"       => sub{Usage()} );


@paired_files = split(/,/,join(',',@paired_files));
@unpaired_files = split(/,/,join(',',@unpaired_files));
@mapping_files = split(/,/,join(',',@mapping_files));
@barcode_files = split(/,/,join(',',@barcode_files));
@demultiplex_fa_files = split(/,/,join(',',@demultiplex_fa_files));
$pe_orientation = lc $pe_orientation;

unless (@demultiplex_fa_files or $demultiplex_dir){
	&Usage("Missing input files.") unless @unpaired_files or @paired_files or @barcode_files;
}
&Usage("Missing output directory at flag -o .") unless $outDir;
&Usage("Missing mapping file.") unless @mapping_files;

my $id_file;
if ($demultiplex_dir){
    $demultiplex_dir=Cwd::abs_path("$demultiplex_dir");
    my ($pe_files, $se_files)=&parse_mapping_files($demultiplex_dir,\@mapping_files);
    @paired_files = @{$pe_files};
    @unpaired_files = @{$se_files};
}


my @make_paired_paired_files;
my %file;
if (@paired_files)
{
    if (scalar(@paired_files) % 2) { Usage("Please check paired data input are even file numbers\n") ;}
    map { if(is_file_empty($_)){ Usage("Please check paired data input at flag -p.\n $_ doesn't not exist or empty.");}
	  if( ! is_fastq($_)){Usage("$_ is not in fastq format"); } 
	  $file{basename($_)}=1;} @paired_files;
    #make pair in a new array 'read1_1 read1_2', 'read2_1 read2_2' ...
    for(my$i=0;$i<=$#paired_files;$i=$i+2)
    {            
        if (&is_paired($paired_files[$i], $paired_files[$i+1]))
        {
            push @make_paired_paired_files, "$paired_files[$i] $paired_files[$i+1]";
        }
        else
        {
            print ("The sequence names of the paired end reads in $paired_files[$i],$paired_files[$i+1] are not matching.\nWill use them as si
ngle end reads\n");
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
	  if ( ! is_fastq($_))
          {
              Usage("$_ is not in fastq format");
          }
        } @unpaired_files;
}

`mkdir -p $outDir`;
$outDir=Cwd::abs_path("$outDir");

######  LOG  ######
    my $error_log_file="$outDir/errorLog.txt";
    my $process_log_file = "$outDir/processLog.txt";
    # capture error log
    open ( LOG, ">>", $process_log_file) or die "Failed to write $process_log_file\n$!";
    open(STDERR, '>&', STDOUT) or die "Can't redirect stderr: $!";
    open(STDERR, '>', $error_log_file) or die "Can't redirect stderr: $!";
    $SIG{__WARN__} = sub {print STDERR @_; &lprint  (@_)};
    $SIG{__DIE__} = sub {print STDERR @_; &lprint  (@_);exit 1};
    # print whole running command and config file to process log
    &lprint("\nProject Start: ".  &getTimeString."\n");
    print LOG qx/ps -o args $$/;
    &lprint("\nPipeline Version: $version based on Qiime $Qiime_version\n");
    &lprint("\nProject Dir: $outDir\n");

######

my $parameter_file = "$outDir/parameter.txt";
my $pick_otu_ref;
my $core_diveristy_options="--recover_from_failure";

`echo \"pick_otus:similarity $similarity\" > $parameter_file`;
`echo \"pick_otus:enable_rev_strand_match True\" >> $parameter_file`;
`echo plot_taxa_summary:chart_type $chartType >>$parameter_file`;

my $silva_version=119;
if ($target_ref =~ /SILVA/i){
	if ($silva_version == 104){
		`echo assign_taxonomy:id_to_taxonomy_fp $RealBin/data/silva_104/Silva_taxa_mapping_104set_97_otus.txt >> $parameter_file`;
		`echo assign_taxonomy:reference_seqs_fp $RealBin/data/silva_104/silva_104_rep_set.fasta >> $parameter_file`;
		`echo align_seqs:template_fp $RealBin/data/silva_104/core_Silva_aligned.fasta >>$parameter_file`;
		$pick_otu_ref=" -r $RealBin/data/silva_104/silva_104_rep_set.fasta ";
	}else{
		`echo assign_taxonomy:id_to_taxonomy_fp $RealBin/data/Silva119_release/taxonomy/97/taxonomy_97_raw_taxa.txt >> $parameter_file`;
		`echo assign_taxonomy:reference_seqs_fp $RealBin/data/Silva119_release/rep_set/97/Silva_119_rep_set97.fna >> $parameter_file`;
		`echo align_seqs:template_fp $RealBin/data/Silva119_release/core_alignment/core_Silva119_alignment.fna >>$parameter_file`;
		$pick_otu_ref=" -r $RealBin/data/Silva119_release/rep_set/97/Silva_119_rep_set97.fna ";
	}
	`echo filter_alignment:suppress_lane_mask_filter True >>$parameter_file`;
	`echo filter_alignment:entropy_threshold 0.10 >>$parameter_file`;
	`echo filter_alignment:allowed_gap_frac 0.80 >>$parameter_file`;
}
if ($target_ref =~ /ITS/i){
	`echo assign_taxonomy:id_to_taxonomy_fp $RealBin/data/its_12_11_otus/taxonomy/97_otu_taxonomy.txt >> $parameter_file`;
	`echo assign_taxonomy:reference_seqs_fp $RealBin/data/its_12_11_otus/rep_set/97_otus.fasta >> $parameter_file`;
	$pick_otu_ref=" -r $RealBin/data/its_12_11_otus/rep_set/97_otus.fasta --suppress_align_and_tree ";
	$core_diveristy_options .= " --nonphylogenetic_diversity ";
}

$parameter_file = " -p $parameter_file ";
#&check_mapping_files(\@mapping_files);

my $before_demultiplex_fastq;
my $before_demultiplex_barcode;
my $before_demultiplex_mapping;
my $merge_mapping_file= &merge_mapping_file;
my @new_fastq;
if(@demultiplex_fa_files){
	my $outputDir="$outDir/slout";
	`mkdir -p $outputDir`;
	&lprint("Concatenate demultiplexed fasta files\n");
	map {`cat $_ >> $outputDir/seqs.fna`; } @demultiplex_fa_files;
}else{
## if provide demultiplex fasta. by pass the split fastq step and barcode handle
	if (!@barcode_files and !$demultiplex_dir)
	{
	        my @old_fastq = (@make_paired_paired_files)?@make_paired_paired_files:@unpaired_files;
	        my ($fastq_files_r,$barcode_files_r)=&get_barcode_internal_fastq(\@old_fastq,$barcode_len);
	        if (@make_paired_paired_files) 
	        {
	            @make_paired_paired_files=@{$fastq_files_r};
	            $barcode_len = $barcode_len * 2;
	        }
	        else
	        {
	            @unpaired_files=@{$fastq_files_r};
	        }
	        @barcode_files=@{$barcode_files_r};
	}
	
	if (@make_paired_paired_files)
	{
		# QC before join
#	FaQCs.pl -p CBM_31-51_ilmnA_Undetermined_000000000-ADAN0_L001_R1_001.fastq.gz CBM_31-51_ilmnA_Undetermined_000000000-ADAN0_L001_R2_001.fastq.gz -d QC2 -t 24 -5trim_off -q 20 -split_size 200000 -n 1000
	    #@make_paired_paired_files=&QC;
	    ($before_demultiplex_fastq,$before_demultiplex_barcode)=&join_paired_ends;
	    #$before_demultiplex_mapping .= "$merge_mapping_file," foreach @make_paired_paired_files;
	    #$before_demultiplex_mapping =~ s/,$//;
	}
	else
	{
	    $before_demultiplex_fastq = join(',',@unpaired_files);
	    $before_demultiplex_barcode = join(',',@barcode_files);
	}
	$before_demultiplex_mapping = join(',',@mapping_files);    

	&split_libraries($before_demultiplex_fastq,$before_demultiplex_barcode,$before_demultiplex_mapping);
}

my ($biom,$rep_set_fna,$rep_set_tre)=&pick_otus;

my ($otu_tables, $num_sample,$sampling_depth);

$biom=&filter_by_taxomony($biom) if ($negative_taxa);
$biom=&substract_no_template_control($biom) if ($ntc_list);

($sampling_depth,$num_sample,$otu_tables)=&determine_sampling_depth_and_otu_table($biom);

#$sampling_depth = $sampling_depth_cutoff if ($sampling_depth_cutoff);

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
if ($num_sample>50){
	# only plot area chart when large samples
	`sed -i 's/plot_taxa_summary:chart_type [a-Z,A-Z]*/plot_taxa_summary:chart_type area/' $parameter_file`;
}
if ($num_sample>1)
{
    &core_diversity_analysis($biom,$sampling_depth,$analysis_outputDir,$rep_set_tre,$core_diveristy_options);
    &heatmap($biom,$analysis_outputDir,$rep_set_tre);
    &update_index_html;
}
else
{
    &heatmap($biom,$analysis_outputDir,$rep_set_tre,$merge_mapping_file);
    &summarize_taxa($biom,$analysis_outputDir);
    &alpha_diversity_analysis($biom,$analysis_outputDir,$rep_set_tre,$merge_mapping_file,$sampling_depth);
    &update_index_html_one_sample;
}



system("rm -rf jobs") if (-e "jobs");

## clean up
unless ($debug){
	system("rm -rf $outDir/fastq-join_joined");
	system("rm -rf $outDir/newFastq");
}

&lprint("\nTotal".&printRunTime($^T)."\n");
close LOG;


### END MAIN ###
### Below are subroutins ###
sub alpha_diversity_analysis{
    my $biom = shift;
    my $dir = shift;
    my $otu_tree=shift;
    my $mapping_file=shift;
    my $depth = shift;
    my $cmd = "alpha_rarefaction.py -f -i $biom -t $otu_tree -e $depth -m $mapping_file -o $dir";
    if (-e "$dir/rarefaction.finish")
    {
        &lprint("\nAlpha Rarefaction Analysis finished\n");
        return 0;
    }
    &process_cmd($cmd,"Alpha Rarefaction Analysis ");
    system("touch $dir/rarefaction.finish");
}
sub summarize_taxa
{
    my $biom = shift;
    my $dir = shift;
    my $outputDir="$dir/taxa_summary";
    my $cmd="summarize_taxa.py -i $biom -o $outputDir"; 
    if (-e "$outputDir/taxa.finish")
    {
        &lprint("\nTaxanomy summarized\n");
        return 0;
    }
    &process_cmd($cmd,"Generating taxanomy summary");

    my @matrix = glob("$outputDir/*.txt");
    my $matrix = join(",",@matrix);
    $cmd = "plot_taxa_summary.py -i $matrix -o $outputDir -c bar";
    &process_cmd($cmd,"Generating taxanomy barplot");
    system("touch $outputDir/taxa.finish");
}

sub heatmap
{
    my $biom = shift;
    my $dir = shift;
    my $otu_tree=shift;
    my $mapping_file=shift;
    my $outputFile="$dir/heatmap.pdf";
    my $cmd="make_otu_heatmap.py -i $biom -o $outputFile "; 
    $cmd .= " -t $otu_tree " if (-e "$otu_tree");
    $cmd .= " -m $mapping_file " if ( -e $mapping_file); 
    if (-e "$outputFile")
    {
        &lprint("\nHeatMap generated\n");
        return 0;
    }
    &process_cmd($cmd,"Generating Heatmap");
    system("touch $dir/heatmap.finish");
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
    &lprint("\nGet barcode internal fastq\n");
    foreach my $file (@{$files_ref})
    {
        if (@make_paired_paired_files)
        {
            my ($R1,$R2)=split /\s+/,$file;
            my ($R1_no_barcodes_fastq,$R1_barcodes)=&split_barcode($R1,$barcode_len,$outputDir);
            my ($R2_no_barcodes_fastq,$R2_barcodes)=&split_barcode($R2,$barcode_len,$outputDir);
            my $barcodes=&catBarcodes($R1_barcodes,$R2_barcodes,$outputDir);
            push @return_fastq_files, "$R2_no_barcodes_fastq $R1_no_barcodes_fastq" if ($pe_orientation eq "rf");
            push @return_fastq_files, "$R1_no_barcodes_fastq $R2_no_barcodes_fastq" if ($pe_orientation eq "fr");
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
        print $fh $return_fastq_files[$_],"\n";
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
    my $pid;
    if ( $file=~/\.gz$/i ) { $pid=open($fh, "gunzip -c $file |") or die ("gunzip -c $file: $!"); }
    else { $pid=open($fh,'<',$file) or die("$file: $!"); }
    return ($fh,$pid);
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
    my ($ib1,$pid1)= open_file($b1);
    my ($ib2,$pid2)= open_file($b2);
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
        print $ofh "$id$seq$b2_seq\n$q_id$q_seq$b2_q_seq\n" if ($pe_orientation eq "fr");
        print $ofh "$b2_id$seq$b2_seq\n$b2_q_id$q_seq$b2_q_seq\n" if ($pe_orientation eq "rf");
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
    my ($fh,$pid)= open_file($file);
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
    my $validation_options_r = shift;
    my $validation_options="";
    $validation_options = join (" ", keys %$validation_options_r) if $validation_options_r;
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
        my $cmd="validate_mapping_file.py $validation_options -m $file -o $outputDir";
        &process_cmd($cmd,"Checking Mapping File");
        my $ret=`grep Error $log`;
        if ($ret) 
        {
	    &lprint($ret."\n");
	    `cp $outputDir/combined_mapping_corrected.txt $file` if (-e "$outputDir/combined_mapping_corrected.txt");
        #    exit;
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
        &lprint("\nJoining Paired-end Files Done\n");
        open (my $fh, "$outputDir/join.finished") or die "Cannot open $outputDir/join.finished\n";
        $joined_file = <$fh>;   chomp $joined_file;
        if ($demultiplex_dir){
            foreach my $file_i (0..$#make_paired_paired_files){
                my ($R1,$R2)=split /\s+/,$make_paired_paired_files[$file_i];
		my $basename=basename($R1);
                map { if($id_file->{$_} =~ /$basename/){$id_file->{$_}="$tmp_dir$file_i/fastqjoin.join.fastq"} } keys %{$id_file} ;
            }
        }
	if (!eof($fh)){
            $joined_barcode = <$fh>;  chomp $joined_barcode;
        }
        return($joined_file,$joined_barcode);
    }
    foreach my $file_i (0..$#make_paired_paired_files)
    {
        my ($R1,$R2)=split /\s+/,$make_paired_paired_files[$file_i];   
	if ($demultiplex_dir){
            my $basename=basename($R1);
            map { if($id_file->{$_} =~ /$basename/){$id_file->{$_}="$tmp_dir$file_i/fastqjoin.join.fastq"} } keys %{$id_file} ;
        }
        $R1=Cwd::abs_path("$R1");
        $R2=Cwd::abs_path("$R2");
	my $barcode_options = ($demultiplex_dir)?"":"-b $barcode_files[$file_i]";
        my $cmd="join_paired_ends.py -f $R1 -r $R2 $barcode_options -o $tmp_dir$file_i";
        &process_cmd($cmd,"Joining Paired-end Reads");
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
        &lprint("\nMerging Mapping Files Done\n");
        &check_mapping_files(\@tmp) unless (@demultiplex_fa_files);
        return $merged_mapping_file;
    }
    my ($mapping_files_r, $validation_options_r) = &fix_mapping_file(\@mapping_files);
    @mapping_files = @$mapping_files_r;
    my $cmd="merge_mapping_files.py -o $merged_mapping_file -m ". join( ',' ,@$mapping_files_r);
  #  my $validation_options_r;
 #   my $cmd="merge_mapping_files.py -o $merged_mapping_file -m ". join( ',' ,@mapping_files);
    if (scalar(@mapping_files)>1)
    {
        &process_cmd($cmd,"Merge Mapping Files");
    }
    else
    {
        system("cp @$mapping_files_r[0] $merged_mapping_file");
    }
    &check_mapping_files(\@tmp,$validation_options_r) unless (@demultiplex_fa_files);
    #unlink @$mapping_files_r;
    return $merged_mapping_file;
}

sub fix_mapping_file {
    my $files = shift;
    my @fix_mapping_files;
    my %validate_options;
    foreach my $f (@$files){
        my ($basename,$path,$suffix) = fileparse($f,qr/\.[^.]*/);
        my $fh;
        if ($f =~ /xlsx$/){
            open ($fh, "xlsx2csv -d tab $f | ") or die "Cannot read $f\n";
        }else{
            open ($fh, $f) or die "Cannot open $f\n";
        }
        open (my $ofh, ">$outDir/$basename.txt") or die "Cannot write $outDir/$basename.txt\n";
        push @fix_mapping_files, "$outDir/$basename.txt";
	my @header; #SampleID  BarcodeSequence  LinkerPrimerSequence Description
	my $miss_barcode=0;
	my $miss_primer=0;
	my $miss_desc=0;
        while(<$fh>){
                next if (/^\n/);
                chomp;
                if (/SampleID/){
			@header = split /\t/,$_;
			s/\s+// for (@header);
			if (! /BarcodeSequence/){splice @header, 1, 0, 'BarcodeSequence'; $miss_barcode=1;}
			if (! /LinkerPrimerSequence/){splice @header, 2, 0, 'LinkerPrimerSequence'; $miss_primer=1;}
			if (! /Description/){splice @header, @header, 0, 'Description'; $miss_desc=1;}
			print $ofh join("\t",@header),"\n";
		}elsif(/^#/){
			print $ofh $_,"\n";
		}else{
			my @array = split /\t/,$_;
			$array[0] =~ s/[^a-zA-Z0-9]/\./g;
			if ($miss_barcode){ splice @array, 1, 0, ''; $validate_options{"-b"}=1;}
					else{ $array[1] =~ s/\s+//g;}
			if ($miss_primer){ splice @array, 2, 0, ''; $validate_options{"-p"}=1;}
					else{ $array[2] =~ s/\s+//g;}
			if ($miss_desc){ splice @array, @array, 0, 'NA';}
			print $ofh join("\t",@array),"\n";
		}
	}
	close $fh;
	close $ofh;
    }
    return (\@fix_mapping_files, \%validate_options);
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
        &lprint("\nDe-multiplexing Done\n");
        return 0;
     }
     my $offset = ($phred_offset)?"--phred_offset $phred_offset":"";
     my $cmd;
     $cmd="split_libraries_fastq.py -o $outputDir -i $fastq_file -b $barcode_file -m $mapping_file --barcode_type $barcode_len -q $quality_cutoff -p $min_reads_q_fraction -n $seq_max_n $offset";
     if ($demultiplex_dir){
         $fastq_file = join(",",values %{$id_file});
         my $sample_ids = join(",",keys %{$id_file});
	 $cmd="split_libraries_fastq.py -o $outputDir -i $fastq_file --sample_ids $sample_ids --barcode_type 'not-barcoded' -q $quality_cutoff -p $min_reads_q_fraction -n $seq_max_n $offset";
     }
     $cmd .= "--store_demultiplexed_fastq " if ($UPARSE_opt);
     &process_cmd($cmd,"De-multiplexing samples");
     
     system("touch $outputDir/split.finished");
}


sub pick_otus
{
    my $outputDir="$outDir/otus";
    my $biom="$outputDir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom";
    $biom = "$outputDir/otu_table_mc${min_otu_size}_w_tax.biom" if ($target_ref =~ /ITS/i);
    my $rep_set_fna="$outputDir/rep_set.fna";
    my $rep_set_tre="$outputDir/rep_set.tre";
    my $seq = "$outDir/slout/seqs.fna";
    my $seq_fastq = "$outDir/slout/seqs.fastq";
    my $UPARSE_log = "$outputDir/log.txt";
    if (-e "$outputDir/pickOTU.finished")
    {
        &lprint("\nPicking OTUs Done\n");
        return ($biom,$rep_set_fna,$rep_set_tre);
    }
    #my $cmd="pick_open_reference_otus.py -f -i $seq -r $reference_seqs -o $outputDir -s 0.01 $parameter_file --min_otu_size $min_otu_size -aO $CPUs";
    my $cmd;
    $cmd="pick_open_reference_otus.py -f -i $seq -o $outputDir -s 0.01 $pick_otu_ref $parameter_file --min_otu_size $min_otu_size -aO $CPUs";
    $cmd="UPARSE_pick_otu.pl -i $seq_fastq -min_otu_size $min_otu_size -t $CPUs -similarity $similarity -o $outputDir 1>$UPARSE_log 2>\&1 " if ($UPARSE_opt);
    &process_cmd($cmd,"Pick OTUs");
    system("touch $outputDir/pickOTU.finished");
    system("rm -rf $outputDir/tmp");
    system("rm -rf $outputDir/step*otus");
    return($biom,$rep_set_fna,$rep_set_tre);
}

sub filter_by_taxomony {
	my $biom=shift;
	my $output = "$outDir/otus/otu_table_filter_by_taxonomy.biom";
	my $cmd="filter_taxa_from_otu_table.py -i $biom -o $output -p $negative_taxa";
	if (-e "$output")
    	{
        	&lprint("\nFilter by Taxonomy Done\n");
		return $output;
    	}
	&process_cmd($cmd,"Filter by Taxonomy");
	return $output;
}

sub substract_no_template_control {
	my $biom=shift;
	my $biom_summary="$outDir/otus/biom_table_summary_with_ntc.txt";
	my $otu_table="$outDir/otus/otu_table_tabseparated.txt";
	my $new_otu_table="$outDir/otus/otu_table_tabseparated_no_ntc.txt";
	my $old_otu_table="$outDir/otus/otu_table_tabseparated_with_ntc.txt";
	my $new_biom = "$outDir/otus/otu_table_substract_by_ntc.biom";
	my $otu_table_cmd="biom convert -i $biom -o $otu_table --to-tsv --header-key=\"taxonomy\"";
	unlink ($otu_table,$new_otu_table,$new_biom);
	&process_cmd($otu_table_cmd,"Generating OTUs table");
	my $cmd = "substract_no_template_control.pl -i $otu_table -ntc $ntc_list -o $new_otu_table";
	&process_cmd($cmd,"Substract No Template Control");
	$cmd = "biom convert -i $new_otu_table -o $new_biom --to-hdf5 --table-type=\"OTU table\" --process-obs-metadata taxonomy";
	&process_cmd($cmd,"Convert OTU Table to Biom");	
	$cmd="biom summarize-table -i $biom -o $biom_summary ";
	&process_cmd($cmd,"Generating Biom Table Summmary") if (! -e $biom_summary);
	$cmd="biom summarize-table --qualitative -i $biom -o $biom_summary.otu ";
	&process_cmd($cmd,"Generating Biom Table Summmary") if (! -e "$biom_summary.otu");
	system("mv $otu_table $old_otu_table");
	system("mv $new_otu_table $otu_table");
	return ($new_biom);
}

sub determine_sampling_depth_and_otu_table
{
    my $biom=shift;;
    my $biom_summary="$outDir/otus/biom_table_summary.txt";
    my $otu_tables="$outDir/otus/otu_table_tabseparated.txt";
    my $sampling_depth;
    my $num_sample;
    my $cmd="biom summarize-table -i $biom -o $biom_summary ";
    my $otu_table_cmd="biom convert -i $biom -o $otu_tables --to-tsv --header-key=\"taxonomy\"";
    &process_cmd($otu_table_cmd,"Generating OTUs table") if (! -e $otu_tables);
    &process_cmd($cmd,"Generating Biom Table Summmary") if (! -e $biom_summary);
    $cmd="biom summarize-table --qualitative -i $biom -o $biom_summary.otu ";
    &process_cmd($cmd,"Generating Biom Table Summmary") if (! -e "$biom_summary.otu");
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
	if ( $value > $sampling_depth_cutoff){
		$pass_cutoff++;
        	push @sample_counts, int($value);
	}
    }
    close $fh;
    if ($pass_cutoff == 0 ){
	my $msg = "ERROR: No sample size  > $sampling_depth_cutoff. Stop process. Please adjust sampling depth based on file otus/". &basename($biom_summary)."\n"; 
	&lprint($msg);
	exit(1);
    }
    $sampling_depth=  (sort {$a<=>$b} @sample_counts)[0];
    $sampling_depth = $sampling_depth - 1; 
    #my ($mean_depth,$std,$median_depth,$q1,$q3) = &statistical_calculation_on_array(@sample_counts);
    #$sampling_depth = $mean_depth - 2*$std;
    #if ($mean_depth < 100)
    #{
	#die "Average sample size < 100. Too small to proceed. See file $biom_summary\n";
    #}
    #if ($sampling_depth<100)
    #{
    #    my $outlier_cutoff = $q1 - 1.5*($q3-$q1);
    #    $sampling_depth = ($outlier_cutoff>100)? $outlier_cutoff:$q1-100;
    #}
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
    my $biom=shift;
    my $sampling_depth=shift;
    my $outputDir=shift;
    my $rep_set_tre=shift;
    my $options=shift;
    my $DieCatch=1;
    if (-e "$outputDir/diversity.finished")
    {
        &lprint("\nDiversity Analysis Done\n");
        return 0;
    }
    my ($mapping_file_r,$categories)=&get_catetory_from_mapping_file($merge_mapping_file);
    $categories = "-c $categories " if ($categories);
    #my $cmd = "core_diversity_analyses.py --recover_from_failure -i $biom -o $outputDir -m $merge_mapping_file -e $sampling_depth -t $reference_tree $parameter_file -aO $CPUs";
    my $cmd = "core_diversity_analyses.py $options -i $biom -o $outputDir -m $merge_mapping_file -e $sampling_depth $categories $parameter_file -aO $CPUs";
    $cmd .= " -t $rep_set_tre " if ($options !~ /nonphylogenetic_diversity/i);
    my $return=&process_cmd($cmd,"Perform Diversity Analysis and Taxanomy Summary",$DieCatch);
    if ($return)
    {
	$cmd .= " -w > core_analysis_commands.sh";
	&process_cmd($cmd,"Perform Diversity Analysis and Taxanomy Summary");
	&process_cmd("sh core_analysis_commands.sh", "Skip error step and run remaining analysis");
	unlink "core_analysis_commands.sh";
    }
    system("touch $outputDir/diversity.finished");
    system("cp $outputDir/index.html $outputDir/index.html.org");
}

sub get_catetory_from_mapping_file{
	# structur
	# $hash->{sampleID}->{feature}= Value
	my $mapping_file=shift;
	my %hash;
	my %unique;
	my $sample_num=0;
	my $ntc_sample_num=0;
	open (my $fh, $mapping_file) or die "Cannot read $mapping_file $!\n";
	my @header;
	while(<$fh>){
		chomp;
		@header = split /\t/,$_ if (/SampleID/);
		next if (/^#/);	
		next if (/^\n/);	
		my $ntc=0;
		my @array = split /\t/,$_;
		for my $i (1..$#array){
			if ($header[$i] !~ /SampleID|Barcode|Linker|Description/){
				$unique{$header[$i]}->{$array[$i]}++ if ( $array[$i] ne 'NTC' && $array[$i] !~ /no template/i );
				$ntc=1 if ( $array[$i] eq 'NTC' || $array[$i] =~ /no template/i );
			}
			$hash{$array[0]}->{$header[$i]}=$array[$i];
		}
		$sample_num++;
		$ntc_sample_num++ if (  $ntc );
	}
	close $fh;
	my @category_for_analysis;
	foreach my $feature (keys %unique){
		my $unique_feature_num = scalar (keys %{$unique{$feature}});
		push @category_for_analysis, $feature if ($unique_feature_num > 1 && $unique_feature_num < ($sample_num - $ntc_sample_num));
	}
	my $categories = join (",",@category_for_analysis);
	return (\%hash,$categories);
}

sub update_index_html_one_sample
{
    my $index = "$outDir/otus/index.html";
    system ("mv $outDir/checkMappingFile $outDir/analysis/checkMappingFile");
    system ("cp $outDir/combined_mapping.txt $outDir/analysis/combined_mapping.txt");
    copy($index,"$index.org");
    open (my $fh, "<", "$index.org") or die "Cannot open $index.org\n";
    open (my $ofh, ">","$index") or die "Cannot write $index\n";
    while(<$fh>)
    {
        chomp;
        if (/^<table/){
            print $ofh "<h1>Project: $title</h1>\n";
            print $ofh $_,"\n";
        }
        elsif(/td>Run summary/){
            print $ofh $_;
            if ($negative_taxa){
                print $ofh "<tr><td>Mapping File Check</td><td> <a href=\"../analysis/combined_mapping.txt\" target=\"_blank\">combined_mapping.txt</a></td></tr>\n";
	    }else{
                print $ofh "<tr><td>Mapping File Check</td><td> <a href=\"../analysis/checkMappingFile/combined_mapping.html\" target=\"_blank\">combined_mapping.html</a></td></tr>\n";
	    }
            print $ofh "<tr><td>OTUs HeatMap</td><td> <a href=\"../analysis/heatmap.pdf\" target=\"_blank\">heatmap.pdf</a></td></tr>\n";
        }elsif(/Taxonomy assignments/){
		print $ofh $_,"\n";
		print $ofh "<tr><td>Bar plot</td><td><a href=\"../analysis/taxa_summary/bar_charts.html\" target=\"_blank\">bar_charts.html</a></td></tr>\n";
	}elsif(/er>Trees/){
		print $ofh "<tr><td>BIOM table statistics</td><td> <a href=\"./biom_table_summary.txt\" target=\"_blank\">biom_table_summary.txt</a></td></tr>\n";
		print $ofh $_,"\n";
        }elsif(/\<\/table/){
		print $ofh "<tr colspan=2 align=center bgcolor=#e8e8e8><td colspan=2 align=center>Alpha diversity results</td></tr>\n";
		print $ofh "<tr><td>Rarefaction plots</td><td><a href=\"../analysis/alpha_rarefaction_plots/rarefaction_plots.html\" target=\"_blank\">rarefaction_plots.html</a></td></tr>\n";
            	print $ofh $_,"\n";
	}else{
            print $ofh $_,"\n";
	}
    }
    close $fh;
    close $ofh;
}

sub update_index_html
{
    my $index = "$outDir/analysis/index.html";
    system ("rm -rf $outDir/analysis/checkMappingFile");
    system ("mv $outDir/checkMappingFile $outDir/analysis/checkMappingFile");
    system ("cp $outDir/combined_mapping.txt $outDir/analysis/combined_mapping.txt");
    open (my $fh, "$index.org") or die "Cannot open $index.org\n";
    open (my $ofh, ">$index") or die "Cannot write $index\n";
    while(<$fh>)
    {
        chomp;
        if (/^<table/)
        {
            print $ofh "<h1>Project: $title</h1>\n";
            print $ofh $_,"\n";
        }
        elsif(/Master/)
        {
            print $ofh $_;
            if ($negative_taxa){
                print $ofh "<tr><td>Mapping File Check</td><td> <a href=\"./combined_mapping.txt\" target=\"_blank\">combined_mapping.txt</a></td></tr>\n";
	    }else{
                print $ofh "<tr><td>Mapping File Check</td><td> <a href=\"./checkMappingFile/combined_mapping.html\" target=\"_blank\">combined_mapping.html</a></td></tr>\n";
	    }
            print $ofh "<tr><td>Representative Sequences</td><td> <a href=\"./rep_set.fna\" target=\"_blank\">rep_set.fna</a></td></tr>\n";
            print $ofh "<tr><td>Representative Sequences Tree</td><td> <a href=\"./rep_set.tre\" target=\"_blank\">rep_set.tre</a></td></tr>\n";
            print $ofh "<tr><td>OTUs Table</td><td> <a href=\"./OTUs.table.txt\" target=\"_blank\">OTUs.table.txt</a></td></tr>\n";
            print $ofh "<tr><td>Biom Table</td><td> <a href=\"./table.biom.gz\" target=\"_blank\">table.biom.gz</a></td></tr>\n";
            print $ofh "<tr><td>OTUs HeatMap</td><td> <a href=\"./heatmap.pdf\" target=\"_blank\">heatmap.pdf</a></td></tr>\n";
        }
        else
        {
            if ($_ =~ /bar_charts|area_charts/){
                 my ($file)= $_ =~ /href=\"\.(\S+)\"\s/; 
                 $file = "$outDir/analysis$file";
                 next if ( ! -e "$file");
            }
            print $ofh $_,"\n";
        }
    }
    close $fh;
    close $ofh;
}

sub parse_mapping_files{
    my $dir=shift;
    my $mapping_files=shift;
    my @pe_files;
    my @se_files;
    if ( ! -d $dir ){
        my $msg = "ERROR: the input directroy does not exist or isn't a directory\n";
        &lprint($msg);
        exit(1);
    }
    foreach my $f (@{$mapping_files}){
        my $file_column_index;
        my $fh;
	if ($f =~ /xlsx$/){
            open ($fh, "xlsx2csv -d tab $f | ") or die "Cannot read $f\n";
        }else{
            open ($fh, $f) or die "Cannot read $f\n";
        }
	while(<$fh>){
            chomp;
            next if (/^\n/);
            next unless (/\S/);
            if (/SampleID/){
                my @header = split /\t/,$_;
                ( $file_column_index )= grep { $header[$_] =~ /files/i } 0..$#header;
            }elsif(! /^#/){
                my @array = split /\t/,$_;
		$array[0] =~ s/[^a-zA-Z0-9]/\./g;
		$array[$file_column_index] =~ tr/"//d;
                $id_file->{$array[0]}=$array[$file_column_index];
		my @files = map { "$dir/$_" } split /,|\s+/,$array[$file_column_index];
		$id_file->{$array[0]}="$dir/$array[$file_column_index]" if (scalar(@files)==1);
                if (scalar(@files) % 2){
                    push @se_files,@files;
                }else{
                    push @pe_files,@files;
                }
            }
        }
        close $fh;
    }
    return (\@pe_files,\@se_files);
}

sub process_cmd {
    my ($cmd, $msg, $dieCatch) = @_;


    if ($msg) {
        &lprint("\n\n");
        &lprint("###########################################################################\n");
        &lprint("Qiime ".&getTimeString."  $msg\n");
        &lprint("###########################################################################\n");
    }
    
    &lprint("Qiime CMD: $cmd\n");
    if ($debug) {
        print STDERR "\n\n-WAITING, PRESS RETURN TO CONTINUE ...";
        #my $wait = <STDIN>;
        print STDERR "executing cmd.\n\n";
    }
    

    my $time_start = time();
    
    my $ret = system($cmd);
    #my $time_end = time();

    if ($ret) {
	if ($dieCatch)
        {
	    &lprint($ret);
        }else
	{
	    die "Error, CMD: $cmd died with ret $ret";
        }
    }
 
    #my $number_minutes = sprintf("%.1f", ($time_end - $time_start) / 60);
    &lprint("\nQiime".&printRunTime($time_start)."\n");
 #   print "TIME: $number_minutes min. for $cmd\n";
    #&lprint("TIME: $number_minutes min.\n");
    

    return $ret;
}

sub printRunTime {
  my $time=shift;
  my $runTime = time() - $time;
  my $time_string =sprintf(" Running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime % 3600) / 60), 
  int($runTime % 60));
  return $time_string;
}

sub getTimeString
{
    my $now_string = strftime "[%Y %b %e %H:%M:%S]", localtime;
    #$now_string =~ s/\s//g;
    return $now_string;
}

sub lprint {
      my ($line) = @_;
      print LOG $line;  
      print $line;
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


sub is_fastq
{
    $SIG{'PIPE'}=sub{};
    my $file=shift;
    my ($fh,$pid)= open_file($file);
    my $head=<$fh>;
    close $fh; 
    kill 9, $pid; # avoid gunzip broken pipe
    
    $SIG{'PIPE'} = 'DEFAULT';
    ($head =~/^@/)?
        return 1:
        return 0;
}
