#! /usr/bin/env perl
# required: 1. R
#           2. samtools >1.1
#           3. bwa 0.6 
#           4. bowtie2
#           5. bcftools > 1.1 
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
# 20180123 use samtools 1.6 and bcftools 1.6 

use Getopt::Long;
use File::Basename;
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use fastq_utility;
use Parallel::ForkManager;

my $debug=0;

$|=1;
my ($file1, $file2, @ref_files, $outDir, $pacbio, $offset);
my $file_long="";
my $singleton="";
my $paired_files="";
my $numCPU=4;
my $bwa_options ="-t $numCPU ";
my $minimap2_options ="-t $numCPU ";
my $bowtie_options ="-p $numCPU -a ";
my $snap_options ="-t $numCPU -M ";
my $cov_cut_off=80;
my $aligner="bwa";
my $gen_consensus=1;
my ($window_size, $step_size);
#my ($window_size, $step_size)=(1000,200);
my $pacbio_bwa_option="-b5 -q2 -r1 -z10 "; 
my $prefix="ReadsMapping";
my $plot_only=0;
my $skip_aln=0;
my $no_plot=0;
my $no_snp=0;
my $min_indel_candidate_depth=3;  #minimum number gapped reads for indel candidates
# varinat filter
my $min_alt_bases=3;  # minimum number of alternate bases
my $min_alt_ratio=0.3; #  minimum ratio of alternate bases
my $max_depth=1000000; # maximum read depth
my $min_depth=7; #minimum read depth
my $snp_gap_filter=3; #SNP within INT bp around a gap to be filtered


$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";
 
GetOptions( 
            'aligner=s' => \$aligner,
            'p=s'       => \$paired_files,
            'ref=s{,}' => \@ref_files, # reference/contigs file
            'pre=s' => \$prefix,
            'long=s' =>  \$file_long,
            'u=s' => \$singleton, # illumina singleton 
#            'window_size=i' => \$window_size,  # for coverage plot
#            'step_size=i' => \$step_size,  # for coverage plot
            'd=s'   => \$outDir,
            'bwa_options=s' => \$bwa_options,
            'bowtie_options=s' => \$bowtie_options,
            'snap_options=s'  => \$snap_options,
            'minimap2_options=s'  => \$minimap2_options,
            'pacbio' => \$pacbio,
            'consensus=i' => \$gen_consensus,
            'min_indel_candidate_depth' => \$min_indel_candidate_depth,
            'min_alt_bases' => \$min_alt_bases,
	    'min_alt_ratio' => \$min_alt_ratio,
            'max_depth' => \$max_depth,
            'min_depth' => \$min_depth,
            'snp_gap_filter' => \$snp_gap_filter,
            'cpu=i' => \$numCPU,
            'plot_only' => \$plot_only,
            'skip_aln'  => \$skip_aln,
            'no_plot'   => \$no_plot,
            'no_snp'    => \$no_snp,
            'debug'     => \$debug,
            'help|?',  sub {Usage()}
);

my $tmp = "$outDir/tmp";

## input check ##
unless (@ref_files) { &Usage("No Reference files.");}
map{ if ( ! -e $_ ){&Usage("The reference file $_ not exist.");} } @ref_files;
unless ( $outDir) { &Usage("No output directory specified.");}
unless ( $paired_files or -e $file_long or -e $singleton) { &Usage("No Reads input."); }
if ($paired_files){
  ($file1, $file2) = split /\s+/,$paired_files;
  unless (-e $file1 && -e $file2) {print "$file1 or $file2 not exists\n";&Usage();}
}

`mkdir -p $tmp`;

#if ($step_size > $window_size) {die "The step_size ($step_size) should be less than window_size ($window_size)\n&Usage";}

## output file variable initialized ##
my $final_stats_output="$outDir/$prefix.alnstats.txt";
my $plotsPdf="$outDir/${prefix}_plots.pdf";
#my $final_vcf_output="$outDir/$prefix.vcf";
#my $bcf_output="$outDir/$prefix.raw.bcf";
#my $final_bam_output="$outDir/$prefix.sort.bam";
#my $final_bam_index_output="$outDir/$prefix.sort.bam.bai";
#my $pileup_output="$outDir/$prefix.pileup";
#my $ref_window_gc="$outDir/$prefix.ref_windows_gc.txt";
#my $final_consensusSeq="$outDir/${prefix}.consensus.fasta";
#unlink $final_bam_output if (!$plot_only);
#unlink $final_bam_index_output if (!$plot_only);
unlink $final_stats_output;

if (! -e $outDir)
{
     mkdir $outDir;
}

if ($bwa_options =~ /-t\s+\d+/) { $bwa_options =~ s/-t\s+\d+/-t $numCPU/; } else { $bwa_options .= " -t $numCPU ";}
if ($minimap2_options =~ /-t\s+\d+/){$minimap2_options =~ s/-t\s+\d+/-t $numCPU/;}else{$minimap2_options .= " -t $numCPU ";};
if ($bowtie_options =~ /-p\s+\d+/){$bowtie_options =~ s/-p\s+\d+/-p $numCPU/ ;}else{$bowtie_options .= " -p $numCPU ";}
if ($snap_options =~ /-t\s+\d+/){$snap_options =~ s/-t\s+\d+/-t $numCPU/;}else{$snap_options .= " -t $numCPU ";}
my $samtools_threads=$numCPU;

my @bam_outputs;
for my $ref_file_i ( 0..$#ref_files ){
	my $ref_file=$ref_files[$ref_file_i];	
	my ($ref_file_name, $ref_file_path, $ref_file_suffix)=fileparse("$ref_file", qr/\.[^.]*/);
	my $bam_output = "$outDir/$ref_file_name.sort.bam";
	my $bam_index_output = "$outDir/$ref_file_name.sort.bam.bai";
	push @bam_outputs, $bam_output;
	unless ($plot_only){  # skip the alignment steps, SNP steps, assume bam and pileup files were generated.
	unless ($skip_aln){ # skip the alignment steps

	# index reference
	if ( $aligner =~ /bowtie/i and ! -e "$ref_file.1.bt2"){
		# fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
		$ref_file=&fold($ref_file);
		`bowtie2-build $ref_file $ref_file`;
	}
	elsif ($aligner =~ /bwa/i and ! -e "$ref_file.bwt"){
    		# fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    		$ref_file=&fold($ref_file);
    		`bwa index $ref_file 2>/dev/null`; 
	}
	elsif ($aligner =~ /snap/i ){
    		# fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    		$ref_file=&fold($ref_file);
    		`snap-aligner index $ref_file $ref_file.snap -bSpace `;
	}
	elsif ($aligner =~ /minimap2/i ){
    		# fold sequence in 100 bp per line (samtools cannot accept > 65535 bp one line sequence)
    		$ref_file=&fold($ref_file);
    		`minimap2 -d $tmp/$ref_file_name.mmi $ref_file `;
	}
	$ref_files[$ref_file_i]=$ref_file;
	## index reference sequence 
	`samtools faidx $ref_file`;

	if ($file_long)
	{
   		print "Mapping long reads to $ref_file_name\n";
   		if ($aligner =~ /bowtie/i){
     			`bowtie2 -a --local $bowtie_options -x $ref_file -fU $file_long -S $outDir/LongReads$$.sam`;
   		}
   		elsif($aligner =~ /bwa/i)
  		{
     			if ($pacbio){
				# `bwa bwasw -M -H $pacbio_bwa_option -t $bwa_threads $ref_file $file_long -f $outDir/LongReads$$.sam`;
				`bwa mem -x pacbio $bwa_options $ref_file $file_long | samtools view -@ $samtools_threads -ubS - | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/LongReads$$.bam - `;
     			}
			else{
				#`bwa bwasw -M -H -t $bwa_threads $ref_file $file_long -f $outDir/LongReads$$.sam`;
				`bwa mem $bwa_options $ref_file $file_long | samtools view -@ $samtools_threads -ubS - |  samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/LongReads$$.bam`;
			}
		#my $mapped_Long_reads=`awk '\$3 !~/*/ && \$1 !~/\@SQ/ {print \$1}' $tmp/LongReads$$.sam | uniq - | wc -l`;
  		#`echo -e "Mapped_reads_number:\t$mapped_Long_reads" >>$outDir/LongReads_aln_stats.txt`;
		}
		elsif ($aligner =~ /snap/i){
			`snap-aligner single $ref_file.snap $file_long -o $outDir/LongReads$$.sam $snap_options`;
		}
		elsif ($aligner =~ /minimap2/i){
			`minimap2 -La $minimap2_options  $tmp/$ref_file_name.mmi $file_long > $outDir/LongReads$$.sam`;
		}


		`samtools view -@ $samtools_threads -t $ref_file.fai -uhS $outDir/LongReads$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/LongReads$$.bam` if ( -s "$outDir/LongReads$$.sam");
	}
	if ($paired_files){
		print "Mapping paired end reads to $ref_file_name\n";
		$offset = fastq_utility::checkQualityFormat($file1);
		my $quality_options="";
		if ($aligner =~ /bowtie/i){
			$quality_options = " --phred64 " if ($offset==64);
			`bowtie2 $bowtie_options $quality_options -x $ref_file -1 $file1 -2 $file2 -S $outDir/paired$$.sam`;
		}
		elsif ($aligner =~ /bwa_short/i){
			$quality_options = " -I " if ($offset==64);
			`bwa aln $bwa_options $quality_options $ref_file $file1 > $tmp/reads_1_$$.sai`;
			`bwa aln $bwa_options $quality_options $ref_file $file2 > $tmp/reads_2_$$.sai`;
			`bwa sampe -a 100000 $ref_file $tmp/reads_1_$$.sai $tmp/reads_2_$$.sai $file1 $file2 > $outDir/paired$$.sam`;
		}
		elsif ($aligner =~ /bwa/i){
			`bwa mem $bwa_options $ref_file $file1 $file2 | samtools view -@ $samtools_threads -ubS -| samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/paired$$.bam `;
		}
		elsif ($aligner =~ /snap/i){
			`snap-aligner paired $ref_file.snap $file1 $file2 -o $outDir/paired$$.sam $snap_options`;
		}
		elsif ($aligner =~ /minimap2/i){
			`minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $file1 $file2 > $outDir/paired$$.sam`;
		}
		`samtools view -@ $samtools_threads -t $ref_file.fai -uhS $outDir/paired$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/paired$$.bam` if (-s "$outDir/paired$$.sam");
	}

	if ($singleton){
		print "Mapping single end reads to $ref_file_name\n";
		$offset = fastq_utility::checkQualityFormat($singleton);
		my $quality_options="";
		if ($aligner =~ /bowtie/i){
			$quality_options = " --phred64 " if ($offset==64);
			`bowtie2 $bowtie_options $quality_options -x $ref_file -U $singleton -S $outDir/singleton$$.sam`;
		}
		elsif($aligner =~ /bwa_short/i){
			$quality_options = " -I " if ($offset==64);
			`bwa aln $bwa_options $quality_options $ref_file $singleton > $tmp/singleton$$.sai`;
			`bwa samse -n 50 $ref_file $tmp/singleton$$.sai $singleton > $outDir/singleton$$.sam`;
		}
		elsif ($aligner =~ /bwa/i){
			`bwa mem $bwa_options $ref_file $singleton |samtools view -@ $samtools_threads -ubS - | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/singleton$$.bam `;
		}
		elsif($aligner =~ /snap/i){
			`snap-aligner single $ref_file.snap $singleton -o $outDir/singleton$$.sam $snap_options`;
		}
		elsif ($aligner =~ /minimap2/i){
			`minimap2  $minimap2_options -ax sr $tmp/$ref_file_name.mmi $singleton> $outDir/singleton$$.sam`;
		}
		`samtools view -@ $samtools_threads -t $ref_file.fai -uhS $outDir/singleton$$.sam | samtools sort -T $tmp -@ $samtools_threads -O BAM -o $outDir/singleton$$.bam` if (-s "$outDir/singleton$$.sam");
	}

	# merge bam files if there are different file type, paired, single end, long..
	if ($file_long and $paired_files and $singleton){
		`samtools merge -f -h $outDir/paired$$.bam -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam $outDir/LongReads$$.bam`;
	}
	elsif($file_long and $paired_files){
		`samtools merge -f -h $outDir/paired$$.bam -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/LongReads$$.bam`;
	}
	elsif($paired_files and $singleton){
		`samtools merge -f -h $outDir/paired$$.bam -@ $samtools_threads $bam_output $outDir/paired$$.bam $outDir/singleton$$.bam`;
	}
	elsif($singleton and $file_long){
		`samtools merge -f -h $outDir/singleton$$.bam -@ $samtools_threads $bam_output $outDir/singleton$$.bam $outDir/LongReads$$.bam`;
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

	} # unless ($skip_aln);


	} # unless ($plot_only)
} # foreach $ref_file



my $Rscript = "$tmp/Rscript$$";
open (my $pdf_fh, ">$Rscript") or die "Cannot write $Rscript\n";
print $pdf_fh "pdf(file=\"$plotsPdf\",width=10,height=8); \n";
my $stats_print_string_head = "Ref\tRef_len\tRef_GC%\tMapped_reads\tRef_recovery%\tAvg_fold(x)\tFold_std\tNum_of_Gap\tTotal_Gap_bases";
$stats_print_string_head .= "\tNum_of_SNPs\tNum_of_INDELs" if (!$no_snp);
`echo  "$stats_print_string_head" > $final_stats_output`;

my $pm = new Parallel::ForkManager($samtools_threads);

$pm -> run_on_finish ( # called BEFORE the first call to start()
		sub {
        		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
			
			if (defined($data)) {  # children are not forced to send anything
				print $pdf_fh $data->{pdfRscript},"\n";
				`echo "$data->{stats_print_string}" >> $final_stats_output`;
			}
		}
	);

for my $ref_file_i ( 0..$#ref_files){
	my $hash_ref;
	my $ref_file = $ref_files[$ref_file_i];
	my ($ref_file_name, $ref_file_path, $ref_file_suffix)=fileparse("$ref_file", qr/\.[^.]*/);
	$pm->start($ref_file_name) and next;
	my $bam_output = "$outDir/$ref_file_name.sort.bam";
	my $bam_index_output = "$outDir/$ref_file_name.sort.bam.bai";
	my $pileup_output="$outDir/$ref_file_name.pileup";
	my $bcf_output = "$outDir/$ref_file_name.raw.bcf";
	my $vcf_output="$outDir/$ref_file_name.vcf";
	my $stats_output="$outDir/$ref_file_name.alnstats.txt";
	my $consensusSeq="$outDir/$ref_file_name.consensus.fasta";
	my $ref_window_gc="$outDir/$ref_file_name.windows_gc.txt";
	
	unless ($plot_only){ # skip the alignment steps, SNP steps, assume bam and pileup files were generated.
		## SNP call
		if (!$no_snp){
			print "SNPs/Indels call on $ref_file_name\n";
			`bcftools mpileup -d $max_depth -L $max_depth -m $min_indel_candidate_depth -Ov -f $ref_file $bam_output | bcftools call -cO b - > $bcf_output 2>/dev/null`;
			`bcftools view -v snps,indels,mnps,ref,bnd,other -Ov $bcf_output | vcfutils.pl varFilter -a$min_alt_bases -d$min_depth -D$max_depth -r $min_alt_ratio > $vcf_output`;
		}

		## index BAM file 
		`samtools index $bam_output $bam_index_output`; 

		## generate statistical numbers 
		print "Generate alignment statistical numbers \n";
		`samtools flagstat $bam_output > $stats_output`; 
 
		## derived chimera info 
		if ($aligner=~ /bwa/i and $paired_files){ 
			my $proper_paired = `grep "properly paired"  $stats_output | awk '{print \$1}' `;
			my $all_mapped_paired = `grep "with itself and mate mapped"  $stats_output | awk '{print \$1}' `;
			my $chimera = $all_mapped_paired - $proper_paired;
			chomp $chimera;
			`echo  "Chimera:\t$chimera" >>$stats_output`;
		}
	}
	my %proper_base_hash;
	my %unproper_base_hash;
	if (0) { # disable for now
		my $proper_paired_coverage="$outDir/$ref_file_name.proper_paired$$.coverage";
		my $unproper_paired_coverage="$outDir/$ref_file_name.unproper_paired$$.coverage";

		## generate proper-paired reads coverage
		`samtools view -@ $samtools_threads -u -h -f 2 $bam_output | samtools mpileup -BQ0 -d10000000 -f $ref_file - | awk '{print \$1"\\t"\$2"\\t"\$4}'  > $proper_paired_coverage`;

		## generate non-proper-paired reads coverage 2 (properpaired)+4(query unmapped)+8(mate unmapped)
		`samtools view -@ $samtools_threads -u -h -F 14 $bam_output | samtools mpileup -ABQ0 -d10000000 -f $ref_file - | awk '{print \$1"\\t"\$2"\\t"\$4}'  > $unproper_paired_coverage`;
		# build proper_paired mapped reads base coverage hash
		open (my $fh1,"$proper_paired_coverage");
		while (<$fh1>)
		{
			chomp;
			my ($id ,$pos, $cov)=split /\t/ ;
			$proper_base_hash{$id}->{$pos}=$cov;
		}
		close $fh1;

		# build unproper_paired mapped reads base coverage hash
		open (my $fh2,"$unproper_paired_coverage");
		while (<$fh2>)
		{
			chomp;
			my ($id ,$pos, $cov)=split /\t/;
			$unproper_base_hash{$id}->{$pos}=$cov;
		}
		close $fh2;
	}

	## generate genome coverage plots and histograms 
	print "Generate genome coverage plots and histograms...\n";




	# get reference informaiton
	my $num_ref=0;
	my $ref_hash=&get_ref_info($ref_file,$ref_window_gc);
	$ref_hash=&mapped_reads_per_contigs($bam_output,$ref_hash);
	&get_consensus($bcf_output, $ref_hash ,$consensusSeq) if ( -e "$bcf_output" && $gen_consensus);


	system("mkdir -p $outDir/Coverage_plots") if (! $no_plot);
	foreach my $ref_name (sort {$ref_hash->{$b}->{reads} <=> $ref_hash->{$a}->{reads} } keys %{$ref_hash})
	{      
		$num_ref++;
		my ($snp_num , $indel_num);
		my $ref_len = $ref_hash->{$ref_name}->{len};
		my $ref_GC = $ref_hash->{$ref_name}->{GC};
		my $ref_desc = $ref_hash->{$ref_name}->{desc};
		my $mapped_reads = $ref_hash->{$ref_name}->{reads} || "0";
		my $stats_print_string = $ref_name."\t".$ref_len."\t".$ref_GC."\t".$mapped_reads."\t";
		# generate coverage file
		my $coverage_output="$outDir/${prefix}_${ref_name}.coverage";
		my $WindowCoverage_output="$outDir/${prefix}_${ref_name}.window_size_coverage";
		my $gap_output="$outDir/${prefix}_${ref_name}.gap.coords";
		my $coverage_plot="$outDir/Coverage_plots/${prefix}_${ref_name}_base_coverage.png";
		my $histogram="$outDir/Coverage_plots/${prefix}_${ref_name}_coverage_histogram.png";

		
		my $pileup_cmd = "samtools mpileup -A -BQ0 -d10000000 -r $ref_name -f  $ref_file $bam_output ";
		# build base coverage hash
		my @base_array= (0) x $ref_len;
		open (my $pileup_fh,"$pileup_cmd | ") or die "$! no $pileup_cmd";
		while (<$pileup_fh>)
		{
			chomp;
			my ($id ,$pos,$ref_base, $cov, $seq, $qual)=split /\t/;
			$base_array[$pos-1]=$cov;
		}
		close $pileup_fh;

		$stats_print_string .= &window_size_coverage($coverage_output,$WindowCoverage_output,\@base_array,$gap_output,$ref_name,$ref_len);
  
		my $properpair_coverage_output;
		my $unproperpair_coverage_output;
		my $other_coverage_plot;
		#if ($paired_files){
		if (0) { # disable 
			$properpair_coverage_output="$outDir/${prefix}_${ref_name}.p$$.window_size_coverage";
			$unproperpair_coverage_output="$outDir/${prefix}_${ref_name}.up$$.window_size_coverage";
			$other_coverage_plot="$outDir/${prefix}_${ref_name}_coverage_comparison.png";
			&window_size_coverage("",$properpair_coverage_output,\%proper_base_hash,"",$ref_name,$ref_len);
			&window_size_coverage("",$unproperpair_coverage_output,\%unproper_base_hash,"",$ref_name,$ref_len);
		}
		if (!$no_snp){
			($snp_num , $indel_num)= &SNP_INDEL_COUNT("$vcf_output","$ref_name");
			$stats_print_string .= $snp_num ."\t". $indel_num;
		}
		if ($num_ref>1){
			`echo  -e "$stats_print_string" >> $stats_output`;
		}else{
			`echo  -e "\n${stats_print_string_head}\n$stats_print_string" >> $stats_output`;
		}
		#$stats_print_string="";  
		# pdf
		my $R_pdf_script=&plot_coverage($coverage_output,$WindowCoverage_output,$ref_window_gc,$gap_output,$prefix,$ref_name,$ref_desc,"","");
 
		$hash_ref->{pdfRscript} .= $R_pdf_script;
		$hash_ref->{stats_print_string} .= $stats_print_string."\n";
		# png
		&plot_coverage($coverage_output,$WindowCoverage_output,$ref_window_gc,$gap_output,$prefix,$ref_name,$ref_desc,$histogram,$coverage_plot);

		unless ($debug){
			#unlink $WindowCoverage_output;
			#unlink $unproperpair_coverage_output;
 			#unlink $properpair_coverage_output;
			#unlink $coverage_output;
		}
	} #foreach ref segement

	#   if ($num_ref>10) {print "There are more than 10 reference sequences, the covearge plot will only be generated first 10 sequences\n";}
	$pm->finish(0,$hash_ref);
} # foreach my $ref_file

$pm->wait_all_children;

print $pdf_fh "\ntmp<-dev.off()\nquit()\n";
close $pdf_fh;
system ("R --vanilla --slave --silent < $Rscript 2>/dev/null") if (!$no_plot);
unlink "Rplots.pdf" if ( -e "Rplots.pdf");

#clean up
unless ($debug){
	`rm -rf $tmp`;
	`rm -rf $outDir/*$$* $outDir/*window_size_coverage $outDir/*.windows_gc.txt`;
}
### END ###

sub plot_coverage
{
   my $coverge_file = shift;
   my $window_coverage_file= shift;
   my $ref_window_gc = shift;
   my $gap_file=shift;
   my $perfix = shift;
   my $ref_name= shift;
   my $ref_desc=shift;
   my $histogram_png=shift;
   my $coverage_png=shift;
   my $coverage_xlab = ($ref_desc)? $ref_desc:$ref_name;
   my $png_Rscript= "$tmp/Rscript_png$$";
   open (my $png_fh, ">$png_Rscript") or die "Cannot write $png_Rscript\n" if ($histogram_png);
   my $print_string;
   $print_string = "bitmap(file=\"$histogram_png\",width=1024,height=640,units=\"px\")\n" if ($histogram_png);
   $print_string .=  "
# histogram
# read file
a<-read.table(file=\"$coverge_file\")
mean_cov<-mean(a\$V2)
std_cov<-sd(a\$V2)
b<-round (5*std_cov)
#c<-a\$V2[a\$V2<(mean_cov+b)]
#d<-length(a\$V2[a\$V2>=(mean_cov+b)])
reflen<-length(a\$V2)
# for coverage plot
coverage<-sprintf(\"Coverage: %.2f %%\", (length(a\$V2[a\$V2>0])/length(a\$V2))*100)
par(mar=c(5,6,4,2))
if (mean_cov < 1)
{
	hist(a\$V2,main=\"Mapping Reads To Reference ${ref_name}: Coverage Depth Histogram\",xlab=\"Coverage(fold)\",ylab=\'Frequency\')
}else{
	h<-hist(a\$V2,breaks=c(0:round(mean_cov+b),max(a\$V2)),plot=FALSE)
	plot(h\$count[2:length(h\$count)],type=\'h\',lwd=3, col=\'black\',main=\"Mapping Reads To Reference ${ref_name}: Coverage Depth Histogram\",xlab=\"Coverage(fold)\",ylab=\'Frequency\',xaxt=\"n\",xlim=c(1,length(h\$count)))
	x<-seq(1,round(length(h\$count)-std_cov) , round(mean_cov/5))
	axis(1,labels=x,at=x,tick=TRUE)
	axis(1,labels=paste(\">\",round(mean_cov+b)),at=round(mean_cov+b),tick=TRUE,las=2)
}
leg.txt<-paste(\"Average fold: \",format(mean_cov,digit=4),\"sd\", format(std_cov,digit=4));
legend(\"topright\",leg.txt)
";

$print_string .= "\nbitmap(file=\"$coverage_png\",width=1024,height=640,units=\"px\")\n" if ($histogram_png);
$print_string .= "
# coverage plot
# init device
#png(filename=\"$coverage_png\",width=1024,height=640)
#
def.par <- par(no.readonly = TRUE) # get default parameters

# setup plotting area
par(mar=c(5,6,4,2))
#par(mar = c(5, 5, 5, 5), xpd=TRUE, cex.main=1.2, cex.lab=1.2, cex.axis=1.2)
a<-read.table(file=\"$window_coverage_file\")
refGC<-read.table(file=\"$ref_window_gc\")
refGC_coord<-refGC\$V2[refGC\$V1==\"$ref_name\"]
refGC_percetage<-refGC\$V3[refGC\$V1==\"$ref_name\"]
data.gaps<-read.table(file=\"$gap_file\",header=TRUE)
gapBp<-sum(data.gaps\$Length);
gapNum<-length(data.gaps\$Length);

par(fig=c(0,1,0,0.75),mar=c(5, 6, 1, 2),cex.main=1.2)
plot(a\$V1,a\$V2,type=\"l\",col=\"blue\",cex=2,xlab=\"$coverage_xlab\",ylab=\"Coverage (fold)\",main=\"\",xlim=c(0,reflen))
leg.txt<-paste(coverage,\";\", sprintf(\"Average fold: %.2fx +/- %.2f\",mean_cov,std_cov),\";\",sprintf(\"Gaps: %d (%d bp)\",gapNum,gapBp),\";\")
mtext(leg.txt,3,adj=0.05,line=-1,cex=0.9)

# get margin coordiates
pa<-par('usr');
# plot gap regions
if (gapNum > 0){
  for(i in 1:dim(data.gaps)[1]){
    rect(data.gaps[i,1],round(pa[3]),data.gaps[i,2],round(pa[3])+(round(pa[4])-round(pa[3]))/100,col=\"black\",border=NA)
  }
}
par(fig=c(0,1,0.75,1),mar=c(0, 6, 3, 2),new=TRUE,cex.main=1.2)
plot(refGC_coord,refGC_percetage,type=\"l\",cex=2,xaxt=\'n\',xlab=\"\",ylab=\"GC %\",main=\"Mapping Reads To Reference ${ref_name}: Genome Coverage\",xlim=c(0,reflen),bty='n')
mtext(\"Reference GC%\",3,adj=0.05,line=-1,cex=0.8)
par(def.par)#- reset to default
";

if ($histogram_png){
   
  print $png_fh "$print_string\ntmp<-dev.off()\nquit()\n";
  close $png_fh;
  system ("R --vanilla --slave --silent < $png_Rscript 2>/dev/null") if (!$no_plot);
  unlink $png_Rscript;
}
#if ($no_snp == 0)
#{
#  leg.txt<-c(leg.txt,paste(\"# of SNPs: \", $snp_num),paste(\"# of INDELs: \", $indel_num))
#}
    return $print_string;
}

sub mapped_reads_per_contigs {
 	my $bam_output=shift;
 	my $ref_hash_r=shift;
  	my %filter;
	open (my $sam_fh, "samtools view -F 4 $bam_output | ") or die "Cannot read $bam_output";
	while(<$sam_fh>){
		chomp;
		my @samFields=split /\t/,$_;
		my $R1_R2 = 1;
		$R1_R2 = 2 if ($samFields[1] & 128);
		my $unique_id=$samFields[0]."_$R1_R2";
		my $ref_id=$samFields[2];
		$ref_id=~ s/\//_/g;
		next if ($ref_id eq '*');
		next if ($filter{$ref_id}{$unique_id});
		$filter{$ref_id}{$unique_id}++;
		$ref_hash_r->{$ref_id}->{reads}++;
	}
	close $sam_fh;
	return $ref_hash_r;
}

sub get_ref_info 
{
    # Given reference file
    # return hash refernece for id as key and len and GC content.
    my $file=shift;
    my $ref_window_gc_file=shift;
    my %hash;
    my $id;
    my $desc;
    my $seq;
    my $seq_len;
    my $GC_content;
    my $avg_pos;
    open (my $out_fh,">$ref_window_gc_file") or die "$!\n";
    open (my $ref_fh,$file) or die "$!\n";
    while (<$ref_fh>)
    {
       chomp;
       if (/>(\S+)\s*(.*)/)
       {
          if ($seq)
          {
             $seq_len = length $seq;
             my $GC_num = $seq=~ tr/GCgc/GCgc/;
             $GC_content = sprintf ("%.2f",$GC_num/$seq_len*100);
             $hash{$id}->{desc}=$desc;
             $hash{$id}->{len}= $seq_len;
             $hash{$id}->{GC}=$GC_content;
             #$window_size= ($seq_len>1000)? int($seq_len/1000)+10:10;
             $window_size= int($seq_len/500)||2;
             $step_size = int($window_size/5)||1;
             $avg_pos = $window_size/2;
             for (my $i=0; $i<=$seq_len-$window_size;$i=$i+$step_size)
             {
                  my $window_seq=substr($seq,$i,$window_size);  
                  $GC_num = $window_seq=~ tr/GCgc/GCgc/;  
                  $GC_content = $GC_num/$window_size*100;
                  if ($i==0)
                  {
                      print $out_fh $id,"\t",$avg_pos,"\t",$GC_content,"\n";
                  }
                  else
                  {
                      print $out_fh $id,"\t",$avg_pos+$i,"\t",$GC_content,"\n";
                  }

             }

          }
          $id=$1;
	  $desc=$2;
          $id =~ s/\//_/g;
          $seq="";
       } 
       else
       {
         $seq.=$_;
       }
    }
          if ($seq)
          {
             $seq_len = length $seq;
             my $GC_num = $seq=~ tr/GCgc/GCgc/;
             $GC_content = sprintf ("%.2f",$GC_num/$seq_len*100);
             $hash{$id}->{desc}=$desc;
             $hash{$id}->{len}= $seq_len;
             $hash{$id}->{GC}=$GC_content;
             #$window_size= ($seq_len>1000)? int($seq_len/1000)+10:10;
             $window_size= int($seq_len/500)||2;
             $step_size = int($window_size/5)||1;
             $avg_pos = $window_size/2;
             for (my $i=0; $i<=$seq_len-$window_size;$i=$i+$step_size)
             {
                  my $window_seq=substr($seq,$i,$window_size);  
                  $GC_num = $window_seq=~ tr/GCgc/GCgc/;  
                  $GC_content = $GC_num/$window_size*100;
                  if ($i==0)
                  {
                      print $out_fh $id,"\t",$avg_pos,"\t",$GC_content,"\n";
                  }
                  else
                  {
                      print $out_fh $id,"\t",$avg_pos+$i,"\t",$GC_content,"\n";
                  }
             }

          }
    close $ref_fh;
    close $out_fh;
    return \%hash;
}

sub window_size_coverage {
   # given output files names and a hash for each ref name and for each base;
   # output coverage per base/window_size per ref. output gap per ref. 
   # return statistiacl numbers, genome recovery, fold coveage, fold coverage std, gap number, gap total bases.  
   my ($coverage_output,$WindowCoverage_output,$base_array,$gap_output,$ref_name,$ref_len)=@_;
   #$window_size= ($ref_len>1000)? int($ref_len/1000)+10:10;
   $window_size= int($ref_len/500)||2;
   $step_size = int($window_size/5)||1;
   my $pos_cov;
   my $cov_sum;
   my $step=0;
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
   my $covered_base_num=0;
   my @cov_array;
   my $stats_return;
   my $cov_out_fh;
   my $window_cov_out_fh;
   my $gap_fh;
   if ($coverage_output)
   {
      open ($cov_out_fh, ">$coverage_output") or die "$! $coverage_output\n";
      open ($gap_fh, ">$gap_output" ) or die "$! $gap_output\n";
      print $gap_fh "Start\tEnd\tLength\tRef_ID\n";
   }
 #  print $window_size," window\t step ",$step_size,"\n";
   open ($window_cov_out_fh, ">$WindowCoverage_output") or die "$! $WindowCoverage_output\n";
   for (1..$ref_len)
   {
      if ($base_array->[$_ - 1]){
         $pos_cov=$base_array->[$_ - 1];
         if ($coverage_output)
         {
            print $cov_out_fh $_,"\t",$pos_cov,"\n";
            if (@gap_array)
            {
               $gap_length = $gap_array[-1] - $gap_array[0]+1;
               print $gap_fh $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t",$ref_name,"\n";
               $gap_count++;
               $gap_total_len += $gap_length;
               @gap_array=();
            }
            $covered_base_num++;
         }
      }else{
         $pos_cov=0;
         if ($coverage_output)
         {
            print $cov_out_fh $_,"\t",$pos_cov,"\n";
            push @gap_array, $_;
         }
      }
      #push @cov_array,$pos_cov;
      $cov_sum += $pos_cov;
      $step_sum += $pos_cov;
      if (($_ % $step_size)==0)
      {
          push @step_sum, $step_sum;
          $step_sum=0;
      }
      if ($_ == $window_size)
      {
          $step=1;
          $window_sum = $cov_sum;
          print $window_cov_out_fh $avg_pos,"\t",$window_sum/$window_size,"\n";
      }
       
      if ($_ > $window_size){
         $step_sum2 += $pos_cov;
         if (($_-$window_size)%$step_size == 0)
         {
            push @step_sum2, $step_sum2;
            $step_sum2=0;
         }
      }
      if ($_ == ($window_size+$step_size*$step))
      {
          my $previous_step_sum = shift @step_sum;
          my $after_step_sum = shift @step_sum2;
          $window_sum = $window_sum + $after_step_sum - $previous_step_sum;
          $avg_pos = $avg_pos + $step_size; 
          print $window_cov_out_fh $avg_pos,"\t",$window_sum/$window_size,"\n";
          $step++;
      }
   }
   if ($coverage_output)
   {
       if (@gap_array){
           $gap_length = $gap_array[-1] - $gap_array[0]+1;
           print $gap_fh $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t",$ref_name,"\n";
           $gap_count++;
           $gap_total_len += $gap_length;
       }
       my ($std_cov,$avg_cov)= &standard_deviation($base_array);
       my $percent_genome_coverage = sprintf ("%.4f",$covered_base_num/$ref_len*100);
       my $fold = sprintf ("%.2f",$avg_cov);
       my $fold_std = sprintf ("%.2f",$std_cov);
       $stats_return = $percent_genome_coverage."\t".$fold."\t".$fold_std."\t".$gap_count."\t".$gap_total_len."\t";      
       close $cov_out_fh;
   }
   close $gap_fh;
   close $window_cov_out_fh;  
   return ($stats_return);
}

sub standard_deviation {
  my($numbers) = @_;
  # Step 1, find the mean of the numbers
  my $total1 = 0;
  foreach my $num (@$numbers) {
    $total1 += $num;
  }
  my $mean1 = $total1 / (scalar @$numbers);

  # Step 2, find the mean of the squares of the differences
  # between each number and the mean
  my $total2 = 0;
  foreach my $num (@$numbers) {
    $total2 += ($mean1-$num)**2;
  }
  my $mean2 = $total2 / (scalar @$numbers);

  # Step 3, standard deviation is the square root of the
  # above mean
  my $std_dev = sqrt($mean2);
  return ($std_dev,$mean1);
}

sub fold {
    # fold and filter reads length by 200 bp.
    my $file=$_[0];
    my $seq;
    my $seq_name;
    my $seq_desc;
    my $len_cutoff=0;
    my $seq_num;
    my ($file_name, $file_path, $file_suffix)=fileparse("$file", qr/\.[^.]*/);
    my $fold_seq_file="$tmp/${file_name}$file_suffix";
    open (my $in_fh,$file);
    open (my $out_fh,">$fold_seq_file");
    while(<$in_fh>){
      chomp;
      if(/>(\S+)\s*(.*)/)
      {
         if ($seq and length($seq)>$len_cutoff)
         {
           $seq =~ s/ //g;
           $seq =~ s/(.{100})/$1\n/g;
           chomp $seq;
           print $out_fh ">","$seq_name $seq_desc","\n",$seq,"\n";
         }
         $seq_name=$1;
         $seq_desc=$2;
         $seq_name =~ s/\W/\_/g;
         $seq="";
         $seq_num++;
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
         print $out_fh ">","$seq_name $seq_desc","\n",$seq,"\n";
    }
    close $in_fh;
    close $out_fh;
    if ($seq_num<1){die "No sequence in your reference file\n";}
    return ($fold_seq_file);
}

sub SNP_INDEL_COUNT
{
   my  $file=shift;
   my  $ref=shift;
   open (my $fh,$file) or die "$!";
   my $indel_count=0;
   my $SNPs_count=0;
   $ref =~ s/\|/\\\|/g;
   while(<$fh>)
   {
       chomp;
       next if (/^#/);
       next if ($_ !~ /$ref/);
       if (/INDEL/)
       {
           $indel_count++;
       }
       else
       {
           $SNPs_count++;
       }
   }
   close $fh;
   return ($SNPs_count,$indel_count);
}

sub get_consensus
{
  my $rawbcf=shift;
  my $refHash=shift;
  my $outputFile=shift;
  open (my $o_fh, "> $outputFile") or die "$!\n";
  open (my $fh, "bcftools view $rawbcf|") or die "$!\n";
  my ($last_chr, $seq, $qual, $last_pos, @gaps);
  my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
           GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');

  $last_chr = '';
  while (<$fh>) {
        next if (/^#/);
        my @t = split;
        if ($last_chr ne $t[0]) {
          if ($last_chr){
            if ($refHash->{$last_chr}->{len} - $last_pos > 1)
            {
               $seq .= 'N' x ($refHash->{$last_chr}->{len} - $last_pos);
            }
            $seq = &fold_str($seq);
            print $o_fh "\>Consensus_To_Ref_$last_chr\n$seq"; 
          }
          ($last_chr, $last_pos) = ($t[0], 0);
          $seq = $qual = '';
          @gaps = ();
        }
     
        die("[vcf2fq] unsorted input\n") if ($t[1] - $last_pos < 0);
        if ($t[1] - $last_pos > 1) {
          $seq .= 'N' x ($t[1] - $last_pos - 1);
#          $qual .= '!' x ($t[1] - $last_pos - 1);
        }
        if (length($t[3]) == 1 && $t[7] !~ /INDEL/ && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
          my ($ref, $alt) = ($t[3], $1);
          my ($b, $q);
          $q = $1 if ($t[7] =~ /FQ=(-?[\d\.]+)/);
          if ($q < 0) {
		$_ = ($t[7] =~ /AF1=([\d\.]+)/)? $1 : 0;
                $b = ($_ < .5 || $alt eq '.')? $ref : $alt;
                $q = -$q;
          } else {
                $b = $het{"$ref$alt"};
                $b ||= 'N';
          }
#          $b = lc($b);
#          $b = uc($b) if (($t[7] =~ /MQ=(\d+)/ && $1 >= $_Q) && ($t[7] =~ /DP=(\d+)/ && $1 >= $_d && $1 <= $_D));
#          $q = int($q + 33 + .499);
#          $q = chr($q <= 126? $q : 126);
          $seq .= $b;
#          $qual .= $q;
        } elsif ($t[4] ne '.') { # an INDEL not deal with it yet
          push(@gaps, [$t[1], length($t[3])]);
        }
        $last_pos = $t[1];
  }
  close $fh;
  if ($seq)
  {
      if ($refHash->{$last_chr}->{len} - $last_pos >= 1)
      {
          $seq .= 'N' x ($refHash->{$last_chr}->{len} - $last_pos);
      }
      $seq = &fold_str($seq);
      print $o_fh "\>Consensus_To_Ref_$last_chr\n$seq";
  }
}

sub fold_str {
  my ($s) = @_;
  $s =~ s/(.{60})/$1\n/g;
  unless ($s =~ /\n$/){ $s= $s ."\n";}
  return $s;
}

sub Usage 
{
my $msg=shift;
print $msg."\n\n" if $msg;
print <<"END";
Usage: perl $0 
               -p                        'leftSequenceFile rightSequenceFile' 
                                         Space-separated paired-end reads in quote
               -u                        sequenceFile
                                         Provides a file containing single-end reads.
               -long                     long reads file in fasta format.  
                                         --pacbio   <bool> using this flag combined with -long for Pacbio long reads (bwa only) 
               -ref                      reference sequences file in fasta format
               -pre                      output files' prefix (default "ReadsMapping")
               -d                        output directory
               -consensus                <bool> output consensus fasta file (default: on, set 0 to turn off)
               -aligner                  bwa or bowtie or snap or minimap2 (default: bwa)
               -bwa_options <String>     bwa options
                                         type "bwa mem" to see options
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
               -minimap2_options         type "minimap2" to see options
                                         default: "-t 4 "
               -skip_aln                 <bool> skip the alignment steps, assume bam files were generated 
                                         and with proper prefix,outpurDir.  default: off
               -no_plot                  <bool> default: off
               -no_snp                   <bool> default: off
               -debug                    <bool> default: off
               -cpu                      number of CPUs [4]. will overwrite aligner options. 
	
               # Variant Filter parameters
               -min_indel_candidate_depth minimum number gapped reads for indel candidates [3]
               -min_alt_bases            minimum number of alternate bases [3]
               -min_alt_ratio            minimum ratio of alternate bases [0.3]
               -max_depth                maximum read depth [1000000]
               -min_depth                minimum read depth [7]
               -snp_gap_filter           SNP within INT bp around a gap to be filtered [3]
		 

Synopsis:
      perl $0 -p 'reads1.fastq reads2.fastq' -u sinlgeton.fastq -long pyroSeq.fasta -ref reference.fasta -pre ReadsMapping -d /outputPath/

END

               #-window_size              genome coverage plot (default: 1000 bp)
               #-step_size                genome coverage plot (default: 200 bp)
exit;
}

