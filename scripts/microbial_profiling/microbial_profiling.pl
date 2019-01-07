#!/usr/bin/env perl
# Paul Li, B-11  (April 2013)
#

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use POSIX qw(strftime);
use File::Basename;
my $script_dirname = dirname(__FILE__);

$|=1;
$ENV{PATH} = "$Bin:$Bin/../:$Bin/script/:$Bin/bin/:$ENV{PATH}:/cm/shared/apps/uge/8.4.4/bin/lx-amd64/";
my $main_pid = $$;

my %opt;
my $threads=4;
my $res=GetOptions(\%opt,
                   'list|l=s',
                   'setting|s=s',
                   'cpu|c=i',
                   'output|o=s',
                   'debug|d',
                   'help|h|?') || &usage();

if ( $opt{help} ) { &usage(); }
unless ( $opt{setting} ) { print "ERROR: No settings file provided.\n"; &usage(); }

#read settings
my $ini_file = $opt{setting};
my $tools = &restore_settings( $ini_file );
my $pid = $$;

#default settings
$threads  = $tools->{system}->{THREADS};
$threads = $opt{cpu} if defined $opt{cpu};
my $extract  = $tools->{system}->{EXTRACT_NUM};
my $p_outdir = $tools->{system}->{OUTDIR};
$p_outdir    = $opt{output} if defined $opt{output};
$tools->{system}->{OUTDIR} = $p_outdir;
my $p_seqdir = $p_outdir."/".$tools->{system}->{SEQDIR};
my $p_logdir = $p_outdir."/".$tools->{system}->{LOGDIR};
my $p_repdir = $p_outdir."/".$tools->{system}->{REPDIR};
my $rep_top  = $tools->{system}->{REPORT_TOP};
my $hl_list  = $tools->{system}->{HIGHLIGHT_LIST};
my $max_process_num = $tools->{system}->{MAX_PROCESS_NUM};

$rep_top = 5   unless defined $rep_top;
#$threads = 4   unless defined $tools->{system}->{THREADS};
$extract = 0.1 unless defined $tools->{system}->{EXTRACT_NUM};
$tools->{system}->{THREADS}=$threads;

#prepare output directory
`mkdir -p $p_outdir`;
`mkdir -p $p_outdir/script`;
`mkdir -p $p_seqdir`;
`mkdir -p $p_logdir`;
`mkdir -p $p_repdir`;

#output file
my $post_script   = "$p_outdir/script/run_post_script.sh";
my $sqdx_script   = "$p_outdir/script/run_sequedex.sh";
my $fileinfo_out  = "$p_repdir/".$tools->{system}->{FILEINFO_OUT};
my $filelist_out  = "$p_repdir/sequence_list.txt";
my $res_usage_out = "$p_repdir/".$tools->{system}->{RESUSAGE_OUT};
my $summary_out   = "$p_repdir/".$tools->{system}->{SUMMARY_OUT};
my $logfile       = "$p_outdir/$tools->{system}->{LOGFILE}";

#sequedex out script
my $sqdx_fh;
my $post_fh;
open($post_fh, ">$post_script") or die "ERROR: Can't create post-processing script file $post_script: $!\n";
my $info_fh;
unless( -s $fileinfo_out ){
	open($info_fh, ">$fileinfo_out") or die "ERROR: Can't create fileinfo file $fileinfo_out: $!\n";
}
my $list_fh;
unless( -e $opt{list} ){
	open($list_fh, ">$filelist_out") or die "ERROR: Can't create filelist $filelist_out: $!\n";
}

my $log_fh;
if (defined $logfile) {
	$opt{verbose} = 1;
	open($log_fh, ">$logfile") || die "ERROR: Can't create log file $logfile: $!\n";
}

# retrieve input files
my $file_info;
my $cmd;
my $count=1;

my @files = @ARGV;

my @filelist_header;
if( defined $opt{list} ){
	open LIST, $opt{list} || die "Can't open $opt{list}\n";
	while(<LIST>){
		chomp;
		next if /^--/;
		next if /^\s*$/;
		next if /^#/;

		if (/^PREFIX/){
			@filelist_header = split /\t/, $_;
			next;
		}

		if( scalar @filelist_header ){
			my @fields = split /\t/, $_;
			for (my $i=0; $i<=$#filelist_header; $i++){
				$file_info->{$count}->{$filelist_header[$i]} = $fields[$i];
			}
		}
		else{
			push @files, $_;
		}

		$count++;
	}
	close LIST;
}

foreach my $file ( @files ) {
	if ( $file =~ /,/ ){
		$file_info->{$count}->{FASTQPE} = $file;
	}
	else{
		$file_info->{$count}->{FASTQ} = $file;
	}
	$count++;
}
my $num = scalar keys %$file_info;
die "No input FASTQ files." if $num < 1;

&prepSequence($file_info, $tools);

my $out = "";

&_notify("\n[FILE PATH]\n\n");

$out .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\n", "PREFIX", "FASTQ", "FASTQSE", "FASTQPE", "FASTA", "FASTA_EXTRACT";
$out .= sprintf "-------------------------------------------------------------------------------------------\n";

foreach my $idx ( sort {$a<=>$b} keys %$file_info ) {
	$out .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\n",
			$file_info->{$idx}->{PREFIX},
			$file_info->{$idx}->{FASTQ},
			defined $file_info->{$idx}->{FASTQSE} ? $file_info->{$idx}->{FASTQSE} : "",
			defined $file_info->{$idx}->{FASTQPE} ? $file_info->{$idx}->{FASTQPE} : "",
			$file_info->{$idx}->{FASTA},
			$file_info->{$idx}->{FASTA_EXTRACT};
}
&_notify($out);
print $list_fh "$out" unless $opt{list};

&_notify("\n[FASTQ STATS]\n\n");

if( $info_fh && $tools->{system}->{RUN_SEQ_STATS} ){
	foreach my $idx ( sort {$a<=>$b} keys %$file_info )
	{
		die "[FASTQ_STATS] FATAL: Can't read FASTQ: ".$file_info->{$idx}->{FASTQ}."\n" if !-e $file_info->{$idx}->{FASTQ};
		my $file = $file_info->{$idx}->{FASTQ};
		my ($rc,$tl) = &countFastq_exe($file);
		$file_info->{$idx}->{INFO}->{TOL_READS} = $rc;
		$file_info->{$idx}->{INFO}->{TOL_BASES} = $tl;
	}
	
	$out = "";
	$out .= sprintf "%14s%14s%15s%14s%14s%11s%8s  %s\n","TOL_READS","PROCESSED","PLATFORM","TOL_BASES","AVG_LENGTH","AVG_SCORE","OFFSET","DATASET";
	$out .= sprintf "---------------------------------------------------------------------------------------------------------------------\n";
	
	my $info_print = "DATASET\tTOL_READS\tPROCESSED\tPLATFORM\tTOL_BASES\tAVG_LENGTH\tAVG_SCORE\tOFFSET\n";
		
	foreach my $idx ( sort {$a<=>$b} keys %$file_info )
	{
		$out .= sprintf "%14s%15s%14s  %s\n",
				$file_info->{$idx}->{INFO}->{TOL_READS},
				$file_info->{$idx}->{INFO}->{TOL_BASES},
				$file_info->{$idx}->{PREFIX};
	
		$info_print .= sprintf "%s\t%.2f\t%.2f\n",
				$file_info->{$idx}->{PREFIX},
				$file_info->{$idx}->{INFO}->{TOL_READS},
				$file_info->{$idx}->{INFO}->{TOL_BASES};
	}
	&_notify("$out");
	print $info_fh "$info_print";
}
else{
	 &_notify("[FASTQ_STATS] File information exists. Skipping step!\n") unless $info_fh;
	 &_notify("[FASTQ_STATS] Skipping calculating sequence stats!\n") unless $tools->{system}->{CAL_SEQ_STATS};
}

#run tools
&_notify("\n[TOOLS]\n\n");
my @childs;
my %job_ids;
my @toolnames = sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools;

if( $tools->{system}->{RUN_TOOLS} ){
	foreach my $idx ( sort {$a<=>$b} keys %$file_info )
	{
		foreach my $tool ( @toolnames )
		{
			next if $tool eq 'system';
			my $input  = $file_info->{$idx}->{FASTQ};
			my $fnb    = $file_info->{$idx}->{PREFIX};
			my $outdir = "$p_outdir/$idx\_$fnb/$tool";
			my $tool_rep_dir = "$p_repdir/$idx\_$fnb/$tool";
			my $prefix = "$fnb";
			my $log    = "$p_logdir/$fnb-$tool.log";

			print "Tool ($tool) - PID: $$, starting...\n";

			# skip if tool output exists
			if( -e "$tool_rep_dir/.finished" ){
				&_notify("[RUN_TOOL] [$tool] Result exists. Skipping tool $tool!\n");
				&_notify("[RUN_TOOL] [$tool] Running time: 00:00:00\n");
				next;
			}

			# prepare command
			my $code;
			my $time = time;
			my $cmd = $tools->{$tool}->{COMMAND};
			my $qsub_cmd = $tools->{system}->{QSUB_COMMAND};
			$cmd = &param_replace( $cmd, $file_info, $tools, $idx, $tool );
			$qsub_cmd = &param_replace( $qsub_cmd, $file_info, $tools, $idx, $tool );
			if ($tools->{system}->{RUN_TOOL_AS_JOB}){
				unlink $log;
				&_notify("[RUN_TOOL] [$tool] COMMAND: $qsub_cmd\n");
				&_notify("[RUN_TOOL] [$tool] Logfile: $log\n");
				#&_notify("$qsub_cmd -v EDGE_HOME=$ENV{EDGE_HOME} $script_dirname/script/$cmd\n");
				$code = `$qsub_cmd -v EDGE_HOME=$ENV{EDGE_HOME} $script_dirname/script/$cmd 2>&1 `;
				my ($job_id) = $code =~ /Your job (\d+)/;
				#&_notify("$code\n");
				&_notify("[RUN_TOOL] [$tool] cluster info: $code\n");
				if ($job_id){
					$job_ids{$job_id}->{tool}=$tool;
					$job_ids{$job_id}->{time}=$time;
					$job_ids{$job_id}->{outdir}=$outdir;
					$job_ids{$job_id}->{log}=$log;
				}else{
					&_notify("[RUN_TOOL] [$tool] Error occured: qsub error\n");
				}
			}else{
				&_notify("[RUN_TOOL] [$tool] COMMAND: $cmd\n");
				&_notify("[RUN_TOOL] [$tool] Logfile: $log\n");
				$code = system("$cmd > $log 2>&1");
				&_notify("[RUN_TOOL] [$tool] Error occured.\n") if $code;
				my $runningtime = &timeInterval($time);
				&_notify("[RUN_TOOL] [$tool] Running time: $runningtime\n");
				`touch "$outdir/.finished"` if !$code;
			}

		}
	}	
}
else{
	&_notify("[RUN_TOOL] Skipped.\n");
}

if ($tools->{system}->{RUN_TOOL_AS_JOB}){
	my $timelimit = $tools->{system}->{QSUB_TIMELIMIT};
	my $wait=1;
	my $qstat_cmd = $tools->{system}->{QSTAT_COMMAND};
	my $qdel_cmd = $tools->{system}->{QDEL_COMMAND};
	while ($wait){
		my $qsub_job_exit=0;
		foreach my $id (keys %job_ids){
			if ($job_ids{$id}->{finish}){
				$qsub_job_exit++;
				next;
			}
			my $job_status=`$qstat_cmd $id 2>&1`;
			my $time = $job_ids{$id}->{time};
			my $tool = $job_ids{$id}->{tool};
			my $outdir = $job_ids{$id}->{outdir};
			my $log = $job_ids{$id}->{log};
			if ($job_status =~ /do not exist/){
				$qsub_job_exit++;
				my $runningtime = &timeInterval($time);
				if (`grep -a -i ERROR $log`){
					&_notify("[RUN_TOOL] [$tool] Error occured.\n");
				}else{
					`touch "$outdir/.finished"`;
				}
				&_notify("[RUN_TOOL] [$tool] Running time: $runningtime\n");
				$job_ids{$id}->{finish}=1;
			}
			if ((time-$time) > $timelimit){
				`$qdel_cmd $id`;
				$qsub_job_exit++;
				$job_ids{$id}->{finish}=1;
				&_notify("[RUN_TOOL] [$tool] Error: TIME OUT\n");
			}
		}
		last if scalar(keys %job_ids) ==  $qsub_job_exit; 
		last if ! %job_ids;
		sleep 5;
	}		
}
foreach (@childs) {
	my $tmp = waitpid($_, 0);
}

#generate post-process script
my $heatmap_scale = $tools->{system}->{HEATMAP_SCALE} ? $tools->{system}->{HEATMAP_SCALE} : 'log';
my $heatmap_top   = $tools->{system}->{HEATMAP_DISPLAY_TOP} ? $tools->{system}->{HEATMAP_DISPLAY_TOP} : 0;

foreach my $idx ( sort {$a<=>$b} keys %$file_info ){
	my $fa = $file_info->{$idx}->{FASTA};
	my $fnb = $file_info->{$idx}->{PREFIX};
	my $tmpdir = "$p_outdir/temp";
	my $gottcha_present=0;

	my $pwd = `pwd`;
	chomp $pwd;

	print $post_fh "
      export PATH=$Bin:$Bin/../:$Bin/script/:$Bin/bin/:$ENV{EDGE_HOME}/thirdParty/Anaconda2/bin:\$PATH;
	  cd $pwd;

      echo \"[Post-processing #$idx $fnb]\";
      mkdir -p $tmpdir;
	";

	foreach my $tool ( sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools )
	{
		next if $tool eq 'system';
		
		print $post_fh "\n(\n";
	
		$gottcha_present=1 if $tool =~ /gottcha/i;

		my $outdir = "$p_outdir/$idx\_$fnb/$tool";
		my $tool_rep_dir = "$p_repdir/$idx\_$fnb/$tool";
		my $prefix = "$fnb";

		print $post_fh "echo \"==> processing result: $tool\";\n";

		# copy output list & krona files
		print $post_fh "
          mkdir -p $tool_rep_dir
          echo \"====> Copying result list to report directory...\";
          [[ -e \"$outdir/$prefix.out.list\" ]] && cp  $outdir/$prefix.out.list $tool_rep_dir/$fnb-$tool.list.txt;
          [[ -e \"$outdir/$prefix.out.megan\" ]] && cp  $outdir/$prefix.out.list $tool_rep_dir/$fnb-$tool.megan;
          [[ -e \"$outdir/$prefix.full.tsv\" ]] && cp  $outdir/$prefix.full.tsv $tool_rep_dir/$fnb-$tool.full.tsv;
          [[ -e \"$outdir/$prefix.krona.html\" ]] && cp  $outdir/$prefix.krona.html $tool_rep_dir/$fnb-$tool.krona.html;
          [[ -e \"$outdir/$prefix.classification.csv\" ]] && cp $outdir/$prefix.classification.csv $tool_rep_dir/$fnb-$tool.classification.csv;
          [[ -e \"$outdir/.finished\" ]] && cp  $outdir/.finished $tool_rep_dir/.finished;
          if [ -e \"$outdir/$prefix.sam\" ]
          then
            samtools view -b -@ $threads -S $outdir/$prefix.sam -o $tool_rep_dir/$fnb-$tool.bam 2>/dev/null;
          fi 
          if [ -e \"$outdir/$prefix.gottcha.sam\" ]
          then
            samtools view -b -@ $threads -S $outdir/$prefix.gottcha.sam -o $tool_rep_dir/$fnb-$tool.bam 2>/dev/null;
          fi
          if [ -e \"$outdir/$prefix.pangia.sam\" ] 
          then
            cp -rf $outdir/* $tool_rep_dir/;
	  fi
          if ls $outdir/$prefix.gottcha_*.sam 1>/dev/null 2>&1
          then 
            cp -f $outdir/$prefix.gottcha_*.sam $tool_rep_dir/;
          fi
          if [ -e \"$outdir/$prefix.out.read_classification\" ]
          then
            cp $outdir/$prefix.out.read_classification $tool_rep_dir/$fnb-$tool.read_classification
          fi

          echo \"====> Generating phylo_dot_plot for each tool...\";
          if [ -e \"$outdir/$prefix.out.tab_tree\" ]
          then
            phylo_dot_plot.pl -i $outdir/$prefix.out.tab_tree --score $outdir/$prefix.out.tab_tree.score -p $outdir/$prefix.tree -t $fnb-$tool
          fi
		";

        print $post_fh "
		  if [ -e \"$outdir/$prefix.tree.svg\" ]
		  then
              cp $outdir/$prefix.tree.svg $tool_rep_dir/$fnb-$tool.tree.svg
		  fi
		";

		#need to be done once.
		if( $idx == 1 ){
			foreach my $rank (("genus","species","strain")){
				if ( ! -e "$p_repdir/heatmap_TOOL-$tool.$rank.pdf" ){
					print $post_fh "merge_list_specTaxa_by_tool.pl $p_outdir/*/$tool/*.list -p $fnb --top $heatmap_top -l $rank -otu $tmpdir/$tool.$rank.otu.txt -output $tmpdir/$tool.$rank.heatmap.matrix;\n";
					print $post_fh "biom convert -i $tmpdir/$tool.$rank.otu.txt -o $tool_rep_dir/$fnb-$tool-$rank.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null\n";
					print $post_fh "heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s $heatmap_scale --in $tmpdir/$tool.$rank.heatmap.matrix --out $p_repdir/heatmap_TOOL-$tool.$rank.pdf  2>/dev/null\n";
				}
			}
		}
		
		print $post_fh ")&\n";
	}

	print $post_fh "\nwait\n";

	my $fnb_rep_dir = "$p_repdir/$idx\_$fnb";

	print $post_fh "
echo \"==> Generating Radar Chart...\";
convert_list2radarChart.pl --level genus   --outdir $p_repdir --outprefix radarchart_DATASET $fnb_rep_dir/*/*.list.txt &
convert_list2radarChart.pl --level species --outdir $p_repdir --outprefix radarchart_DATASET $fnb_rep_dir/*/*.list.txt &
convert_list2radarChart.pl --level strain  --outdir $p_repdir --outprefix radarchart_DATASET $fnb_rep_dir/*/*.list.txt &

echo \"==> Generating matrix for heatmap by tools...\";
merge_list_specTaxa_by_dataset.pl $fnb_rep_dir/*/*.list.txt --top $heatmap_top -l genus   --otu $tmpdir/$fnb.genus.otu.txt --output $tmpdir/$fnb.genus.heatmap.matrix & 
merge_list_specTaxa_by_dataset.pl $fnb_rep_dir/*/*.list.txt --top $heatmap_top -l species --otu $tmpdir/$fnb.species.otu.txt --output $tmpdir/$fnb.species.heatmap.matrix &
merge_list_specTaxa_by_dataset.pl $fnb_rep_dir/*/*.list.txt --top $heatmap_top -l strain  --otu $tmpdir/$fnb.strain.otu.txt --output $tmpdir/$fnb.strain.heatmap.matrix &

wait 

biom convert -i $tmpdir/$fnb.genus.otu.txt -o $p_repdir/$fnb.genus.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &
biom convert -i $tmpdir/$fnb.species.otu.txt -o $p_repdir/$fnb.species.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &
biom convert -i $tmpdir/$fnb.strain.otu.txt -o $p_repdir/$fnb.strain.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &

wait

echo \"==> Generating heatmaps...\";
heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s $heatmap_scale --in $tmpdir/$fnb.genus.heatmap.matrix   --out $p_repdir/heatmap_DATASET-$fnb.genus.pdf --title $fnb.genus 2>/dev/null &
heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s $heatmap_scale --in $tmpdir/$fnb.species.heatmap.matrix --out $p_repdir/heatmap_DATASET-$fnb.species.pdf --title $fnb.species 2>/dev/null  &
heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s $heatmap_scale --in $tmpdir/$fnb.strain.heatmap.matrix  --out $p_repdir/heatmap_DATASET-$fnb.strain.pdf --title $fnb.strain 2>/dev/null &
";
	print $post_fh "\nwait\n";
	print $post_fh "echo \"[END #$idx $fnb]\"\n\n";
}


my $hl_flag= "--highlight_list=$hl_list" if ( defined $hl_list && -e $hl_list );
print $post_fh "\necho \"==> Generating TOP$rep_top contamination report...\"
convert_list2report.pl --top $rep_top --list $filelist_out --setting $opt{setting} --output $p_outdir > $summary_out &
\necho \"==> Generating resource usage report...\"
uge_helper -l $logfile > $res_usage_out &
\necho \"==> Producing report XLSX file...\"
generate_xlsx_report.pl $hl_flag --list $filelist_out --setting $opt{setting} --output $p_outdir &
wait
";

close $log_fh  if $log_fh;
close $post_fh if $post_fh;
close $sqdx_fh if $sqdx_fh;
close $info_fh if $info_fh;
close $list_fh if $list_fh;

if( $tools->{system}->{RUN_POST_PROCESS} ){
	&_notify("\n[POST PROCESS] Generate report...\n");
	`bash $post_script`;
	&_notify("\n[POST PROCESS] Done.\n");
}

#clean up output directory
unless( $opt{debug} ){
	`rm -rf $p_outdir/script`;
	`rm -rf $p_outdir/temp`;
	`rm -rf $p_outdir/*_*`;
	`rm -rf $p_seqdir`;
#	`rm -rf $p_logdir`;
}

###############################################################

sub prepSequence {
	my $file_info = shift;
	my $tools = shift;

	#flags
	my ($FASTQ, $FASTQSE, $FASTQPE, $FASTA, $FASTA_EXTRACT, $SPLITRIM_DIR) = (0,0,0,0,0,0);
	
	foreach my $tool ( keys %$tools ){
		my $cmd = $tools->{$tool}->{COMMAND};
		$FASTQ        = 1 if $cmd =~ /%FASTQ%/;
		$FASTQSE      = 1 if $cmd =~ /%FASTQSE%/;
		$FASTQPE      = 1 if $cmd =~ /%FASTQPE%/;
		$FASTA        = 1 if $cmd =~ /%FASTA%/;
		$FASTA_EXTRACT = 1 if $cmd =~ /%FASTA_EXTRACT%/;
		$SPLITRIM_DIR = 1 if $cmd =~ /%SPLITRIM_DIR%/;
	}
	$FASTA = 1 if $FASTA_EXTRACT;

	foreach my $count ( sort {$a<=>$b} keys %$file_info ){
		&_verbose("\n[PREP_SEQ] Processing #$count sequence...\n");

		my $fnb = $file_info->{$count}->{PREFIX};
		unless( $fnb ){
			($fnb) = $file_info->{$count}->{FASTQ}   =~ /([^\/]+)\.[^\.]+$/ if defined $file_info->{$count}->{FASTQ};
			($fnb) = $file_info->{$count}->{FASTQPE} =~ /([^\/]+)\.[^\.]+$/ if defined $file_info->{$count}->{FASTQPE};
			$fnb = "Dataset$count" unless $fnb;
			$file_info->{$count}->{PREFIX} = $fnb;
			&_verbose("[PREP_SEQ] Generate filename base for output prefix: $fnb\n");
		}
		else{
			&_verbose("[PREP_SEQ] PREFIX: ".$file_info->{$count}->{PREFIX}."\n");
		}
	
		if( !-e $file_info->{$count}->{FASTQ} && $FASTQ ){
			&_verbose("[PREP_SEQ] FASTQ not found! Generate FASTQ seq: $p_seqdir/$fnb.fastq\n");
			my @pe = split /,/, $file_info->{$count}->{FASTQPE};
			if( -e $pe[0] && -e $pe[1] ){
				my $cmd = "cat $pe[0] $pe[1] > $p_seqdir/$fnb.fastq";
				`$cmd`;
				$file_info->{$count}->{FASTQ} = "$p_seqdir/$fnb.fastq";
			}
		}
		else{
			&_verbose("[PREP_SEQ] FASTQ: ".$file_info->{$count}->{FASTQ}."\n");
		}

		my $file = $file_info->{$count}->{FASTQ};
		my ($path) = $file =~ /(.*)\/[^\/]+$/;
		$path = "." unless $path;

		#splitrim
		#my ($TRIM_FIXL, $TRIM_MINQ) = (30,20);
		#if( $SPLITRIM_DIR ){
		#	if( !-s "$p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}/${fnb}_splitrim.fastq" ){
		#		&_verbose("[PREP_SEQ] SPLITRIM_DIR not found! Generate SPLITRIM_DIR: $p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}/\n");
		#		my $fastq =  $file_info->{$count}->{FASTQ};
		#		my $cmd = "$ENV{EDGE_HOME}/thirdParty/gottcha/bin/splitrim --inFile=$fastq --fixL=$TRIM_FIXL --recycle --minQ=$TRIM_MINQ --prefix=$fnb --outPath=$p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ} > /dev/null 2>&1";
		#		my $exe = system($cmd);
		#		die "Can't not generate splitrim directory.\n" if $exe;
		#		$file_info->{$count}->{SPLITRIM_DIR} = "$p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}";
		#	}	
		#	else{
		#		&_verbose("[PREP_SEQ] SPLITRIM_DIR: "."$p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}"."\n");
		#	}
		#	$file_info->{$count}->{SPLITRIM_DIR} = "$p_seqdir/splitrim_fixL${TRIM_FIXL}Q${TRIM_MINQ}";
		#}

		#convert FASTQ to FASTA
		if( !defined $file_info->{$count}->{FASTA} && $FASTA ){
			$file_info->{$count}->{FASTA} = "$path/$fnb.fa"        if -e "$path/$fnb.fa";
			$file_info->{$count}->{FASTA} = "$path/$fnb.fasta"     if -e "$path/$fnb.fasta";
			$file_info->{$count}->{FASTA} = "$p_seqdir/$fnb.fa"    if -e "$p_seqdir/$fnb.fa";
			$file_info->{$count}->{FASTA} = "$p_seqdir/$fnb.fasta" if -e "$p_seqdir/$fnb.fasta";
			
			unless( defined $file_info->{$count}->{FASTA} ){
				&_verbose("[PREP_SEQ] FASTA not found. Generate FASTA sequence: $p_seqdir/$fnb.fasta\n");
				# generate fasta
				$cmd = "fastq_to_fasta_fast < $file > $p_seqdir/$fnb.fasta";
				`$cmd`;
				$file_info->{$count}->{FASTA} = "$p_seqdir/$fnb.fasta";
			}
			else{
				&_verbose("[PREP_SEQ] FASTA found: ".$file_info->{$count}->{FASTA}."\n");
			}
		}
		else{
			&_verbose("[PREP_SEQ] FASTA: ".$file_info->{$count}->{FASTA}."\n");
		}
	
		#extract FASTA
		if( !defined $file_info->{$count}->{FASTA_EXTRACT} && $FASTA_EXTRACT ){
			$file_info->{$count}->{FASTA_EXTRACT} = "$path/$fnb.extract.fa"        if -e "$path/$fnb.extract.fa";
			$file_info->{$count}->{FASTA_EXTRACT} = "$path/$fnb.extract.fasta"     if -e "$path/$fnb.extract.fasta";
			$file_info->{$count}->{FASTA_EXTRACT} = "$p_seqdir/$fnb.extract.fa"    if -e "$p_seqdir/$fnb.extract.fa";
			$file_info->{$count}->{FASTA_EXTRACT} = "$p_seqdir/$fnb.extract.fasta" if -e "$p_seqdir/$fnb.extract.fasta";

			unless( defined $file_info->{$count}->{FASTA_EXTRACT} ){
				# generate fasta
				&_verbose("[PREP_SEQ] Extracted FASTA not found. Extracting FASTA ($extract): $p_seqdir/$fnb.extract.fasta\n");
				$cmd = "extract_random_sequences.pl -n $extract -i $file_info->{$count}->{FASTA} -o $p_seqdir/$fnb.extract.fasta";
				&_log("[PREP_SEQ] CMD=$cmd\n");
				`$cmd`;
				$file_info->{$count}->{FASTA_EXTRACT} = "$p_seqdir/$fnb.extract.fasta";
			}
			else{
				&_verbose("[PREP_SEQ] FASTA_EXTRACT found: ".$file_info->{$count}->{FASTA_EXTRACT}."\n");
			}
		}
		else{
			&_verbose("[PREP_SEQ] FASTA_EXTRACT: ".$file_info->{$count}->{FASTA_EXTRACT}."\n");
		}
	}
}

sub param_replace {
	my ($cmd, $info, $tool_ref, $idx, $tool) = @_;

	#variables defined in tool sections
	foreach my $key ( keys %{$tool_ref->{$tool}} ){
		my $val = $tool_ref->{$tool}->{$key};
		$cmd =~ s/%$key%/$val/g;
	}

	#variables defined in [system]
	my $system = $tool_ref->{system};
	foreach my $key ( keys %$system ){
		my $val = $system->{$key};
		$cmd =~ s/%$key%/$val/g;
	}

	#pre-fefined variables
	$cmd =~ s/%PREFIX%/$info->{$idx}->{PREFIX}/g;
	$cmd =~ s/%FASTQ%/$info->{$idx}->{FASTQ}/g;
	$cmd =~ s/%FASTQSE%/$info->{$idx}->{FASTQSE}/g;
	$cmd =~ s/%FASTQPE%/$info->{$idx}->{FASTQPE}/g;
	$cmd =~ s/%FASTA%/$info->{$idx}->{FASTA}/g;
	$cmd =~ s/%SPLITRIM_DIR%/$info->{$idx}->{SPLITRIM_DIR}/g;
	$cmd =~ s/%FASTA_EXTRACT%/$info->{$idx}->{FASTA_EXTRACT}/g;
	$cmd =~ s/%TOOL%/$tool/g;
	$cmd =~ s/%SERIAL%/$idx/g;

	return $cmd;
}

sub restore_settings {
	my ( $file ) = @_;
	open FILE, $file || die "Can't open settings: $file\n";
	my $set;
	my $section;
	my $count=0;
	while(<FILE>){
		chomp;
		next if /^$/;
		next if /^#/;
		next if /^;;/;
		if ( /^\[(.+)\]$/ ){ #new section
			$section = $1;
			die "Section \"$section\" existed. Please check the setting file.\n" if defined $set->{$section};
			$set->{$section}->{ORDER} = $count++;
			next;
		}
		my ($key, $val);
		($key, $val) = $_ =~ /^([^=]+)\s*=\s*(.*)$/;
		$key = uc($key);

		$set->{$section}->{$key} = $val;
	}
	close FILE;
	return $set;
}

sub _log {
    my $msg = shift;
    print $log_fh $msg if $log_fh;
}

sub _verbose {
    my $msg = shift;
    print $msg if $opt{verbose};
    &_log($msg);
}

sub _notify {
    my $msg = shift;
    print $msg;
    &_log($msg);
}

sub _notifyError {
    my $msg = shift;
    &_log($msg);
    die $msg;
}

sub countFastq_exe  
{ 
	my $file=shift; 
	my $seq_count; 
	my $total_length; 
	open my $fh, $file or die "Can't open file $file.";
	while (<$fh>) 
	{  
		my $id=$_; 
		my $seq=<$fh>; 
		chomp $seq; 
		my $q_id=<$fh>; 
		my $q_seq=<$fh>; 
		my $len = length $seq; 
		$seq_count ++; 
		$total_length +=$len; 
	} 
	close $fh; 
	return ($seq_count,$total_length); 
}

#sub getCpuUsage {
#	my $cpu = `top -bn1 | grep load | awk '{printf "%.1f", \$(NF-2)}'`;
#	return $cpu;
#}

sub getCpuUsage {
	my $pid = shift;
	$pid ||= $$;
	my $pstree = `pstree -lp $pid`;
	my @pids = $pstree =~ /\((\d+)\)/mg;
	my $cpu=0;
	foreach my $pid ( @pids ){
		my $usage = `ps -p $pid -o c | tail -n1`;
		chomp $usage;
		$cpu += $usage if $usage =~ /^\d+$/;
	}
	return $cpu/100;
}

sub timeInterval{
	my $now = shift;
	$now = time - $now;
	return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}

sub usage {
print STDERR "

USAGE: perl sample_contamination.pl -s [INI_FILES] [INPUT FASTQ] ([INPUT FASTQ2] [INPUT FASTQ3]...) (-l [FILELIST]) 

    [INPUT FASTQ]        One or more input FASTAQ file need to be provided. Wildcard is allowed. User can also
	                     provide an filelist. Check \"-l\" option for more detail.

    -l [FILELIST]        This argument is optional. The filelist can be a simple FASTQ path list or in a certain
	                     format like /users/218817/bin/sample_contamination.filelist.txt.

    -s [SETTINGS_FILE]   A setting file with certain format is required for the wrapper. An example ini file can be 
                         found in /users/218817/bin/sample_contamination.settings.ini.
 
    -c [NUM]             Number of CPU (default: 4)
    
    -d                   Debug mode. Keep output and all temporary files.

    -o [OUTPUT_DIRECTORY]  

=================================================================================================================

Example:

sample_contamination.pl -s sample_contamination.settings.ini -o example_out sequence/*.fastq  
sample_contamination.pl -s sample_contamination.settings.ini -o example_out -l sequence_list.txt
sample_contamination.pl -s settings.txt -c 4 -o example_out seq1.fastq seq2.fastq seq3.fastq

";
exit(1);
}
