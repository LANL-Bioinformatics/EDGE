#!/usr/bin/env perl
#
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-09-05
#

use strict;
use FindBin qw($RealBin);
use Cwd 'abs_path';
use lib "$RealBin/../../lib";
use CGI qw(:standard);
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use POSIX qw(strftime);
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use Tie::File;

require "edge_user_session.cgi";
require "../cluster/clusterWrapper.pl";
##sample metadata
require "../metadata_scripts/metadata_api.pl";
#END

my $cgi    = CGI->new;
my %opt    = $cgi->Vars();
my $pname  = $opt{proj};
my $username = $opt{username};
my $password = $opt{password};
my $new_proj_name = $opt{rename_project};  #getting project name input value from edge.js
my $new_proj_desc = $opt{project_description};  #getting project description input value from edge.js
my $action = lc($opt{action});
my $shareEmail = $opt{shareEmail};
my $userType = $opt{userType}||"user";
my $protocol = $opt{protocol}||"http:";
my $sid = $opt{sid};
my $taxa_for_contig_extract = $opt{taxa};
my $cptool_for_reads_extract = $opt{cptool};
my $contig_id = $opt{contigID};
my $blast_params = $opt{"edge-contig-blast-params"} || " -num_alignments 10 -num_descriptions 10 -evalue 1e-10 " ;
my $domain	= $ENV{'HTTP_HOST'}|| 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";
$ENV{PATH} = "$EDGE_HOME/bin:$ENV{PATH}";

$pname      ||= $ARGV[0];
$action     ||= $ARGV[1];
$username   ||= $ARGV[2];
$password   ||= $ARGV[3];
$shareEmail ||= $ARGV[4];
$sid        ||= $ARGV[5];
$domain     ||= $ARGV[6];
my $umSystemStatus = $ARGV[7];

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
$sys->{edgeui_input} = "$sys->{edgeui_input}"."/$webhostname" if ( -d "$sys->{edgeui_input}/$webhostname");
my $out_dir     = $sys->{edgeui_output};
my $input_dir   = $sys->{edgeui_input};
my $www_root	= $sys->{edgeui_wwwroot};
my $um_url      = $sys->{edge_user_management_url};
my $keep_days	= $sys->{edgeui_proj_store_days};
$domain       ||= "edgeset.lanl.gov";
$um_url	      ||= "$protocol//$domain/userManagement";
$out_dir      ||= "/tmp"; #for security
$umSystemStatus ||= $sys->{user_management} if (! @ARGV);
my $info;
my $proj_dir    = abs_path("$out_dir/$pname");
my $list;
my $permission;

#cluster
my $cluster 	= $sys->{cluster};
my $cluster_job_prefix = $sys->{cluster_job_prefix};
my $cluster_qsub_options= $sys->{cluster_qsub_options};
&LoadSGEenv($sys) if ($cluster);

#check projects vital
my ($vital, $name2pid, $error);
if($cluster) {
	($vital, $name2pid, $error) = checkProjVital_cluster($cluster_job_prefix);
	if($error) {
		$info->{INFO} = "ERROR: $error";
	}
} else {
	($vital, $name2pid) = &checkProjVital();
}


my $time = strftime "%F %X", localtime;

my ($memUsage, $cpuUsage, $diskUsage) = &getSystemUsage();
$info->{STATUS} = "FAILURE";
#$info->{INFO}   = "Project $pname not found.";
if ( ($memUsage > 99 or $cpuUsage > 99) and $action ne 'interrupt' and !$cluster){
        $info->{INFO}   =  "No enough CPU/MEM resource to perform action. Please wait or contact system administrator.";
        &returnStatus();
}

#session check
my $real_name = $pname;
my $projCode;
my $projStatus;
my @projCodes = split /,/,$opt{proj} if ($action eq 'compare');
my $user_proj_dir = "$input_dir/tmp";
if ( $umSystemStatus )
{
	my $valid = verifySession($sid);
	unless($valid){
		$info->{INFO} = "ERROR: Invalid session found.";
		&returnStatus() if (!@ARGV);
	}
	else{
		($username,$password) = getCredentialsFromSession($sid);
	}
	
	$list = &getUserProjFromDB("owner");

	($real_name,$projCode,$projStatus)= &getProjNameFromDB($pname) if ($action ne 'compare');
	
	$user_proj_dir = "$input_dir/". md5_hex($username)."/MyProjects/$real_name"."_".$pname;
	#separate permission for future uses. A permission module can be added potentially..
	if( defined $list->{$pname} || $userType =~ /admin/){
		$permission->{empty} = 1;
		$permission->{remove} = 1;
		$permission->{delete} = 1;
		$permission->{interrupt} = 1;
		$permission->{rerun} = 1;
		$permission->{archive} = 1;
		$permission->{share} = 1;
		$permission->{unshare} = 1;
		$permission->{publish} = 1;
		$permission->{unpublish} = 1;
		$permission->{tarproj} = 1;
		$permission->{getcontigbytaxa} = 1;
		$permission->{getreadsbytaxa} = 1;
		$permission->{metadata} = 1;
	}
	#print STDERR "User: $username; Sid: $sid; Valid: $valid; Pname: $pname; Realname: $real_name; List:",Dumper($list),"\n";
}else{
	($real_name,$projCode,$projStatus)= &scanProjToList($out_dir,$pname) if ($action ne 'compare');
	if (!$real_name){
		$info->{INFO} = "ERROR: No project with ID $pname.";
		&returnStatus();
	}
}
	$proj_dir = abs_path("$out_dir/$projCode") if ( -d "$out_dir/$projCode");

if ($action eq 'rename' ){
	renameProject($new_proj_name,$new_proj_desc);
	#edgeDB: update run name
	my $runFile = "$proj_dir/metadata_run.txt";
	if(-e $runFile) {
		my $runId = `grep -a "edge-run-id=" $runFile | awk -F'=' '{print \$2}'`;
		chomp $runId;
		`perl edge_db.cgi run-update "$runId" "$new_proj_name"`;
	}
	#edgeDB
}

if( $action eq 'empty' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	if( $name2pid->{$pname} || $name2pid->{$projCode}){
		$info->{INFO} = "ERROR: Project $real_name is running.";
		&returnStatus();
	}

	if( -d $proj_dir ){
		opendir(BIN, $proj_dir) or die "Can't open $proj_dir: $!";
		while( defined (my $file = readdir BIN) ) {
			next if $file eq '.' or $file eq '..';
			`rm -rf $proj_dir/$file` if -d "$proj_dir/$file";
			`rm -f $proj_dir/\.run*`;
		}
		closedir(BIN);

		#prepare new preject log
		`cp $proj_dir/process.log $proj_dir/process.log.bak`;
		`echo "\n*** [$time] EDGE_UI: This project has been emptied ***\n" |tee $proj_dir/process.log > $proj_dir/process_current.log`;
		`grep -a "runPipeline -c" $proj_dir/process.log.bak >> $proj_dir/process.log`;
		`echo "*** [$time] EDGE_UI: project unstarted ***" >> $proj_dir/process.log`;

		opendir(BIN, $proj_dir) or die "Can't open $proj_dir: $!";
		while( defined (my $file = readdir BIN) ) {
			next if $file eq '.' or $file eq '..';
			if( -d "$proj_dir/$file" ){
				$info->{STATUS} = "FAILURE";
				$info->{INFO} = "ERROR: Not able to delete $proj_dir/$file.";
				&returnStatus();
			}
		}
		closedir(BIN);

		$info->{STATUS} = "SUCCESS";
		$info->{INFO} = "Project output has been emptied.";
	}
	else{
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Output directory not found.";
	}
}
elsif( $action eq 'remove' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Cannot remove $real_name from project list.";

	if( -s "$proj_dir/process.log" ){
		`mv $proj_dir/process.log $proj_dir/process.log.bak`;
		if( !-e "$proj_dir/process.log" ){
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "$real_name has been removed from project list.";
		}
	}
	&updateDBProjectStatus($pname,"not list") if ($username && $password);
}
elsif( $action eq 'delete' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	if( -d $proj_dir ){
		#update project list
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Failed to delete the output directory.";

		my $pid = $name2pid->{$pname} || $name2pid->{$projCode};
		if( $pid ){
			my $invalid;
			if($cluster) {
				$invalid = clusterDeleteJob($pid);
				if( $invalid ){
					$info->{INFO} = "Failed to kill the running cluster job. (Job ID: $pid)";
				}
			} else {
				$invalid = &killProcess($pid);
				if( $invalid ){
					$info->{INFO} = "Failed to kill the running process. (PID: $pid)";
				}
			}
		}
	
		#edgeDB: delete run name
		my $runFile = "$proj_dir/metadata_run.txt";
		if(-e $runFile) {
			my $runId = `grep -a "edge-run-id=" $runFile | awk -F'=' '{print \$2}'`;
			chomp $runId;
			`perl edge_db.cgi run-delete "$runId" "$new_proj_name"`;
		}
		#edgeDB
	
		if ($username && $password){
			&updateDBProjectStatus($pname,"delete");
			`rm -f $user_proj_dir`;
			`rm -f $input_dir/public/projects/${real_name}_$pname`;
			`rm -f $input_dir/*/SharedProjects/${real_name}_$pname`;
		}
		`rm -rf $proj_dir`;
		`rm -rf $out_dir/$pname`;
		`rm -f $input_dir/../JBrowse/data/$pname $input_dir/../JBrowse/data/$projCode`;
		if( !-e $proj_dir && !-e "$out_dir/$pname" ){
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "Project $real_name has been deleted.";
		}

	}
	else{
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Output directory not found.";
	}
}
elsif( $action eq 'interrupt' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to stop EDGE process.";
	
	my $pid = $name2pid->{$pname} || $name2pid->{$projCode} || $projStatus eq "unstarted"; 

	if( $pid ){
		my $invalid;
		if($cluster) {
			$invalid = clusterDeleteJob($pid);
			if( !$invalid ){
				`echo "\n*** [$time] EDGE_UI: This project has been interrupted. ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
				$info->{STATUS} = "SUCCESS";
				$info->{INFO}   = "The cluster job (JOB ID: $pid) has been stopped.";
			}
		} else {
			$invalid = &killProcess($pid);
			if( !$invalid  || $projStatus eq "unstarted"){
				`echo "\n*** [$time] EDGE_UI: This project has been interrupted. ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
				$info->{STATUS} = "SUCCESS";
				$info->{INFO}   = "The process (PID: $pid) has been stopped.";
			}
		}
	}
	else{
		$info->{INFO} = "Project $real_name is not running.";
	}
}
elsif( $action eq 'rerun' ){
	if( $umSystemStatus && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to rerun project $real_name.";

	if( -e "$proj_dir/config.txt.bak"){
		cleanProjectForNewConfig();
	}

	my $pid = $name2pid->{$pname} || $name2pid->{$projCode};

	if( ! defined $pid ){
		if($cluster) {
			my $cluster_job_script = "$proj_dir/clusterSubmit.sh";
			if(!-e $cluster_job_script) {
				$info->{INFO} = "Failed to restart this project. File $cluster_job_script not found.";
			} else {
				&updateDBProjectStatus($pname,"running") if ($username && $password);
				my ($job_id,$error) = clusterSubmitJob($cluster_job_script,$cluster_qsub_options);
				if($error) {
					$info->{INFO} = "Failed to restart this project: $error";
				} else {
					$info->{STATUS} = "SUCCESS";
					$info->{INFO}   = "Project $real_name has been restarted (JOB ID: $job_id).";
					$info->{PID}    = $job_id;
					`echo "\n*** [$time] EDGE_UI: This project has been restarted. (unstarted) ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
				}
			}
		} else {
			my $cmd = "";
			my $process_log=(-e "$proj_dir/process.log")? "$proj_dir/process.log":"$proj_dir/process.log.bak";
			open LOG, "$process_log" or die "Can't open process log:$!.";
			foreach(<LOG>){
				chomp;
				if( /runPipeline -c / ){
					$cmd = $_;
				}
			}
			close LOG;
			my ($numcpu) = $cmd =~ /-cpu (\d+)/;
			my $run = &availableToRun($numcpu);
			if (!$run){
				 my $time = strftime "%F %X", localtime;
				`echo "\n*** [$time] EDGE_UI: This project is queued. ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
				`echo "$cmd" >> $proj_dir/process.log`;
				`echo "\n*** [$time] EDGE_UI: Project unstarted ***" >> $proj_dir/process.log`;
				$info->{INFO} = "The server does not have enough CPU available to run this job. The job is queued";
				&updateDBProjectStatus($pname,"unstarted") if ($username && $password);
				&returnStatus();
			}
			if( $cmd ){
				chdir($proj_dir);
				#remove cached report/status
				#`rm -f $proj_dir/.run.complete.status.json`;
				`rm -f $proj_dir/HTML_Report/.complete_report_web`;
				my $newpid = open RUNPIPLINE, "-|", "$cmd > $proj_dir/process_current.log 2>&1 &" or die $!;
				close RUNPIPLINE;
				if( $newpid ){
					$newpid++;
					$info->{STATUS} = "SUCCESS";
					$info->{INFO}   = "Project $real_name has been restarted (PID: $newpid).";
					$info->{PID}    = $newpid;
					&updateDBProjectStatus($pname,"running") if ($username && $password);
				}
				else{
					$info->{INFO} = "Failed to restart this project.";
				}
			}
			else{
				$info->{INFO} = "Failed to restart this project. No runPipeline command found.";
			}
		}
	}
	else{
		$info->{INFO} = "Project $real_name can't be restarted because it's still running (PID: $pid).";
	}
}
elsif( $action eq 'archive' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	my $adir = $sys->{edgeui_archive};

	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to move project $real_name to $adir";
	
	if( ! defined $adir ){
		$info->{INFO} = "No archive directory configured. Please check EDGE_UI settings.";
	}
	elsif( ! -w $adir ){
		$info->{INFO} = "Failed to archive this project. Archive directory is not writable by the user of webserver.";
	}
	elsif( -d "$adir/$pname" ){
		$info->{INFO} = "Archive directory ($adir/$pname) existed. Action aborted.";
	}
	else{
		my $cmd = "$RealBin/edge_archive.pl $proj_dir $adir/$pname &";
		my $pid = open ARCHIVE, "-|", $cmd or die $!;
		close ARCHIVE;

		$pid++;
		if( $pid ){
			open LOG, ">>$proj_dir/process.log" or die "Can't open process log:$!.";
			print LOG "\n[Archive Project]\nDoArchive=1\n[Archive Project]\n Running\n $cmd\n";
			close LOG;
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "Start archiving project $real_name.";
		}
		if ($username && $password){
			&updateDBProjectStatus($pname,"archived");
			`rm -f $user_proj_dir`;
			`rm -f $input_dir/public/projects/${real_name}_$pname`;
			`rm -f $input_dir/*/SharedProjects/${real_name}_$pname`;
		}
	}
}
elsif( $action eq 'tarproj'){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	my $tarFile = "$proj_dir/$real_name.tgz";
	my $tarDir =  "$proj_dir";
	if ($username && $password){
		$tarFile = "$proj_dir/${real_name}_$pname.tgz";
		$tarDir = "$out_dir/${real_name}_$pname";
	}
	$tarDir =~ s/$out_dir//;
	$tarDir =~ s/^\///;
	(my $tarLink =  $tarFile ) =~ s/$www_root//;
	$tarLink =~ s/^\///;
	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to tar project $real_name $tarLink";
	chdir $out_dir;
	if ( ! -e "$proj_dir/.tarfinished" || ! -e $tarLink){
		`ln -s $proj_dir $tarDir` if ($username && $password);
		my $cmd = "tar --exclude=\"*gz\" --exclude=\"*.sam\" --exclude=\"*.bam\" --exclude=\"*.fastq\" -cvzf $tarFile  $tarDir/* ;  touch $proj_dir/.tarfinished ";	
		my $pid;
		if (@ARGV){
			$pid=`$cmd`;
		}else{
			$pid = open TARPROJ, "-|", $cmd or die $!;
			close TARPROJ;
			$pid++;
		}
		if( $pid ){
			unlink $tarDir if ($username && $password);
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "$real_name compressed file is ready to ";
			$info->{LINK}	= "<a data-ajax='false' id='ddownload_link'  href=\"$tarLink\">download</a>";
		}else{
			$info->{INFO}   = "Project $real_name tar file existed";
		}
	}
}
elsif( $action eq 'getcontigbytaxa'){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}	
	my $assemble_outdir="$proj_dir/AssemblyBasedAnalysis";
	my $taxa_outdir="$assemble_outdir/Taxonomy";
	(my $relative_taxa_outdir=$taxa_outdir) =~ s/$www_root//;
	(my $out_fasta_name = $taxa_for_contig_extract) =~ s/[ .']/_/;
	$out_fasta_name = "$real_name"."_"."$out_fasta_name.fasta";
	my $cmd = "$EDGE_HOME/scripts/contig_classifier_by_bwa/extract_fasta_by_taxa.pl -fasta $assemble_outdir/${real_name}_contigs.fa -csv $taxa_outdir/$real_name.ctg_class.top.csv -taxa \"$taxa_for_contig_extract\" -rank genus > $taxa_outdir/$out_fasta_name";
	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to extract $taxa_for_contig_extract contig fasta";
	
	if (  -s  "$taxa_outdir/$out_fasta_name.fasta"){
		$info->{STATUS} = "SUCCESS";
		$info->{PATH} = "$relative_taxa_outdir/$out_fasta_name";
	}else{
		my $pid = open EXTRACTCONTIG, "-|", $cmd or die $!;
		close EXTRACTCONTIG;
		$pid++;

		if( $pid ){
			$info->{STATUS} = "SUCCESS";
			$info->{PATH} = "$relative_taxa_outdir/$out_fasta_name";
		}
	}
}
elsif( $action eq 'getreadsbytaxa'){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	my $read_type="allReads";
	my $reads_fastq="$proj_dir/ReadsBasedAnalysis/Taxonomy/$read_type.fastq"; 
	my $readstaxa_outdir="$proj_dir/ReadsBasedAnalysis/Taxonomy/report/1_$read_type/$cptool_for_reads_extract";
     
	if ( -e "$proj_dir/ReadsBasedAnalysis/UnmappedReads/Taxonomy"){
		$read_type="UnmappedReads";
		$reads_fastq="$proj_dir/ReadsBasedAnalysis/$read_type/Taxonomy/$read_type.fastq";
		$readstaxa_outdir="$proj_dir/ReadsBasedAnalysis/$read_type/Taxonomy/report/1_$read_type/$cptool_for_reads_extract";

	}
	(my $relative_taxa_outdir=$readstaxa_outdir) =~ s/$www_root//;
	(my $out_fasta_name = $taxa_for_contig_extract) =~ s/[ .']/_/g;
	my $extract_from_original_fastq = ($cptool_for_reads_extract =~ /gottcha/i)? " -fastq $reads_fastq " : "";
	$out_fasta_name = "$real_name"."_"."$cptool_for_reads_extract"."_"."$out_fasta_name";
	my $cmd = "$EDGE_HOME/scripts/microbial_profiling/script/bam_to_fastq_by_taxa.pl -rank species  -name \"$taxa_for_contig_extract\" -prefix $readstaxa_outdir/$out_fasta_name -se -zip $extract_from_original_fastq $readstaxa_outdir/${read_type}-$cptool_for_reads_extract.bam 1>>$readstaxa_outdir/ReadsExtractLog.txt";
	#GOTTCHA2 Only
	if( $cptool_for_reads_extract =~ /gottcha2/i ){
		$cmd = "$EDGE_HOME/thirdParty/gottcha2/gottcha.py -s $readstaxa_outdir/*.sam -m extract -x $taxa_for_contig_extract -c > $readstaxa_outdir/$out_fasta_name.fastq; cd $readstaxa_outdir; zip $out_fasta_name.fastq.zip $out_fasta_name.fastq; rm $out_fasta_name.fastq";
	}
	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to extract $taxa_for_contig_extract reads fastq";
	
	if (  -s  "$readstaxa_outdir/$out_fasta_name.fastq.zip"){
		$info->{STATUS} = "SUCCESS";
		$info->{PATH} = "$relative_taxa_outdir/$out_fasta_name.fastq.zip";
	}elsif ( ! -e "$readstaxa_outdir/${read_type}-$cptool_for_reads_extract.bam" && ! glob("$readstaxa_outdir/*.sam") ){
		$info->{INFO}   = "The result bam does not exist.";
		$info->{INFO}   .= "If the project is older than $keep_days days, it has been deleted." if ($keep_days);
	}else
	{
		my $pid = open EXTRACTREADS, "-|", $cmd or die $!;
		close EXTRACTREADS;
		$pid++;

		if( $pid ){
			$info->{STATUS} = "SUCCESS";
			$info->{PATH} = "$relative_taxa_outdir/$out_fasta_name.fastq.zip";
		}
	}
}
elsif( $action eq 'share' || $action eq 'unshare' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	&shareProject($pname,$proj_dir,$shareEmail,$action);
	my $owner = $list->{$pname}->{OWNER};
	if ($action eq 'share'){
    		my $msg = "$owner has shared EDGE project $real_name to you. You can login to $protocol//$domain/edge_ui/ and see the project. Or click link below.\n\n $protocol//$domain/edge_ui/?proj=$projCode\n";
		my $subject = "EDGE project $real_name";
		&sendMail($username,$shareEmail,$subject,$msg);
	}
}
elsif( $action eq 'publish' || $action eq 'unpublish'){
	#print STDERR "USERMANAGMENT: $sys->{user_management}; $action: $permission->{$action}";
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	&publishProject($pname,$action);
	my $public_proj_dir = "$input_dir/public/projects/${real_name}_$pname";
	`ln -sf $proj_dir $public_proj_dir` if ($action eq 'publish' && ! -e "$public_proj_dir");
	`rm -f $public_proj_dir` if ($action eq 'unpublish');
}
elsif( $action eq 'disable-project-display' || $action eq 'enable-project-display'){
	&updateDBProjectDisplay($pname,$action);
}
elsif( $action eq 'compare'){
	my $compare_out_dir = "$out_dir/ProjectComparison/". md5_hex(join ('',@projCodes));
	my $projects = join(",",map { "$out_dir/$_" } @projCodes);
	(my $relative_outdir=$compare_out_dir) =~ s/$www_root//;
	$info->{PATH} = "$relative_outdir/compare_project.html";
	$info->{INFO} = "The comparison result is available <a target='_blank' href=\'$protocol//$domain/edge_ui/$relative_outdir/compare_project.html\'>here</a>";
	if ( -s "$compare_out_dir/compare_project.html"){
		$info->{STATUS} = "SUCCESS";
	}else{
		my $cmd = "$EDGE_HOME/scripts/compare_projects/compare_projects.pl -out_dir $compare_out_dir -projects $projects";
		my $pid = open COMPARE, "-|", $cmd or die $!;
		close COMPARE;
		$pid++;

		if( $pid ){
			my $err = `grep "No Taxonomy Classification" $compare_out_dir/log.txt`;
			if ($err){
				$info->{INFO} = "Error: $err";
				`rm -rf $compare_out_dir`;
				&returnStatus();
			}else{
				$info->{STATUS} = "SUCCESS";
			}
		}
	}
}elsif($action eq 'contigblast'){
	my $blast_out_dir="$proj_dir/AssemblyBasedAnalysis/ContigBlast";
	my $contig_file="$proj_dir/AssemblyBasedAnalysis/${real_name}_contigs.fa"; 
	my $nt_db="$EDGE_HOME/database/nt/nt"; 
	my $cpu = `grep -a "cpu=" $proj_dir/config.txt | awk -F"=" '{print \$2}'`;
	chomp $cpu;
	$blast_params =~ s/-num_threads\s+\d+//;
	`mkdir -p $blast_out_dir`;
	(my $relative_outdir=$blast_out_dir) =~ s/$www_root//;
	$info->{PATH} = "$relative_outdir/$contig_id.blastNT.html";
	$info->{INFO} = "The comparison result is available <a target='_blank' href=\'$protocol//$domain/edge_ui/$relative_outdir/$contig_id.blastNT.html\'>here</a>";
	if ( -s "$blast_out_dir/$contig_id.blastNT.html"){
		$info->{STATUS} = "SUCCESS";
	}else{
		my $cmd = "$EDGE_HOME/scripts/get_seqs.pl $contig_id $contig_file | blastn -query - -db $nt_db $blast_params -out $blast_out_dir/$contig_id.blastNT.html -num_threads $cpu -html ";

		my $pid = open BLAST, "-|", $cmd or die $!;
		close BLAST;

		if( $pid ){
			#parent
			$pid++;
			$info->{STATUS} = "SUCCESS";
		}else{
			#child
			close STDOUT;
		}
	}
}elsif( $action eq 'metadata-delete' ){
##sample metatdata
	$info->{STATUS} = "SUCCESS";
	$info->{INFO}   = "Project $real_name sample metadata has been deleted.";

	my $metadata = "$proj_dir/metadata_sample.txt";
	my $traveldata = "$proj_dir/metadata_travels.txt";
	my $symptomdata = "$proj_dir/metadata_symptoms.txt";
	
	if( -e $metadata ){
		if(-w $metadata) {
			`rm -f $metadata`;
			`rm -f $traveldata`;
			`rm -f $symptomdata`;
		} else {
			$info->{STATUS} = "FAILURE";
			$info->{INFO}   = "Failed to delete the sample metadata. Permission denied.";
		}
	}
	#delete HTML_Report/.complete_report_web
	`rm -f $proj_dir/HTML_Report/.complete_report_web`;
	
} elsif($action eq 'metadata-bsveadd') { 
	if( $sys->{user_management} && !$permission->{metadata} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	if(pushSampleMetadata("add", $proj_dir, $sys)) {
		$info->{STATUS} = "SUCCESS";
		$info->{INFO}   = "Project $real_name sample metadata has been submitted to the BSVE server.";
	} else {
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Failed to submit the sample metadata to the BSVE server";
	}

} elsif($action eq 'metadata-bsveupdate') {
	if( $sys->{user_management} && !$permission->{metadata} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	if(pushSampleMetadata("update", $proj_dir, $sys)) {
		$info->{STATUS} = "SUCCESS";
		$info->{INFO}   = "Project $real_name sample metadata has been updated in the BSVE server.";
	} else {
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Failed to update the sample metadata in the BSVE server";
	}

} 
#END sample metadata
elsif($action eq 'define-gap-depth'){
	my $gap_depth_cutoff =  ($opt{"gap-depth-cutoff"})? $opt{"gap-depth-cutoff"}:0;
	my $gap_out_dir="$proj_dir/ReferenceBasedAnalysis/readsMappingToRef";
	my $gap_outfile="$gap_out_dir/readsToRef_d$gap_depth_cutoff.gaps";
	my $gff_file="$proj_dir/Reference/reference.gff";
	my $gap_analysisOutfile="$gap_out_dir/Gap_d${gap_depth_cutoff}VSReference.report.txt";
	my $gap_analysisOutfile_json="$gap_out_dir/Gap_d${gap_depth_cutoff}VSReference.report.json";
	(my $relative_gap_out_dir=$gap_out_dir) =~ s/$www_root//;

	if ( $gap_depth_cutoff == 0 ){
		$info->{STATUS} = "SUCCESS";
		$info->{PATH} = "$relative_gap_out_dir/GapVSReference.report.json";
		&returnStatus();
	}
	if ( -s $gap_analysisOutfile_json){
		$info->{STATUS} = "SUCCESS";
		$info->{PATH} = "$relative_gap_out_dir/Gap_d${gap_depth_cutoff}VSReference.report.json";
		&returnStatus();
	}
	my $cmd;
	unless ( -s $gap_outfile){
		opendir( my $dh, $gap_out_dir);
		my @coverage_file =  grep { /coverage$/ && -f "$gap_out_dir/$_" } readdir($dh);
		closedir $dh;
		foreach my $file(@coverage_file){
			$cmd .= "$EDGE_HOME/scripts/gap_count.pl $gap_out_dir/$file $gap_depth_cutoff >> $gap_outfile ;";
		}
	}
	$cmd .= "$EDGE_HOME/scripts/gap_analysis.pl -gff $gff_file -gap $gap_outfile > $gap_analysisOutfile;";
	$cmd .= "$EDGE_HOME/scripts/tab2Json_for_dataTable.pl -project_dir $proj_dir -limit 0  -mode ref_gap $gap_analysisOutfile > $gap_analysisOutfile_json";
	my $pid = open GAP, "-|",$cmd or die $!;
	close GAP;
	if( $pid ){
		$pid++;
		$info->{STATUS} = "SUCCESS";
		$info->{PATH} = "$relative_gap_out_dir/Gap_d${gap_depth_cutoff}VSReference.report.json";
	}
}
&returnStatus();

######################################################

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			while(<CONF>){
				chomp;
				last if /^\[/;
				if ( /^([^=]+)=([^=]+)/ ){
					$sys->{$1}=$2;
				}
			}
		}
		last;
	}
	close CONF;
	return $sys;
}

sub killProcess {
	my $pid = shift;
	return "invalid pid" unless $pid > 1;
	my $cmd = 'killtree() { local _pid=$1; local _sig=${2:--TERM}; kill -stop ${_pid}; for _child in $(ps -o pid --no-headers --ppid ${_pid}); do killtree ${_child} ${_sig}; done; kill -${_sig} ${_pid}; }';
	`$cmd; killtree $pid 9`;
}

sub checkProjVital {
	my $ps = `ps aux | grep run[P]`;
	my $vital;
	my $name2pid;
	my @line = split "\n", $ps;
	foreach my $line ( @line ){
		my ( $pid, $proj, $numcpu) = $line =~ /^\S+\s+(\d+) .*\/(\S+) -cpu (\d+)/;
		$vital->{$pid}->{PROJ} = $proj;
		$vital->{$pid}->{CPU} = $numcpu;
		$name2pid->{$proj} = $pid;
	}
	return ($vital,$name2pid);
}

sub availableToRun {
        my $num_cpu = shift;
        my $cpu_been_used = 0;
        if( $sys->{edgeui_auto_queue} && $sys->{edgeui_tol_cpu} ){
                foreach my $pid ( keys %$vital ){
                        $cpu_been_used += $vital->{$pid}->{CPU};
                        return 0 if (($cpu_been_used + $num_cpu) > $sys->{edgeui_tol_cpu});
                }
                return 0 if ($num_cpu > $sys->{edgeui_tol_cpu});
        }
        return 1;
}


sub returnStatus {
	my $json = "{}";
	$json = to_json($info) if $info;
	$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
	print $cgi->header('application/json') unless $ARGV[0];
	print $json;
	exit;
}

sub getProjNameFromDB{
	my $project=shift;
        my %data = (
                email => $username,
                password => $password,
		project_id => $project 
        );
        # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/project/getInfo";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
        #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
	#print $result_json if (@ARGV);
        my $result =  from_json($result_json);
	if ($result->{error_msg})
	{
		 $info->{INFO} .= $result->{error_msg}."\n";;
	}
	else{
		return ($result->{name} , $result->{code}, $result->{status});
	}
}

sub renameProject{
	my $project_name = shift;
	my $project_description = shift;
	$project_description =~ s/(['"])/\\$1/g;
	my $pnameID = $pname;
	#adjust txt files. (config.txt and process.log )
	my $config_file = $proj_dir."/config.txt";
	my $process_log = $proj_dir."/process.log";
	my $config_json = $proj_dir."/config.json";
	#tie my @array, 'Tie::File', $config_file or die;
	if ($umSystemStatus){
		my %data = (
			email => $username,
			password => $password,
			project_id => $pnameID,
			new_project_name => $project_name, 
			new_description => $project_description,
		);

		#Encode the data structure to JSON
		##interacts with the java api to access the sql DB tables
		my $data = to_json(\%data);
		#w Set the request parameters
		my $url = $um_url ."WS/project/update";
		my $browser = LWP::UserAgent->new;
		my $req = PUT $url;
		$req->header('Content-Type' => 'application/json');
		$req->header('Accept' => 'application/json');
		#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
		$req->header( "Content-Length" => length($data) );
		$req->content($data);
		my $response = $browser->request($req);
		my $result_json = $response->decoded_content;
		my $result =  from_json($result_json);
		if ($result->{status})
		{
			$info->{STATUS} = "SUCCESS";
			$info->{INFO} .= " Project has been ${action}d to $project_name.";
		}
	}else{
		if ( -d "$out_dir/$project_name"){
			$info->{STATUS} = "FAILURE";
			$info->{INFO} .= " Project name is used.";
		}else{
			`mv $proj_dir $out_dir/$project_name`;
			$info->{STATUS} = "SUCCESS";
			$info->{INFO} .= " Project has been ${action}d to $project_name.";
			$proj_dir = "$out_dir/$project_name";
		}
	}
	if ($info->{STATUS} eq "SUCCESS"){
     		system("sed -i.bak 's/projname=[[:graph:]]*/projname=$project_name/; s/projdesc=[[:graph:]]*/projdesc=$project_description/;' $config_file") if ( -e $config_file); 
     		system("sed -i.bak 's/$real_name/$project_name/; s/projdesc=[[:graph:]]*/projdesc=$project_description/;' $process_log") if ( -e $process_log);
		system("sed -i.bak 's/edge-proj-name\" : \"[[:graph:]]*\"/edge-proj-name\" : \"$project_name\"/; s/edge-proj-desc\" : \"[[:graph:]]*\"/edge-proj-desc\" : \"$project_description\"/;' $config_json") if ( -e $config_json);
		unlink "$proj_dir/HTML_Report/.complete_report_web";
	}
}		

sub updateDBProjectStatus{
        my $project = shift;
        my $status = shift;
        my %data = (
                email => $username,
                password => $password,
                project_id => $project,
                new_project_status => $status
        );
        # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/project/update";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
        #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
        my $result =  from_json($result_json);
        if (! $result->{status})
        {
                $info->{INFO} .= " Update Project status in database failed.";
        }
}

sub updateDBProjectDisplay{
        my $project = shift;
        my $action = shift;
	my $display = "yes";
	if($action eq "disable-project-display") {
		$display = "no";
	}
        my %data = (
                email => $username,
                password => $password,
                project_id => $project,
                project_display => $display,
        );
        # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/user/updateProjectDisplay";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
        #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
        my $result =  from_json($result_json);
        if (! $result->{status})
        {
                $info->{INFO} .= " Update Project display in database failed.";
        } else {
       		$info->{STATUS} = "SUCCESS";
                $info->{INFO} .= " Project display has been updated.";
	}
}

sub sendMail{
  my $sender=shift;
  my $recipients=shift;
  my $subject=shift;
  my $msg=shift;
  $recipients =~ s/ //g;
  $recipients = join(',', grep (!/$sender/, split(',',$recipients)));
  if (`which sendmail`){
    open(MAIL, "|sendmail -t") or die "$!\n";
    print MAIL "To: $recipients\n";
    print MAIL "From: $sender\n";
    print MAIL "Subject: $subject\n\n";
   # print MAIL "Content-Type: text/html; charset=ISO-8859-1\n";
   # print MAIL "Content-Disposition: inline\n";
    print MAIL "$msg";
    close MAIL;
  }
}

sub shareProject{
	my $project=shift;
	my $proj_dir=shift;
	my $email=shift;
	my $action =shift;
	$email =~ s/ //g;
	# avoid share to owner self.
	$email = join(',', grep (!/$username/, split(',',$email)));
	my %data = (
                email => $username,
                password => $password,
                project_id => $project,
        );
	$data{"${action}_to"} = $email;
	 # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/project/${action}_user";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
        #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
        my $result =  from_json($result_json);
        if ($result->{status})
        {
                $info->{STATUS} = "SUCCESS";
                $info->{INFO} .= " Project $real_name has ${action}d to $email.";
		foreach (split(',',$email)){
			my $user_dir =  "$input_dir/". md5_hex($_);
			my $shared_proj_dir = "$user_dir/SharedProjects/${real_name}_$project";
			if ( $action eq "share"){
				`mkdir -p $user_dir/SharedProjects`;
				`ln -sf $proj_dir $shared_proj_dir` if (!-e $shared_proj_dir);
			}else{# unshare
				`rm -f $shared_proj_dir` if ( -e $shared_proj_dir);
			}
		}
        }else{
		$action =~ s/^([a-z])/\u$1/;
                $info->{STATUS} = "FAILURE";
                $info->{INFO} .= " $action Project $real_name to $email failed. $result->{error_msg}";
	}
}

sub publishProject{
	my $project=shift;
	my $action =shift;
	my %data = (
                email => $username,
                password => $password,
                project_id => $project,
        );
	 # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/project/$action";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
        #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
        my $result =  from_json($result_json);
        if ($result->{status})
        {
                $info->{STATUS} = "SUCCESS";
                $info->{INFO} .= " Project $real_name has been ${action}ed.";
        }else{
		$action =~ s/^([a-z])/\u$1/;
                $info->{STATUS} = "FAILURE";
                $info->{INFO} .= " $action Project $real_name failed. $result->{error_msg}";
	}
}

sub sessionCheck{
	return unless $sys->{user_management}; #user management is off
	my $valid = verifySession($sid);
	unless($valid){
	}
}


sub getUserProjFromDB{
	my $project_type = shift;
	my $viewType = "user";
	my $list;
    my %data = (
        email => $username,
        password => $password
    );
    # Encode the data structure to JSON
    #w Set the request parameters
	my $service;
	if ($username && $password){ 
		$service= ($viewType =~ /admin/i)? "WS/user/admin/getProjects" :"WS/user/getProjects";
		$data{project_type} = $project_type if ($viewType =~ /user/i);
	}else{
		$service="WS/user/publishedProjects";
	}
    my $data = to_json(\%data);
    my $url = $um_url .$service;
    my $browser = LWP::UserAgent->new;
    my $req = PUT $url;
    $req->header('Content-Type' => 'application/json');
    $req->header('Accept' => 'application/json');
    #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
    $req->header( "Content-Length" => length($data) );
    $req->content($data);

    my $response = $browser->request($req);
    my $result_json = $response->decoded_content;
	
	if ($result_json =~ /\"error_msg\":"(.*)"/)
    {
            $list->{INFO}->{ERROR}=$1;
            return;
    }
    my $array_ref =  from_json($result_json);
	#print Dumper($array_ref) if @ARGV;
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $project_name = $hash_ref->{name};
		my $projCode = $hash_ref->{code};
		my $status = $hash_ref->{status};
		next if ($status =~ /delete/i);
		next if (! -r "$out_dir/$id/process.log" && ! -r "$out_dir/$projCode/process.log");
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{REAL_PROJNAME} = $project_name;
		$list->{$id}->{PROJCODE} = $projCode;
		$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}"; 
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};
	}
	return $list;
}

sub scanProjToList{
	my $out_dir = shift;
	my $pname = shift;
	my $config_file;
	my $processLog;
	my ($projid,$projCode,$projName,$projStatus);
	if ($pname && -d "$out_dir/$pname"){
		$config_file = "$out_dir/$pname/config.txt";
		$processLog = "$out_dir/$pname/process_current.log";
	}else{
		$config_file = `grep -a "projid=$pname" $out_dir/*/config.txt | awk -F':' '{print \$1}'`;
	}
	chomp $config_file;
	return ($projName,$projCode,$projStatus) if ( ! -e $config_file);
	if ( -r "$processLog"){
		open (my $fh, $processLog);
		while(<$fh>){
			if (/queued/){
				$projStatus="unstarted";
				last;
			}
			if (/^All Done/){
				$projStatus="finished";
			}
		}
		close $fh;
	}
	open (my $fh, $config_file) or die "Cannot read $config_file\n";
	while(<$fh>){
		last if (/^\[Down/);
		$projid=$1 if (/^projid=(\S+)/);
		$projCode=$1 if (/^projcode=(\S+)/);
		$projName=$1 if (/^projname=(\S+)/);

	}
	close $fh;
	return ($projName,$projid,$projStatus);
}

sub getSystemUsage {
	my $mem = `vmstat -s | awk  '\$0 ~/total memory/ {total=\$1 } \$0 ~/free memory/ {free=\$1} \$0 ~/buffer memory/ {buffer=\$1} \$0 ~/cache/ {cache=\$1} END{print (total-free-buffer-cache)/total*100}'`;
        my $cpu = `top -bn1 | grep load | awk '{printf "%.1f", \$(NF-2)}'`;
        my $disk = `df -h $out_dir | tail -1 | awk '{print \$5}'`;
        $disk= `df -h $out_dir | tail -1 | awk '{print \$4}'` if ($disk !~ /\%/);
        $cpu = $cpu/$sys->{edgeui_tol_cpu}*100;
        $disk =~ s/\%//;
        if( $mem || $cpu || $disk ){
                $mem = sprintf "%.1f", $mem;
                $cpu = sprintf "%.1f", $cpu;
                $disk = sprintf "%.1f", $disk;
                return ($mem,$cpu,$disk);
        }
        else{
                return (0,0,0);
        }
}

sub cleanProjectForNewConfig {
	my $module_ctl;
	$module_ctl->{"Qiime analysis"}                 ->{"general"}       = "$proj_dir/QiimeAnalysis/runQiimeAnalysis.finished";
        $module_ctl->{"Download SRA"}                   ->{"general"}       = "$proj_dir/SRA_Download/DownloadSRA.finished";
        $module_ctl->{"ProPhage Detection"}             ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/Prophage/phageFinder.finished";
	$module_ctl->{"Count Fastq"}                    ->{"general"}       = "$proj_dir/QcReads/countFastq.finished";
	$module_ctl->{"Quality Trim and Filter"}        ->{"general"}       = "$proj_dir/QcReads/runQC.finished";
	$module_ctl->{"Host Removal"}                   ->{"general"}       = "$proj_dir/HostRemoval/*/run*.finished"; #  system("rm -f $outputDir/run${prefix}Removal.finished");
	$module_ctl->{"Assembly"}                       ->{"Provided"}      = "$proj_dir/AssemblyBasedAnalysis/processProvideContigs.finished"; 
	$module_ctl->{"Assembly"}                       ->{"SPAdes"}        = "$proj_dir/AssemblyBasedAnalysis/runSPAdesAssembly.finished"; #system("rm -f $outputDir/runAPAdesAssembly.finished");
	$module_ctl->{"Assembly"}                       ->{"Idba"}          = "$proj_dir/AssemblyBasedAnalysis/runIdbaAssembly.finished";
	$module_ctl->{"Assembly"}                       ->{"Idba"}          = "$proj_dir/AssemblyBasedAnalysis/runIdbaAssembly.finished";
	$module_ctl->{"Assembly"}                       ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/runAssembly.finished";
	$module_ctl->{"Reads Mapping To Contigs"}       ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/readsMappingToContig/runReadsToContig.finished";
	$module_ctl->{"Reads Mapping To Reference"}     ->{"general"}       = "$proj_dir/ReadsBasedAnalysis/readsMappingToRef/runReadsToGenome.finished";
        $module_ctl->{"Reads Mapping To Reference"}     ->{"UnmappedReads"} = "$proj_dir/ReferenceBasedAnalysis/UnmappedReads/retrieveUnmappedReads.finished";
	$module_ctl->{"Reads Taxonomy Classification"}  ->{"AllReads"}      = "$proj_dir/ReadsBasedAnalysis/Taxonomy/taxonomyAssignment.finished";
	$module_ctl->{"Reads Taxonomy Classification"}  ->{"UnmappedReads"} = "$proj_dir/ReadsBasedAnalysis/UnmappedReads/Taxonomy/taxonomyAssignment.finished";
	$module_ctl->{"Contigs Mapping To Reference"}   ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/contigMappingToRef/runContigToGenome.finished";
	$module_ctl->{"Variant Analysis"}               ->{"ReadsBased"}    = "$proj_dir/ReadsBasedAnalysis/contigMappingToRef/variantAnalysis.finished";
	$module_ctl->{"Variant Analysis"}               ->{"AssemblyBased"} = "$proj_dir/AssemblyBasedAnalysis/contigMappingToRef/variantAnalysis.finished";
	$module_ctl->{"Contigs Taxonomy Classification"}->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/Taxonomy/ContigsTaxonomy.finished";
	$module_ctl->{"Contigs Blast"}                  ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/Blast/ContigsBlast.finished";
	$module_ctl->{"Contigs Annotation"}             ->{"general"}       = "$proj_dir/AssemblyBasedAnalysis/Annotation/runAnnotation.finished";
	$module_ctl->{"Phylogenetic Analysis"}          ->{"general"}       = "$proj_dir/SNP_Phylogeny/SNPtree.finished";
	$module_ctl->{"Phylogenetic Analysis"}          ->{"SRA"}           = "$proj_dir/SNP_Phylogeny/SRAreads/download.finished";
	$module_ctl->{"Specialty Genes Profiling"}      ->{"general"}       = "$proj_dir/.runSpecialtyGenesProfiling.finished";
	$module_ctl->{"Specialty Genes Profiling"}      ->{"ReadsBased"}    = "$proj_dir/ReadsBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished";
	$module_ctl->{"Specialty Genes Profiling"}      ->{"AssemblyBased"} = "$proj_dir/AssemblyBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished";
	$module_ctl->{"Primer Validation"}              ->{"general"}       = "$proj_dir/AssayCheck/pcrDesign.finished";
	$module_ctl->{"Primer Design"}                  ->{"general"}       = "$proj_dir/AssayCheck/pcrDesign.finished";
	$module_ctl->{"Generate JBrowse Tracks"}        ->{"general"}       = "$proj_dir/JBrowse/writeJBrowseInfo.finished";
	$module_ctl->{"HTML Report"}                    ->{"general"}       = "$proj_dir/HTML_Report/writeHTMLReport.finished";

	my $new_config = &getSysParamFromConfig( "$proj_dir/config.txt" );
	my $old_config = &getSysParamFromConfig( "$proj_dir/config.txt.bak" );

	foreach my $module ( keys %$new_config ){
		foreach my $param ( keys %{$new_config->{$module}} ){
			#if one of the parameter in the module changed, reset the whole module
			if( $new_config->{$module}->{$param} ne $old_config->{$module}->{$param} ){
				foreach my $task ( keys %{$module_ctl->{$module}} ){
					`rm -f $module_ctl->{$module}->{$task}`;
				}
				last;
			}
		}
	}

	#remove old config file
	`rm -f "$proj_dir/config.txt.bak"`;

	return;
}
