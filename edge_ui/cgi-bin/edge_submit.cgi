#!/usr/bin/env perl
#
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-07
#

use strict;
use FindBin qw($RealBin);
use Cwd 'abs_path';
use lib "$RealBin/../../lib";
use HTML::Template;
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use CGI::Session ( '-ip_match' );
use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser);
use POSIX qw(strftime);
use Digest::MD5 qw(md5_hex);
#use Data::Dumper;
require "edge_user_session.cgi";
require "../cluster/clusterWrapper.pl";

my $cgi = CGI->new;
my %opt = $cgi->Vars();
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

# running vars
my $msg;
my $excode;

my $projConfigTmpl  = "$RealBin/edge_config.tmpl";
# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $um_url      = $sys->{edge_user_management_url};
my $protocol    = $opt{protocol};
my $domain      = $ENV{'HTTP_HOST'};
$um_url ||= "$protocol//$domain/userManagement";
my $debug       = $sys->{debug};

#cluster
my $cluster 	= $sys->{cluster};
my $cluster_qsub_options= $sys->{cluster_qsub_options};
my $cluster_job_resource= $sys->{cluster_job_resource};
my $cluster_job_max_cpu= $sys->{cluster_job_max_cpu};
my $cluster_job_notify = $sys->{cluster_job_notify};
my $cluster_job_prefix = $sys->{cluster_job_prefix};
my $cluster_tmpl = "$RealBin/../cluster/clusterSubmit.tmpl";

#init vars
my $pipeline	= $opt{pipeline};
my $pname       = $opt{'edge-proj-name'} || $ARGV[0];
my $edge_input	= $sys->{edgeui_input};
my $input_dir   = "$edge_input/public";
my $edge_output	= $sys->{edgeui_output};
my $out_dir     = $opt{"edge-proj-outpath"};
$out_dir      ||= $sys->{edgeui_output};
my $edge_total_cpu = $sys->{"edgeui_tol_cpu"};
$edge_total_cpu = $cluster_job_max_cpu if ($cluster);
my $num_cpu     = $opt{"edge-proj-cpu"};
$num_cpu      ||= 4;
my $username 	= $opt{'username'} || $ARGV[1];
my $password 	= $opt{'password'} || $ARGV[2];
my $sid         = $opt{'sid'} || $ARGV[3];

#my $proj_dir    = "$out_dir/$opt{'edge-proj-name'}";
#my $config_out  = "$proj_dir/config.txt";
my $ref_list    = "$RealBin/../data/Ref_list.json";
my $host_list   = "$RealBin/../data/host_list.json";
my $pid;
my $reflist        = &readListFromJson($ref_list);
my $hostlist       = &readListFromJson($host_list);
my @edge_input_pe1 = split /[\x0]/, $opt{"edge-input-pe1[]"};
my @edge_input_pe2 = split /[\x0]/, $opt{"edge-input-pe2[]"};
my @edge_input_se  = split /[\x0]/, $opt{"edge-input-se[]"};
my @edge_qiime_mapping_files = split /[\x0]/, $opt{"edge-qiime-mapping-file-input[]"};
my $edge_qiime_input_dir = $opt{"edge-qiime-reads-dir-input"};
my @edge_qiime_barcode_input;
my @edge_phylo_ref_input;
#check session
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	unless($valid){
		addMessage("RUN_PIPELINE","failure","Invalid session found.");
		returnStatus();
	}
	else{
		($username,$password) = getCredentialsFromSession($sid);
		$input_dir="$edge_input/". md5_hex($username);
	}
}

# Qiime dir submit
if ($edge_qiime_input_dir){
	my ($pe1_file_r,$pe2_file_r,$se_file_r) = &parse_qiime_mapping_files($edge_qiime_input_dir,\@edge_qiime_mapping_files);
	@edge_input_pe1 = @{$pe1_file_r};
	@edge_input_pe2 = @{$pe2_file_r};
	@edge_input_se = @{$se_file_r};
	$opt{"edge-qiime-reads-dir-input"}="$input_dir/$edge_qiime_input_dir";
}

# batch submit
my $projlist;
my @pnames;
if ($opt{"edge-batch-text-input"}){
	$projlist=&readBatchInput();
	@pnames = keys %{$projlist};
}else{
	push @pnames, $pname ;
}
my @real_names = @pnames;



#init SUBMISSION_STATUS
$msg->{SUBMISSION_STATUS}="success";

# validating parameters
&checkParams();

if ($msg->{SUBMISSION_STATUS} eq 'success'){
	#if user management is on. will replace pnames to unique IDs
	if ($username && $password){
		&addProjToDB();
	}else{
		foreach (@pnames){
			`ln -sf $edge_output/$_ $input_dir/projects/$_` if ( ! -e "$input_dir/projects/$_");
		}
	}
}

my ($vital, $name2pid, $error);
if($cluster) {
	($vital, $name2pid, $error) = checkProjVital_cluster($cluster_job_prefix);
	if($error) {
		&addMessage("CLUSTER","failure",$error);
	}
} else {
	($vital, $name2pid) = &checkProjVital();
}

# check no running duplications
&checkRunningProject() if $msg->{SUBMISSION_STATUS} eq 'success';

# prepare to run pipeline
&createProjDir() if $msg->{SUBMISSION_STATUS} eq 'success';

# re-use config or generate one
&createConfig() if $msg->{SUBMISSION_STATUS} eq 'success';

# run pipeline
if($cluster) {
	&runPipeline_cluster() if $msg->{SUBMISSION_STATUS} eq 'success';
} else {
	&runPipeline() if $msg->{SUBMISSION_STATUS} eq 'success';
}

# return
&returnStatus();

############################################################################
sub readBatchInput {
	my $list={};
	my $proj_name;
	foreach (split /\n/, $opt{"edge-batch-text-input"}) {
        	if (/^\#/) { next; }
        	unless (/\w/) { next; }
        	s/^\s+|\s+$//g;
        	chomp;
        	if (/^\[(.*)\]/) { $proj_name=$1; next; }
		my ($key, $val) = split(/=/, $_);
		if ($key ne "description"){
			$val =~ s/ //g;
		}
	        $val="" if (!$val);
        	$list->{$proj_name}->{$key} = $val;
    	}
	return $list;
}

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			while(<CONF>){
				chomp;
				last if /^\[/;
				next if(/^#/);
				if ( /^([^=]+)=(.*)/ ){
					$sys->{$1}=$2;
				}
			}
		}
		last;
	}
	close CONF;
	return $sys;
}

sub createProjDir {
	foreach my $pname (@pnames){
    		my $proj_dir = "$out_dir/$pname";
		$proj_dir = "$out_dir/" . $projlist->{$pname}->{projCode} if ($username && $password);
		#clean the empty dir before checking the conflicted output path
    		system("rmdir $proj_dir") if ( -d $proj_dir && is_folder_empty($proj_dir));
    
		#detect conflicted output path
		if( -l "$sys->{edgeui_output}/$pname" || -d "$sys->{edgeui_output}/$pname"){
			if( abs_path("$sys->{edgeui_output}/$pname") eq $proj_dir ){
				&addMessage("CREATE_OUTPUT","failure","Conflicting output path with existing project ($pname).");
				return;
			}
		}
	}
	foreach my $pname (@pnames){
		my $proj_dir = "$out_dir/$pname";
		$proj_dir = "$out_dir/" . $projlist->{$pname}->{projCode} if ($username && $password);
		#init output directory
		$excode = system("mkdir -m 755 -p $proj_dir");

		if( $excode || ! -d $proj_dir ){
			&addMessage("CREATE_OUTPUT","failure","FAILED to create output directory.");
			return;
		}
		else{
			&addMessage("CREATE_OUTPUT","success","Output directory created.");
		}

		#the user specified output path
		if ($opt{"edge-proj-outpath"}){
		#if( $proj_dir ne "$sys->{edgeui_output}/$pname" ){
			$excode = system("ln -sf $proj_dir $sys->{edgeui_output}/$opt{'edge-proj-name'}");
			if( $excode ){
				&addMessage("CREATE_OUTPUT","failure","FAILED to create symlink to $sys->{edgeui_output}.");
				return;
			}
		}
	}
}

sub is_folder_empty {
	my $dirname = shift;
	opendir(my $dh, $dirname) or die "Not a directory";
	return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

sub createConfig {
	foreach my $pname (@pnames){
		my $config_out = "$out_dir/$pname/config.txt";
		$config_out = "$out_dir/" . $projlist->{$pname}->{projCode} . "/config.txt" if ($username && $password);
		if( defined $opt{"edge-input-config"} && -e $opt{"edge-input-config"} ){
			open CFG, $opt{"edge-input-config"};
			open CFG_OUT, ">$config_out";
			while(<CFG>){
				chomp;
				my $line = $_;
				if( $line =~ /^cpu=/ ){
					$line =~ s/^cpu=.*/cpu=$num_cpu/;
				}
				elsif( $line =~ /^outpath=/ ){
					$line =~ s/^outpath=.*/outpath=$out_dir/;
				}
				elsif( $line =~ /^projname=/ ){
					$line =~ s/^projname=.*/projname=$opt{"edge-proj-name"}/;
				}
				elsif( $line =~ /^projdesc=/ ){
					$line =~ s/^projdesc=.*/projdesc=$opt{"edge-proj-desc"}/;
				}
				elsif( $line =~ /^projid=/ ){
					$line =~ s/^projdesc=.*/projdesc=$pname/;
				}
				print CFG_OUT "$line\n";
			}
			close CFG_OUT;
			close CFG;
			&addMessage("CONFIG","info","Use existing config $config_out.");
		}
		else{
			if ( $opt{"edge-hostrm-sw"} )
			{
				my %hosts;
				map { foreach my $h (@{$hostlist->{$_}}){$hosts{$h}=1;} } split /[\x0]/, $opt{"edge-hostrm-file-fromlist"} if defined $opt{"edge-hostrm-file-fromlist"};
				$hosts{$opt{"edge-hostrm-file"}}=1 if (defined $opt{"edge-hostrm-file"} && -e $opt{"edge-hostrm-file"});
				$opt{"edge-hostrm-file"} = join ",", keys %hosts;
			} 

			if ( $opt{"edge-ref-sw"} )
			{
				my (@refs,@refsl);
				if( defined $opt{"edge-ref-file-fromlist"}){
					@refsl = split /[\x0]/, $opt{"edge-ref-file-fromlist"};
					map { 	my @tmp= `ls -d $EDGE_HOME/database/NCBI_genomes/$_*`; 
						chomp @tmp;  
						my @gfiles = `ls $tmp[0]/*gbk`; 
						chomp @gfiles; 
						push @refs,@gfiles;
					    } @refsl;
				}
				push @refs, $opt{"edge-ref-file"} if (defined $opt{"edge-ref-file"} && -e $opt{"edge-ref-file"});
				$opt{"edge-ref-file"} = join ",", @refs;
			}
			if ($opt{"edge-phylo-sw"})
			{
				my @snpPhylo_refs;
				map { push @snpPhylo_refs, $_} split /[\x0]/, $opt{"edge-phylo-ref-select"} if defined $opt{"edge-phylo-ref-select"};
				$opt{"edge-phylo-ref-list"} = join ",",@snpPhylo_refs if @snpPhylo_refs;
				$opt{"edge-phylo-ref-list-file"} = join ",",@edge_phylo_ref_input if @edge_phylo_ref_input;
			}
	   
			$opt{"edge-taxa-enabled-tools"} =~ s/[\x0]/,/g if $opt{"edge-taxa-sw"};
			$opt{'edge-sra-acc'} = uc $opt{'edge-sra-acc'};
			$opt{'edge-phylo-sra-acc'} = uc $opt{'edge-phylo-sra-acc'};
			$opt{'edge-qiime-mapping-files'} = join ",", @edge_qiime_mapping_files if @edge_qiime_mapping_files;
			$opt{'edge-qiime-barcode-fq-files'} = join ",", @edge_qiime_barcode_input if @edge_qiime_barcode_input;

      		        $opt{"edge-proj-desc"} = $projlist->{$pname}->{"description"} if ($opt{"edge-batch-text-input"});
      		        $opt{"edge-proj-name"} = $projlist->{$pname}->{"REALNAME"}||$pname if ($opt{"edge-batch-text-input"});
			$opt{"edge-proj-id"} = $pname;
			$opt{"edge-proj-code"} = $projlist->{$pname}->{projCode};
			$opt{"edge-proj-runhost"} = "$protocol//$domain";
			$opt{"edge-proj-owner"} = $username if ($username);

			eval {
				my $template = HTML::Template->new(filename => $projConfigTmpl, die_on_bad_params => 0 );
				$template->param( %opt );
			
				open CONFIG, ">$config_out" or die "Can't write config file: $!";
				print CONFIG $template->output();
				close CONFIG;
			
			};
			if( $@ ){
				my $error_msg = join " ", $@;
				&addMessage("GENERATE_CONFIG","failure",$error_msg);
			}
			else{
				&addMessage("GENERATE_CONFIG","success","Config file generated.");
			}
		}
	}
}

sub ref_extract
{
	my $glist = shift;
	my $pname = shift;
	my $cmd = "$EDGE_HOME/scripts/extract_fasta_from_BWA_idx.pl --genomesList $glist -o '$out_dir/$pname/Reference'";
	my @files = `$cmd`;
	chomp(@files);
	return @files;
}

sub runPipeline {
	my $proj_count=1;
	foreach my $pname (@pnames){
		my $proj_dir = "$out_dir/$pname";
		$proj_dir = "$out_dir/" . $projlist->{$pname}->{projCode} if ($username && $password);
		my $config_out = "$proj_dir/config.txt";
		my ($paired_files, $single_files) = ("","");
		my $process_parameters = "-c $config_out -o $proj_dir -cpu $num_cpu -noColorLog ";
    	
 	   	if ($opt{"edge-batch-text-input"}){
    			$paired_files = qq|$projlist->{$pname}->{"q1"} $projlist->{$pname}->{"q2"}| if ( -f $projlist->{$pname}->{"q1"});
    			$single_files = $projlist->{$pname}->{"s"} if ( -f $projlist->{$pname}->{"s"});
    		}else{
			for (0..$#edge_input_pe1)
			{
				$paired_files .= "$edge_input_pe1[$_] $edge_input_pe2[$_] ";
			}
		
			$single_files = join " ", @edge_input_se;
		}
		$process_parameters .= " --debug " if ($debug);
		$process_parameters .= " -p $paired_files " if ($paired_files);
		$process_parameters .= " -u $single_files " if ($single_files);

		my $cmd = "$EDGE_HOME/runPipeline $process_parameters > $proj_dir/process_current.log &";

		my $run = &availableToRun($num_cpu*$proj_count);
		if( $run ){
			chdir($proj_dir);
			$pid = open RUNPIPLINE, "-|", $cmd or die $!;
			close RUNPIPLINE;

			print STDERR "runPipeline command: $cmd\n";

			unless( $pid ){
				&addMessage("RUN_PIPELINE","failure","FAILED to start pipeline.");
			}
			else{
				my $real_name = ($username && $password)? $projlist->{$pname}->{REALNAME}:$pname;
				&addMessage("RUN_PIPELINE","success","Project $real_name has been submitted successfully.");
			}
		}
		else{
			my $time = strftime "%F %X", localtime;
			`echo "\n*** [$time] EDGE_UI: This project is queued. ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
			`echo "$cmd" >> $proj_dir/process.log`;
			`echo "\n*** [$time] EDGE_UI: Project unstarted ***" >> $proj_dir/process.log`;
			&addMessage("RUN_PIPELINE","success","Project $opt{'edge-proj-name'} has been added but unstarted.");
			&addMessage("SYSTEM","success","The server does not have enough CPU available to run this job.");
		}
		$proj_count++;
	}
}

sub runPipeline_cluster {
	my $proj_count=1;
	if($num_cpu > $cluster_job_max_cpu) {
		$num_cpu = $cluster_job_max_cpu;
	}

	if($cluster_job_resource =~ /<CPU\/2>/) {
		#uge
		if($num_cpu %2 == 1) {
			$num_cpu --;
		}

		my $binding_cpu = $num_cpu/2;
		$cluster_job_resource =~ s/<CPU\/2>/$binding_cpu/;
	}
	$cluster_job_resource =~ s/<CPU>/$num_cpu/;

	foreach my $pname (@pnames){
		my $job_name = $cluster_job_prefix.$pname;
		my $proj_dir = "$out_dir/$pname";
		$proj_dir = "$out_dir/" . $projlist->{$pname}->{projCode} if ($username && $password);
		my $config_out = "$proj_dir/config.txt";
		my $cluster_job_script = "$proj_dir/clusterSubmit.sh";
		my $cluster_job_log = "$proj_dir/clusterJob.log";
		my ($paired_files, $single_files) = ("","");
		my $process_parameters = "-c $config_out -o $proj_dir -cpu $num_cpu -noColorLog ";
    	
 	   	if ($opt{"edge-batch-text-input"}){
    			$paired_files = qq|$projlist->{$pname}->{"q1"} $projlist->{$pname}->{"q2"}| if ( -f $projlist->{$pname}->{"q1"});
    			$single_files = $projlist->{$pname}->{"s"} if ( -f $projlist->{$pname}->{"s"});
    		}else{
			for (0..$#edge_input_pe1)
			{
				$paired_files .= "$edge_input_pe1[$_] $edge_input_pe2[$_] ";
			}
		
			$single_files = join " ", @edge_input_se;
		}
		$process_parameters .= " --debug " if ($debug);
		$process_parameters .= " -p $paired_files " if ($paired_files);
		$process_parameters .= " -u $single_files " if ($single_files);

		my $cmd = "$EDGE_HOME/runPipeline $process_parameters > $proj_dir/process_current.log";

		open CT, $cluster_tmpl;
		open CT_OUT, ">$cluster_job_script";
		while(<CT>) {
			chomp;	
			if(/<JOB_NAME>/) {
				s/<JOB_NAME>/$job_name/;
			} elsif (/<JOB_RESOURCE_REQUEST>/) {
				s/<JOB_RESOURCE_REQUEST>/$cluster_job_resource/;
			} elsif(/<JOB_NOTIFY>/) {
				s/<JOB_NOTIFY>/$cluster_job_notify/;
			} elsif (/<JOB_LOG>/) {
				s/<JOB_LOG>/$cluster_job_log/;
			} elsif (/<COMMAND>/) {
				s/<COMMAND>/$cmd/;
			}
			print CT_OUT "$_\n";
		}
		close CT_OUT;
		close CT;
		&addMessage("CLUSTER","info","Create job script: $cluster_job_script");
		
		my ($job_id,$error) = clusterSubmitJob($cluster_job_script,$cluster_qsub_options);
		if($error) {
			&addMessage("CLUSTER","failure","FAILED to submit $cluster_job_script: $error");
		} else {
			&addMessage("CLUSTER","info","$cluster_job_script job id: $job_id");
		}
		$proj_count++;
	}
}

sub addProjToDB{
	my $desc; 
	foreach my $pname_index (0..$#pnames){
		my $pname = $pnames[$pname_index];
		if ($opt{"edge-batch-text-input"}){
			$desc = $projlist->{$pname}->{description};
		}else{
			$desc = $opt{"edge-proj-desc"};
		}
		$desc =~ s/(['"])/\\$1/g;

		my %data = (
			email => $username,
   			password => $password,
   			project_name => $pname,
   			description => $desc,
		);
		# Encode the data structure to JSON
		my $data = to_json(\%data);
		#w Set the request parameters
		my $url = $um_url ."WS/project/add";
		my $browser = LWP::UserAgent->new;
		my $req = POST $url;
		$req->header('Content-Type' => 'application/json');
		$req->header('Accept' => 'application/json');
		#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
		$req->header( "Content-Length" => length($data) );
		$req->content($data);

		my $response = $browser->request($req);
		my $result_json = $response->decoded_content;
		#print $result_json,"\n";
		my $info =  from_json($result_json);
		my $new_id =  $info->{"id"};
		my $projCode = &getProjcode($new_id);
		$pnames[$pname_index] = $new_id;
		&addMessage("PROJECT_NAME","info","Assigned project $pname with ID $new_id");
		&addMessage("ASSIGN_PROJ_ID","failure","Database doens't assign $pname a project ID. You may not be logged in properly. Please contact admins.") if (!$new_id);
		# assign value to new project id.
		if ($opt{"edge-batch-text-input"}){
			foreach my $key (keys %{$projlist->{$pname}}){
				$projlist->{$new_id}->{$key} = $projlist->{$pname}->{$key}; 
				$projlist->{$new_id}->{projCode} = $projCode; 
			}
		}

		$projlist->{$new_id}->{REALNAME} = $pname; 
		$projlist->{$new_id}->{projCode} = $projCode; 

		my $user_project_dir = "$input_dir/MyProjects/${pname}_$new_id";
		#`ln -sf $edge_output/$new_id $user_project_dir` if (! -e $user_project_dir);
		`ln -sf $edge_output/$projCode $user_project_dir` if (! -e $user_project_dir);
	}
}

sub getProjcode {
	my $id = shift;
	my %data = (
		email => $username,
		password => $password,
		project_id => $id
	);
	# Encode the data structure to JSON		
	my $data =  encode_json(\%data);

	# Set the request parameters
	my $url = $um_url. 'WS/project/getInfo';
	my $browser = LWP::UserAgent->new;
	my $req = PUT $url;
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	my $info =  from_json($result_json);
	return  $info->{"code"};
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

sub returnStatus {
	# return json
	my $json = encode_json($msg);
	print $cgi->header('application/json'), $json;
	exit 0;
}

sub readListFromJson {
	my $json = shift;
	my $list = {};
	if( -r $json ){
		open JSON, $json;
		flock(JSON, 1);
  		local $/ = undef;
  		$list = decode_json(<JSON>);
  		close JSON;
	}
	return $list;
}

sub checkRunningProject { 
	foreach my $pid ( keys %$vital ){
		my $running_pname = $vital->{$pid}->{PROJ};
		if( $running_pname eq $pname || defined $projlist->{$running_pname}){
			my $real_name = ($username && $password)? $projlist->{$running_pname}->{REALNAME}:$running_pname;
			&addMessage("RUN_PIPELINE","failure","Project $real_name is already running.");
		}
	}
}

sub addMessage {
	my ($job, $status, $note) = @_;
	if( $job eq "PARAMS" ){
		#e.g.: $msg->{"PARAMS"}->{"edge-ref-file"}="File not found";
		$msg->{$job}->{$status}=$note;
		$msg->{SUBMISSION_STATUS}="failure";
	}
	else{
		#&addMessage("PROJECT_NAME","info","$info");
		$msg->{$job}->{STATUS}=$status;
		$msg->{$job}->{NOTE}=$note if defined $note;
		$msg->{SUBMISSION_STATUS}=$status if $status eq "failure";
	}
}

sub checkParams {
	my %files;		
	if ($num_cpu > $edge_total_cpu){
		&addMessage("PARAMS","edge-proj-cpu","The max number of CPU for the EDGE Server is $edge_total_cpu.");
	}
	if ( ($opt{"edge-batch-text-input"}) and ($opt{"edge-proj-desc"} or $opt{"edge-input-pe1[]"} or $opt{"edge-input-pe2[]"} or $opt{"edge-input-se[]"} or $opt{"edge-sra-acc"})){
		&addMessage("PARAMS","edge-input-sequence","Input error. You have both single project input and batch input.");
    	}
      
	if (  $opt{"edge-batch-text-input"} ){ ## batch input
    		my %namesUsed;
		foreach my $pname (keys %{$projlist}){
			$projlist->{$pname}->{"q1"} = "$input_dir/$projlist->{$pname}->{'q1'}";
    			$projlist->{$pname}->{"q2"} = "$input_dir/$projlist->{$pname}->{'q2'}";
    			$projlist->{$pname}->{"s"} = "$input_dir/$projlist->{$pname}->{'s'}";
    			my $pe1=$projlist->{$pname}->{"q1"};
    			my $pe2=$projlist->{$pname}->{"q2"};
    			my $se=$projlist->{$pname}->{"s"};

    			if ($namesUsed{$pname}){
    				&addMessage("PARAMS","edge-batch-text-input", "Duplicate project name found.");
    			}else{
    				$namesUsed{$pname}=1;
    			}
    			&addMessage("PARAMS","edge-batch-text-input","Invalid project name. Only alphabets, numbers and underscore are allowed in project name.") if ($pname =~ /\W/);
    			&addMessage("PARAMS","edge-batch-text-input","Invalid project name. Please input at least 3 characters.") if (length($pname) < 3);
    			&addMessage("PARAMS","edge-batch-text-input","Invalid characters detected in $pe1 of $pname.") if (-f $pe1 and $pe1 =~ /[\<\>\!\~\@\#\$\^\&\;\*\(\)\"\' ]/);
    			&addMessage("PARAMS","edge-batch-text-input","Invalid characters detected in $pe2 of $pname.") if (-f $pe2 and $pe2 =~ /[\<\>\!\~\@\#\$\^\&\;\*\(\)\"\' ]/);
    			&addMessage("PARAMS","edge-batch-text-input","Invalid characters detected in $pe2 of $pname.") if (-f $se and $se =~ /[\<\>\!\~\@\#\$\^\&\;\*\(\)\"\' ]/);
    			&addMessage("PARAMS","edge-batch-text-input","Input error. Please check the q1 file path of $pname.") if (-f $pe1 && $pe1 !~ /^[http|ftp]/i && ! -e $pe1);
    			&addMessage("PARAMS","edge-batch-text-input","Input error. Please check the q2 file path of $pname.") if (-f $pe2 && $pe2 !~ /^[http|ftp]/i && ! -e $pe2);
    			&addMessage("PARAMS","edge-batch-text-input","Input error. q1 and q2 are identical of $pname.") if ( -f $pe1 && $pe1 eq $pe2);
    			&addMessage("PARAMS","edge-batch-text-input","Input error. Please check the s file path of $pname.") if (-f $se && $se !~ /^[http|ftp]/i && ! -e $se);
    			&addMessage("PARAMS","edge-batch-text-input","Input error. Please check the input file path of $pname.") if (! -f $se && ! -f $pe1 && ! -f $pe2);
    		}
	}else{  ## Single project input
		&addMessage("PARAMS","edge-proj-name","Invalid project name. Only alphabets, numbers, dashs, dot and underscore are allowed in project name.") if( $opt{"edge-proj-name"} =~ /[^a-zA-Z0-9\-_\.]/ );
		&addMessage("PARAMS","edge-proj-name","Invalid project name. Please input at least 3 characters.") if( length($opt{"edge-proj-name"}) < 3 );
		#check invalid character
		foreach my $param (keys %opt ){
			next if $param eq "edge-proj-desc";
			next if $param eq "edge-batch-text-input";
			next if $param eq "username";
			next if $param eq "password";
			next if $param =~ /aligner-options/;
			&addMessage("PARAMS","$param","$param Invalid characters detected.") if $opt{$param} =~ /[\<\>\!\~\@\#\$\^\&\;\*\(\)\"\' ]/;
		}
		
		#check input sequence
		foreach my $i (0..$#edge_input_pe1){
			my $id = "edge-input-pe1-". ($i + 1);
			$edge_input_pe1[$i] =~ s/ //g;
			$edge_input_pe1[$i] = "$input_dir/$edge_input_pe1[$i]" if ($edge_input_pe1[$i] =~ /^\w/);
			&addMessage("PARAMS","$id","Error: duplicated input.") if ($files{$edge_input_pe1[$i]});
			$files{$edge_input_pe1[$i]}=1;
			if ($edge_input_pe1[$i] &&  $edge_input_pe1[$i] !~ /^[http|ftp]/i  && ! -e $edge_input_pe1[$i]){
				&addMessage("PARAMS","$id","Input error. Please check the file path.");
			}
			 &addMessage("PARAMS","$id","Input error. Fastq format required ") unless ( is_fastq($edge_input_pe1[$i]) );
		}
		foreach my $i (0..$#edge_input_pe2){
			my $id = "edge-input-pe2-". ($i + 1);
			$edge_input_pe2[$i] =~ s/ //g;
			$edge_input_pe2[$i] = "$input_dir/$edge_input_pe2[$i]" if ($edge_input_pe2[$i] =~ /^\w/);
			&addMessage("PARAMS","$id","Error: duplicated input.") if ($files{$edge_input_pe2[$i]});
			$files{$edge_input_pe2[$i]}=1;
			if ($edge_input_pe2[$i] &&  $edge_input_pe2[$i] !~ /^[http|ftp]/i  && ! -e $edge_input_pe2[$i]){
				&addMessage("PARAMS","$id","Input error. Please check the file path.");
			}
			 &addMessage("PARAMS","$id","Input error. Fastq format required ") unless ( is_fastq($edge_input_pe2[$i]) );
		}
		foreach my $i (0..$#edge_input_se){
			my $id = "edge-input-se". ($i + 1);
			$edge_input_se[$i] =~ s/ //g;
			$edge_input_se[$i] = "$input_dir/$edge_input_se[$i]" if ($edge_input_se[$i] =~ /^\w/);
			&addMessage("PARAMS","$id","Error: duplicated input.") if ($files{$edge_input_se[$i]});
			$files{$edge_input_se[$i]}=1;
			if ($edge_input_se[$i] &&  $edge_input_se[$i] !~ /^[http|ftp]/i  && ! -e $edge_input_se[$i]){
				&addMessage("PARAMS","$id","Input error. Please check the file path.");
			}
			 &addMessage("PARAMS","$id","Input error. Fastq format required ") unless ( is_fastq($edge_input_se[$i]) );
		}
		foreach my $i (0..$#edge_input_pe1){
			my $id = "edge-input-pe-block".($1 + 1);
			if( (-e $edge_input_pe1[$i] && $edge_input_pe2[$i] !~ /^[http|ftp]/i  && ! -e $edge_input_pe2[$i]) || ($edge_input_pe1[$i] !~ /^[http|ftp]/i  && !-e $edge_input_pe1[$i] && -e $edge_input_pe2[$i]) ){
				&addMessage("PARAMS","$id","Input error. Please check the file path.");
			}
			if ($edge_input_pe1[$i] eq $edge_input_pe2[$i]){
				&addMessage("PARAMS","$id","Input error. Pair-1 and Pair-2 are identical.");
			}
		}
		if ( $opt{'edge-sra-sw'}){
			&addMessage("PARAMS","edge-sra-acc","Input error. Please input SRA accession") if ( ! $opt{'edge-sra-acc'});
		}else{
			if (!@edge_input_pe1 && !@edge_input_pe2 && !@edge_input_se){
				&addMessage("PARAMS","edge-input-pe1-1","Input error.");
				&addMessage("PARAMS","edge-input-pe2-1","Input error.");
				&addMessage("PARAMS","edge-input-se1","Input error.");
			}
			if ($pipeline eq "qiime" && @edge_input_pe1 && @edge_input_se){

				&addMessage("PARAMS","edge-input-se1","Input error. Please provide either paired-end Or single-end fastq.");
			}
		}
		if ($pipeline eq "qiime"){
			foreach my $i (0..$#edge_qiime_mapping_files){
				my $id = "edge-qiime-mapping-file-input". ($i + 1);
				$edge_qiime_mapping_files[$i] =~ s/ //g;
				$edge_qiime_mapping_files[$i] = "$input_dir/$edge_qiime_mapping_files[$i]" if ($edge_qiime_mapping_files[$i] =~ /^\w/);
				&addMessage("PARAMS","$id","Error: duplicated input.") if ($files{$edge_qiime_mapping_files[$i]});
				$files{$edge_qiime_mapping_files[$i]}=1;
				if ($edge_qiime_mapping_files[$i] &&  $edge_input_se[$i] !~ /^[http|ftp]/i  && ! -e $edge_qiime_mapping_files[$i]){
					&addMessage("PARAMS","$id","Input error. Please check the file path.");
				}
			}
		}
	}
	
	#tool parameters
	if ( $opt{"edge-ref-sw"}){
		$opt{"edge-ref-file"} = $input_dir."/".$opt{"edge-ref-file"} if ( $opt{"edge-ref-file"} =~ /^\w/ );
		&addMessage("PARAMS","edge-ref-file","Reference not found. Please check the input referecne.") if( !-e $opt{"edge-ref-file"} && !defined $opt{"edge-ref-file-fromlist"});
		&addMessage("PARAMS","edge-ref-file-fromlist","Reference not found. Please check the input referecne.") if( !-e $opt{"edge-ref-file"} && !defined $opt{"edge-ref-file-fromlist"});
		&addMessage("PARAMS","edge-ref-file","Invalid input. Fasta or Genbank format required") if ( -e $opt{"edge-ref-file"} && ! is_fasta($opt{"edge-ref-file"}) && ! is_genbank($opt{"edge-ref-file"}) );
	}
	if ( $opt{"edge-taxa-sw"} && scalar split(/[\x0]/,$opt{"edge-taxa-enabled-tools"}) < 1 ){
		&addMessage("PARAMS","edge-taxa-tools","You need to choose at least one tool.");
	}
	&addMessage("PARAMS","edge-primer-vaild-sw","Input primer fasta file. You need to turn on the Primer Validation") if ($opt{"edge-primer-valid-file"} && !$opt{"edge-primer-valid-sw"});
	if ( $opt{"edge-primer-valid-sw"} ){
		$opt{"edge-primer-valid-file"} = $input_dir."/".$opt{"edge-primer-valid-file"} if ($opt{"edge-primer-valid-file"} =~ /^\w/);
		&addMessage("PARAMS","edge-primer-valid-file","File not found. Please check the primer FASTA file.") unless -e $opt{"edge-primer-valid-file"};
		&addMessage("PARAMS","edge-primer-valid-file","Invalid input. Fasta format required") if ( -e $opt{"edge-primer-valid-file"} && ! is_fasta($opt{"edge-primer-valid-file"}));
		&addMessage("PARAMS","edge-primer-valid-mm","Invalid input. Natural number required.")     unless $opt{"edge-primer-valid-mm"} =~ /^\d+$/;
	}
	if ( $opt{"edge-primer-adj-sw"} ){
		&addMessage("PARAMS", "edge-primer-adj-df", "Invalid input. Natural number required.")     unless $opt{"edge-primer-adj-df"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-primer-adj-num","Invalid input. Natural number required.")     unless $opt{"edge-primer-adj-num"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-primer-adj-tm-opt", "Invalid input. Input should in range 0-100 C.") unless ( $opt{"edge-primer-adj-tm-opt"} >= 0 && $opt{"edge-primer-adj-tm-opt"} <=100 );
		&addMessage("PARAMS", "edge-primer-adj-tm-max", "Invalid input. Input should in range 0-100 C.") unless ( $opt{"edge-primer-adj-tm-max"} >= 0 && $opt{"edge-primer-adj-tm-max"} <=100 );
		&addMessage("PARAMS", "edge-primer-adj-tm-min", "Invalid input. Input should in range 0-100 C.") unless ( $opt{"edge-primer-adj-tm-min"} >= 0 && $opt{"edge-primer-adj-tm-min"} <=100 );
		&addMessage("PARAMS", "edge-primer-adj-tm-opt", "Invalid input. Tm is not in the range.") if $opt{"edge-primer-adj-tm-opt"} > $opt{"edge-primer-adj-tm-max"} || $opt{"edge-primer-adj-tm-opt"} < $opt{"edge-primer-adj-tm-min"};
		&addMessage("PARAMS", "edge-primer-adj-len-opt", "Invalid input. Natural number required.") unless $opt{"edge-primer-adj-len-opt"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-primer-adj-len-max", "Invalid input. Natural number required.") unless $opt{"edge-primer-adj-len-max"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-primer-adj-len-min", "Invalid input. Natural number required.") unless $opt{"edge-primer-adj-len-min"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-primer-adj-len-opt", "Invalid input. Length is not in the range.") if $opt{"edge-primer-adj-len-opt"} > $opt{"edge-primer-adj-len-max"} || $opt{"edge-primer-adj-len-opt"} < $opt{"edge-primer-adj-len-min"};
	}
	if ( $opt{"edge-qc-sw"} ){
		&addMessage("PARAMS", "edge-qc-q",          "Invalid input. Natural number required.")     unless $opt{"edge-qc-q"}    =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-avgq",       "Invalid input. Natural number required.")     unless $opt{"edge-qc-avgq"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-minl",       "Invalid input. Natural number required.")     unless $opt{"edge-qc-minl"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-n",          "Invalid input. Natural number required.")     unless $opt{"edge-qc-n"}    =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-lc",         "Invalid input. Floating number between 0 and 1 required.") unless ( $opt{"edge-qc-lc"} >=0 && $opt{"edge-qc-lc"} <=1 );
		$opt{"edge-qc-adapter"} = $input_dir."/".$opt{"edge-qc-adapter"} if ($opt{"edge-qc-adapter"} =~ /^\w/);
		&addMessage("PARAMS", "edge-qc-adapter",    "File not found. Please check the file path.") if ( $opt{"edge-qc-adapter"} && ! -e $opt{"edge-qc-adapter"} );
		#&addMessage("PARAMS", "edge-qc-phix",       "Invalid input.")                              unless $opt{"edge-qc-phix"} =~ /(0|1)/;
		&addMessage("PARAMS", "edge-qc-5end",       "Invalid input. Natural number required.")     unless $opt{"edge-qc-5end"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-3end",       "Invalid input. Natural number required.")     unless $opt{"edge-qc-3end"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-qc-adapter",    "Invalid input. Fasta format required") if ( -e $opt{"edge-qc-adapter"} and ! is_fasta($opt{"edge-qc-adapter"}) );
	}
	if ( $opt{"edge-hostrm-sw"} ){
		$opt{"edge-hostrm-file"} = $input_dir."/".$opt{"edge-hostrm-file"} if ($opt{"edge-hostrm-file"} =~ /^\w/);
		&addMessage("PARAMS", "edge-hostrm-file",   "File not found. Please check the host input.") if ( !defined $opt{"edge-hostrm-file-fromlist"} && !-e $opt{"edge-hostrm-file"} );
		&addMessage("PARAMS", "edge-hostrm-file",   "Invalid input. Fasta format required") if ( -e $opt{"edge-hostrm-file"} && ! is_fasta($opt{"edge-hostrm-file"}));
		&addMessage("PARAMS", "edge-hostrm-similarity", "Invalid input. Natural number required.")     unless $opt{"edge-hostrm-similarity"} =~ /^\d+$/;
	}
	if ( $opt{"edge-assembly-sw"} ){
		&addMessage("PARAMS", "edge-assembly-mink", "Invalid input. Natural number required.")     unless $opt{"edge-assembly-mink"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-assembly-maxk", "Invalid input. Natural number required.")     unless $opt{"edge-assembly-maxk"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-assembly-step", "Invalid input. Natural number required.")     unless $opt{"edge-assembly-step"} =~ /^\d+$/;
		&addMessage("PARAMS", "edge-assembly-minc", "Invalid input. Natural number required.")     unless $opt{"edge-assembly-minc"} =~ /^\d+$/;
		if ($opt{"edge-assembled-contig-sw"}){
			$opt{"edge-assembled-contig-file"} = $input_dir."/".$opt{"edge-assembled-contig-file"} if ($opt{"edge-assembled-contig-file"} =~ /^\w/);
			&addMessage("PARAMS","edge-assembled-contig-file", "File not found. Please check the contig file input.") if ( ! -e $opt{"edge-assembled-contig-file"});
			&addMessage("PARAMS","edge-assembled-contig-file", "Invalid input. Fasta format required.") if ( -e $opt{"edge-assembled-contig-file"} && ! is_fasta($opt{"edge-assembled-contig-file"}));
		}else{
			$opt{"edge-spades-pacbio-file"} = $input_dir."/".$opt{"edge-spades-pacbio-file"} if ($opt{"edge-spades-pacbio-file"} =~ /^\w/);
			&addMessage("PARAMS", "edge-spades-pacbio-file", "Invalid input. Fasta/q format required..") if ( -e $opt{"edge-spades-pacbio-file"} && ! is_fasta($opt{"edge-spades-pacbio-file"}) && ! is_fastq($opt{"edge-spades-pacbio-file"}));
			 $opt{"edge-spades-nanopore-file"} = $input_dir."/".$opt{"edge-spades-nanopore-file"} if ( $opt{"edge-spades-nanopore-file"} =~ /^\w/);
			&addMessage("PARAMS", "edge-spades-nanoport-file", "Invalid input. Fasta format required..") if ( -e $opt{"edge-spades-nanopore-file"} && ! is_fasta($opt{"edge-spades-nanopore-file"})  && ! is_fastq($opt{"edge-spades-nanopore-file"}));
		}
	}
	if ( $opt{"edge-anno-sw"} ){
		&addMessage("PARAMS", "edge-assembly-sw",  "You must turn on assembly function to do annotation.") unless $opt{"edge-assembly-sw"} >= 0; 
		#&addMessage("PARAMS", "edge-anno-kingdom",  "Invalid input. Natural number required.")     unless $opt{"edge-assembly-mink"} >= 0; 	
		if ($opt{"edge-anno-tool"} =~ /RATT/){
			 $opt{"edge-anno-source-file"}  = $input_dir."/".$opt{"edge-anno-source-file"} if ( $opt{"edge-anno-source-file"} =~ /^\w/);
			&addMessage("PARAMS","edge-anno-source-file","File not found. Please provide the Genbank Source File") if ( ! $opt{"edge-anno-source-file"} );
			&addMessage("PARAMS","edge-anno-source-file","Invalid input. Genbank format required") if ( -e $opt{"edge-anno-source-file"} && ! is_genbank($opt{"edge-anno-source-file"}) );
		}
	}
	if ( $opt{"edge-blast-sw"} ){
		&addMessage("PARAMS", "edge-blast-nt", "NT Database directory not found. Please check the path.") unless ( defined $opt{'edge-blast-nt'} && -d "$opt{'edge-blast-nt'}" ); 
		&addMessage("PARAMS", "edge-blast-nr", "NR Database directory not found. Please check the path.") unless ( defined $opt{'edge-blast-nr'} && -d "$opt{'edge-blast-nr'}" ); 
	}
	if ( $opt{"edge-phylo-sw"} ){
		&addMessage("PARAMS", "edge-phylo-patho", "Invalid input. Please select from precomputed SNP DB or form Genomes list.") if ( !$opt{'edge-phylo-patho'} && !defined $opt{'edge-phylo-ref-select'});
		&addMessage("PARAMS", "edge-phylo-ref-select", "Invalid input. Please select from precomputed SNP DB or form Genomes list.") if ( !$opt{'edge-phylo-patho'} && !defined $opt{'edge-phylo-ref-select'});
		&addMessage("PARAMS", "edge-phylo-patho", "You have both input types. Please select either from precomputed SNPdb OR form Genomes list.") if ( $opt{'edge-phylo-patho'} && ($opt{'edge-phylo-ref-select'} || $opt{'edge-phylo-ref-file'}));
		if ($opt{"edge-phylo-ref-file"}){
			@edge_phylo_ref_input = split /[\x0]/, $opt{"edge-phylo-ref-file"} if defined $opt{"edge-phylo-ref-file"};
			for my $i (0..$#edge_phylo_ref_input){
				$edge_phylo_ref_input[$i] = $input_dir."/".$edge_phylo_ref_input[$i] if ($edge_phylo_ref_input[$i]=~ /^\w/);
				my $id = 'edge-phylo-ref-file-'. ($i + 1);
				&addMessage("PARAMS",$id,"Invalid input. Genbank format required") if ( -e $edge_phylo_ref_input[$i] && ! is_fasta($edge_phylo_ref_input[$i]) );
			}
		}
	}
	if ($pipeline eq "qiime"){
		$opt{"edge-qiime-sw"} =1;
		$opt{"edge-qc-sw"} =0;
		$opt{"edge-hostrm-sw"} =0;
		$opt{"edge-assembly-sw"} = 0;
		$opt{"edge-ref-sw"} = 0;
		$opt{"edge-taxa-sw"} = 0;
		$opt{"edge-contig-taxa-sw"} = 0;
		$opt{"edge-anno-sw"} = 0 ;
		$opt{"edge-phylo-sw"} = 0 ;
		$opt{"edge-primer-valid-sw"} = 0 ;
		$opt{"edge-primer-adj-sw"} = 0 ;
		$opt{"edge-anno-sw"} = 0 ;
		$opt{"edge-jbroswe-sw"} = 0 ;
		@edge_qiime_barcode_input = split /[\x0]/, $opt{"edge-qiime-barcode-fq-file-input"} if defined $opt{"edge-qiime-barcode-fq-file-input"};
		foreach my $i (0..$#edge_qiime_barcode_input){
			my $id = "edge-qiime-barcode-fq-file-input". ($i + 1);
			$edge_qiime_barcode_input[$i] =~ s/ //g;
			$edge_qiime_barcode_input[$i] = "$input_dir/$edge_qiime_barcode_input[$i]" if ($edge_qiime_barcode_input[$i] =~ /^\w/);
			&addMessage("PARAMS","$id","Error: duplicated input.") if ($files{$edge_qiime_barcode_input[$i]});
			$files{$edge_qiime_barcode_input[$i]}=1;
			if ($edge_qiime_barcode_input[$i]  && -e $edge_qiime_barcode_input[$i]){
				&addMessage("PARAMS","$id","Input error. FASTQ format required") if ( ! is_fastq($edge_qiime_barcode_input[$i]));
			}else{
				&addMessage("PARAMS","$id","Input error. Please check the file path.");
				
			}
		}
		&addMessage("PARAMS", "edge-qiime-barcode-length","Invalid input. Natural number required.") unless $opt{"edge-qiime-barcode-length"}=~ /^\d+$/;
		&addMessage("PARAMS", "edge-qiime-phred-quality-threshold", "Invalid input. Input should in range 0-41.") unless ( $opt{"edge-qiime-phred-quality-threshold"} >= 0 && $opt{"edge-qiime-phred-quality-threshold"} <=41 );
		&addMessage("PARAMS", "edge-qiime-max-n","Invalid input. Natural number required.") unless $opt{"edge-qiime-max-n"}=~ /^\d+$/;
		&addMessage("PARAMS", "edge-qiime-min-per-read-length-fraction","Invalid input. Floating number between 0 and 1 required.") unless ( $opt{"edge-qiime-min-per-read-length-fraction"} >=0 && $opt{"edge-qiime-min-per-read-length-fraction"} <=1 );
		&addMessage("PARAMS", "edge-qiime-minimum-otu-size","Invalid input. Natural number required.") unless $opt{"edge-qiime-minimum-otu-size"}=~ /^\d+$/;
		&addMessage("PARAMS", "edge-qiime-similarity","Invalid input. Floating number between 0 and 1 required.") unless ( $opt{"edge-qiime-similarity"} >=0 && $opt{"edge-qiime-similarity"} <=1 );
		&addMessage("PARAMS", "edge-qiime-sampling-depth","Invalid input. Natural number required.") unless $opt{"edge-qiime-sampling-depth"}=~ /^\d+$/;

		
	}
}


sub parse_qiime_mapping_files{
	my $qiime_dir=shift;
	my $mapping_files=shift;
	my @pe1_files;
	my @pe2_files;
	my @se_files;
	if ( ! -d "$input_dir/$qiime_dir" ){
		my $msg = "ERROR: the input $qiime_dir directroy does not exist or isn't a directory\n";
		exit(1);
	}
	foreach my $f (@{$mapping_files}){
		my $file_column_index;
		open (my $fh, "$input_dir/$f") or die "Cannot read $f\n";
		while(<$fh>){
			chomp;
			next if (/^\n/);
			if (/SampleID/){
				my @header = split /\t/,$_;
				( $file_column_index )= grep { $header[$_] =~ /files/i } 0..$#header;
			}elsif(! /^#/){
				my @array = split /\t/,$_;
				my @files = map { "$qiime_dir/$_" } split /,|\s+/,$array[$file_column_index];
				if (scalar(@files) % 2){
					push @se_files,@files;
				}else{
					push @pe1_files,$files[0];
					push @pe2_files,$files[1];
				}
			}
		}
		close $fh;
	}
	return (\@pe1_files,\@pe2_files,\@se_files);
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

sub is_fasta
{
	$SIG{'PIPE'}=sub{};
	my $file=shift;
	my ($fh,$pid)= open_file($file);
	my $head=<$fh>;
	close $fh;
	kill 9, $pid; # avoid gunzip broken pipe
	$SIG{'PIPE'} = 'DEFAULT';
	my $return = ($head =~/^>/)? 1 : 0;
	return $return;
}

sub is_genbank
{
	$SIG{'PIPE'}=sub{};
	my $file=shift;
	my ($fh,$pid) = open_file($file);
	my $head=<$fh>;
	close $fh;
	kill 9, $pid; # avoid gunzip broken pipe
	$SIG{'PIPE'} = 'DEFAULT';
	my $return = ($head =~ /^LOCUS/i)? 1 : 0;
	return $return;
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
