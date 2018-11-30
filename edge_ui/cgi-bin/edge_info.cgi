#!/usr/bin/env perl
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-07
#

use strict;
use FindBin qw($RealBin);
use lib "../../lib";
use JSON;
use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser);
use POSIX qw(strftime);
use Data::Dumper;
use LWP::UserAgent;
use HTTP::Request::Common;
use Digest::MD5 qw(md5_hex);
use Email::Valid;
use Storable 'dclone';

require "./edge_user_session.cgi";
require "../cluster/clusterWrapper.pl";
$ENV{PATH} = "/bin:/usr/bin";

######################################################################################
# DATA STRUCTURE:
#
#     $info->{LIST}->{1}->{NAME}    // project name
#                       ->{TIME}    // submission time
#                       ->{PID}     // runpipeline pid
#                       ->{STATUS}  // status [Complete|running]
#                       ->{DESC}    // description
#                  ->{2}...
#
#     $info->{PROG}->{1}->{NAME}    // step name
#                       ->{DO}      // [0|1|auto]
#                       ->{STATUS}  // [unfinished|skip|already|running|done|failed]
#                  ->{2}...          
#
#     $info->{INFO}->{CPUU} ...
#
######################################################################################

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
my $init  = $opt{init};
$pname ||= $ARGV[0];
my $username    = $opt{'username'}|| $ARGV[1];
my $password    = $opt{'password'}|| $ARGV[2];
my $viewType    = "user";
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[3];
my $protocol = $opt{protocol} || 'http:';
my $sid         = $opt{'sid'}|| $ARGV[4];
my $ip          = $ARGV[5];
my $forceupdate = $opt{'forceupdate'};
$ENV{REMOTE_ADDR} = $ip if $ip;

my $domain      = $ENV{'HTTP_HOST'} || 'edge-dev-master.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;

my $list; # ref for project list
my @projlist; # project list index
my $prog; # progress for latest job
my $info; # info to return

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
$sys->{edgeui_input} = "$sys->{edgeui_input}"."/$webhostname" if ( -d "$sys->{edgeui_input}/$webhostname");
my $um_url      = $sys->{edge_user_management_url};
my $out_dir     = $sys->{edgeui_output};
my $www_root    = $sys->{edgeui_wwwroot};
my $edge_total_cpu = $sys->{"edgeui_tol_cpu"};
my $max_num_jobs = $sys->{"max_num_jobs"};
my $edge_projlist_num = $sys->{"edgeui_project_list_num"};
my $hideProjFromList = 0;
$um_url ||= "$protocol//$domain/userManagement";
$umSystemStatus ||= $sys->{user_management} if (!@ARGV);
$umSystemStatus = ($umSystemStatus eq "false")?0:$umSystemStatus;

#cluster
my $cluster 	= $sys->{cluster};
my $cluster_job_prefix = $sys->{cluster_job_prefix};
my $cluster_job_max_cpu= $sys->{cluster_job_max_cpu};
&LoadSGEenv($sys) if ($cluster);


my ($memUsage, $cpuUsage, $diskUsage) = &getSystemUsage();
$info->{INFO}->{CPUU} = $cpuUsage;
$info->{INFO}->{MEMU} = $memUsage;
$info->{INFO}->{DISKU} = $diskUsage;

my $runcpu = ($cluster)? int(($cluster_job_max_cpu-1)/$max_num_jobs): int(($edge_total_cpu-1)/$max_num_jobs);
$info->{INFO}->{RUNCPU} = ($runcpu>1)? $runcpu :1;

$info->{INFO}->{PROJLISTNUM} = $edge_projlist_num;


if ($init){
	&loadInitSetup();
	&returnStatus();;
}
#($umSystemStatus =~ /true/i)? &getUserProjFromDB():&scanNewProjToList();

#check projects vital
my ($vital, $name2pid, $error);
if($cluster) {
	($vital, $name2pid, $error) = checkProjVital_cluster($cluster_job_prefix);
	if($error) {
		$info->{INFO}->{ERROR}= "CLUSTER ERROR: $error";
	}
} else {
	($vital, $name2pid) = &checkProjVital();
}

# session check
if( $umSystemStatus ){
	my $valid = verifySession($sid);
	if($valid){
		$info->{INFO}->{SESSION_STATUS} = "valid";
		($username,$password,$viewType) = getCredentialsFromSession($sid);
		my $user_dir =  $sys->{edgeui_input}."/". md5_hex(lc($username));
		my $list_json = "$user_dir/.edgeinfo.json";
		my $user_config = $sys->{edgeui_input}."/". md5_hex(lc($username))."/user.properties";
		$info->{INFO}->{FORCE}=$forceupdate;
		if ( $forceupdate ne "true"  && -e $list_json) {
			$list = readListFromJson($list_json);
		}else{
			&getUserProjFromDB();
		}
		&getProjInfoFromDB($pname) if ($pname and ! defined $list->{$pname});
		&cache_user_projects_info($username,$password,$viewType,$list);	
	}
	else{
		&getUserProjFromDB();
		&getProjInfoFromDB($pname) if ($pname and ! defined $list->{$pname});
		$info->{INFO}->{SESSION_STATUS} = "invalid";
	}
}
else{
	&scanNewProjToList();
}

my $time = strftime "%F %X", localtime;
my $idx;
my @delete;
@projlist = sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
# selected project on index 0
if( scalar @projlist ){
	my $progs;
	my $count=0;
	# retrive progress info of a project that is selected by the following priorities:
	#  1. assigned project
	#  2. latest running project
	#  3. lastest project
	$idx= ($pname)? (grep $list->{$_}->{NAME} eq $pname, @projlist)[0] : $projlist[0];
	my @running_idxs = grep { $list->{$_}->{STATUS} eq "running" or $list->{$_}->{STATUS} =~ /unstarted|interrupted|in process|unknown/ and $list->{$_}->{NAME} ne $pname } @projlist;
	$idx = shift @running_idxs if (scalar(@running_idxs) && !$pname);
	$idx = $projlist[0] if (!$idx); # when given $pname does not exist.
	my @projlist_update = ($idx,@running_idxs); # update running projects and focus project program info.
	foreach my $i ( @projlist_update ) {
		last if ($edge_projlist_num && ++$count > $edge_projlist_num);
		my $lproj    = $list->{$i}->{NAME};
		my $lprojc   = $list->{$i}->{PROJCODE};
		my $lcpu     = $list->{$i}->{CPU};
		my $lstatus  = $list->{$i}->{STATUS};
		my $lpid     = $list->{$i}->{PID};
		my $realpid  = $name2pid->{$lproj}|| $name2pid->{$lprojc}; 
		my $proj_dir = "$out_dir/$lproj";
		$proj_dir = "$out_dir/$lprojc" if ( $lprojc && -d "$out_dir/$lprojc");
		my $log      = "$proj_dir/process.log";
		my $sjson    = "$proj_dir/.run.complete.status.json";
		my $current_log      = "$proj_dir/process_current.log";
		my $config   = "$proj_dir/config.txt";
		my $config_json = "$proj_dir/config.json";
		($list->{$i}->{PROJCONFIG} = $config_json) =~ s/$www_root// if ($i == $idx);
		#remove project from list if output directory has been removed	
		unless( -e $log || -e $config){
			delete $list->{$i};
			push @delete, $i;
			next;
		}	
	
		#status JSON
		#if( -e $sjson ){
		#	my $storedStatus = readListFromJson($sjson);
		#	$list->{$i} = $storedStatus->{LIST};
		#	$progs->{$i} = $storedStatus->{PROG};
		#	next;
		#}
		# update current project status
		if( -r $log ){
			my ($p_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id,$proj_runtime,$rnaPipeline) = &parseProcessLog($log);
			$list->{$i}->{TIME} ||= $proj_start;
			$list->{$i}->{TIME} ||= strftime "%F %X", localtime;
			$list->{$i}->{PID} = $realpid;
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{DESC} = $proj_desc;
			($list->{$i}->{PROJLOG} = $current_log) =~ s/$www_root//;
			#for unstarted project, read steps from config file
			$p_status = "unknown" if (!$p_status);
			(my $tmp,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id) = &parseProcessLog($config) if $p_status =~ /unstarted|unknown|interrupted/;
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{PROJNAME} = $proj_name;

			if( defined $name2pid->{$lproj} || defined $name2pid->{$lprojc} ){ #running
				$list->{$i}->{STATUS} = "running";
			}
			elsif( $p_status =~ /running/ ){
				unless ($rnaPipeline){
					# the process log reports it's running, but can't find vital
					# Unexpected exit detecteda
					$list->{$i}->{STATUS} = "failed";
					my $msg = "\n*** [$time] EDGE_UI: Pipeline failed (PID:$realpid $lproj $lprojc). Unexpected exit detected! ***";
					open (my $fh, ">>", "$proj_dir/process_current.log") or die "$!"; print $fh $msg; close $fh;
					open (my $lfh, ">>", "$log") or die "$!"; print $lfh $msg; close $lfh;
				}
			}
			else{
				$list->{$i}->{STATUS} = $p_status;
			}

			$progs->{$i} = $prog;
			if ( $list->{$i}->{DBSTATUS} && ($list->{$i}->{STATUS} ne $list->{$i}->{DBSTATUS})){
				&updateDBProjectStatus($i, $list->{$i}->{STATUS},$proj_start,$proj_runtime);
			}elsif( $p_status ne 'Complete'){
				&updateDBProjectStatus($i, $list->{$i}->{STATUS},$proj_start,$proj_runtime);
			}
			#if( $list->{$i}->{STATUS} eq "finished" && !-e $sjson ){
			#	my $storedStatus;
			#	$storedStatus->{LIST} = $list->{$i};
			#	$storedStatus->{PROG} = $progs->{$i};
			#	saveListToJason($storedStatus, $sjson);
			#}
		}
	}
	$info->{PROG} = $progs->{$idx};
	# with user management, NAME becomes unique project id
	$info->{INFO}->{NAME}   = $list->{$idx}->{NAME};
	$info->{INFO}->{PROJNAME}   = $list->{$idx}->{PROJNAME}; 
	$info->{INFO}->{PROJCODE}   = $list->{$idx}->{PROJCODE};
	$info->{INFO}->{PROJLOG} = $list->{$idx}->{PROJLOG};
	$info->{INFO}->{STATUS} = $list->{$idx}->{STATUS};
	$info->{INFO}->{TIME}   = strftime "%F %X", localtime;
	$info->{INFO}->{PROJTYPE} = $list->{$idx}->{PROJTYPE} if ($list->{$idx}->{PROJTYPE});
	$info->{INFO}->{PROJCONFIG} = $list->{$idx}->{PROJCONFIG} if ($list->{$idx}->{PROJCONFIG});

	## sample metadata
	$info->{INFO}->{SHOWMETA}   = $list->{$idx}->{SHOWMETA} if ($list->{$idx}->{SHOWMETA});
	$info->{INFO}->{ISOWNER} = $list->{$idx}->{ISOWNER} if($list->{$idx}->{ISOWNER});
	$info->{INFO}->{HASMETA}   = $list->{$idx}->{HASMETA} if ($list->{$idx}->{HASMETA});
	$info->{INFO}->{SHAREBSVE}   = $list->{$idx}->{SHAREBSVE} if ($list->{$idx}->{SHAREBSVE});
	$info->{INFO}->{METABSVE}   = $list->{$idx}->{METABSVE} if ($list->{$idx}->{METABSVE});
	## END sample metadata
}
delete $list->{$_} foreach @delete;
$info->{LIST} = $list if $list;
&returnStatus();

######################################################

sub cache_user_projects_info {
	my $user = shift;
	my $pass = shift;
	my $type= shift;
	my $list = shift;
	
	my $user_dir =  $sys->{edgeui_input}."/". md5_hex(lc($user));
	my $upage_json = "$user_dir/.user_pp.json";
	my $apage_json = "$user_dir/.admin_pp.json";
	my $list_json = "$user_dir/.edgeinfo.json";
	my $now = time(); 
	# fork this process
	my $pid = fork();
	die "Fork failed: $!" if !defined $pid;
	if ($pid==0){
		# releases browser from waiting child to finish
		open STDERR, ">/dev/null";
		open STDIN, "</dev/null";
		open STDOUT, ">/dev/null";
		if ( ! -e $list_json || ($now-(stat $list_json)[9]) > 180 || $forceupdate eq "true") {
			if (!$list){
				&getUserProjFromDB();
				&getProjInfoFromDB($pname) if ($pname and ! defined $list->{$pname});
			}
			my $list_from_api = dclone($list);
			saveListToJason($list_from_api, $list_json);
		}
		# update if file age > 5 min = 300 sec
		#if ( ! -e $upage_json || ($now-(stat $upage_json)[9]) > 300 || $forceupdate eq "true") {
		#	if ($type eq "admin"){
		#		system("perl edge_projectspage.cgi $user $pass true $type user \'\' true 2>/dev/null 1>/dev/null");
		#		system("perl edge_projectspage.cgi $user $pass true $type $type \'\' true  2>/dev/null 1>/dev/null");
		#	}else{
		#		system("perl edge_projectspage.cgi $user $pass true $type user \'\' true  2>/dev/null 1>/dev/null");
		#	}
		#}
		exit;
	}
}

sub readListFromJson {
	my $json = shift;
	my $list = {};
	if( -r $json ){
		open (JSON, "<", $json) or die "$!";
		flock(JSON, 1);
  		local $/ = undef;
  		$list = decode_json(<JSON>);
  		close JSON;
	}
	return $list;
}

sub saveListToJason {
	my ($list, $file) = @_;
	open JSON, ">", "$file" or die "Can't write to file: $file\n";
  	my $json = encode_json($list);
	print JSON $json;
  	close JSON;
}



sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open CONF, "<", $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			$flag=1;
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
	die "Incorrect system file\n" if (!$flag);
	return $sys;
}

sub parseProcessLog {
	my ($log)=@_;
	my $cnt=0;
	my $prog;
	my $lastline;
	my $proj_status="unknown";
	my $proj_start;
	my $proj_desc;
	my $proj_name;
	my $proj_id;
	my $numcpu;
	my $tol_running_sec=0;
	my $last_runtime;
	my ($step,$ord,$do,$status);
	my %map;
	my $rnaPipeline=0;

	open LOG, "<", $log or die "Can't open $log.";
	foreach(<LOG>){
		chomp;
		next if /^$/;
		next if /^#/;
		if( /Total Running time: (\d+):(\d+):(\d+)/){
			$last_runtime = "$1h $2m $3s";
			next;
		}
		if( /Project Start: (\d{4}) (\w{3})\s+(\d+)\s+(.*)/){
			my ($yyyy,$mm,$dd,$hms) = ($1,$2,$3,$4);
			my %mon2num = qw(jan 1  feb 2  mar 3  apr 4  may 5  jun 6  jul 7  aug 8  sep 9  oct 10 nov 11 dec 12);
			$mm = $mon2num{ lc substr($mm, 0, 3) };
			$mm = sprintf "%02d", $mm;
			$dd = sprintf "%02d", $dd;
			$proj_start  = "$yyyy-$mm-$dd $hms";
			$proj_status = "unknown";
			$cnt=0;
			$numcpu=0;
			undef %{$prog};
			undef %map;
		}
		elsif(/runPipeline_rRNA/){
			$rnaPipeline=1;
		}
		elsif( /^cpu=(\d+)$/ ){
			$numcpu=$1;
		}
	
		elsif( /^projdesc=(.*)/ ){
			$proj_desc=$1;
		}
		elsif( /^projname=(.*)/){
			$proj_name=$1;
		}
		elsif( /^projid=(.*)/){
			$proj_id=$1;
		}
		elsif( /^\[(.*)\]/ ){
			my $step = $1;
			next if $step eq "system" or $step eq "project";

			if( defined $map{"$step"} ){
				$ord = $map{"$step"};
			}
			else{
				$cnt++;
				my $step = $1;
				$prog->{$cnt}->{NAME}=$step;
				$map{"$step"}=$cnt;
			}
		}
		elsif( /^Do.*=(.*)$/ ){
			next if (/KeggOmicsView/);
			my $do = $1;
			$prog->{$cnt}->{DO}= 'auto' if ($do eq 'auto');
			$prog->{$cnt}->{DO}= 1 if ($do eq 1);
			$prog->{$cnt}->{DO}= 0 if ($do eq 0 && !$prog->{$cnt}->{DO});
			$prog->{$cnt}->{STATUS}="skip";
			$prog->{$cnt}->{STATUS}="unfinished" if ($prog->{$cnt}->{DO});
		}
		elsif( /Finished/ ){
			$prog->{$ord}->{STATUS} = "finished";
		}
		elsif( /Running time: (\d+:\d+:\d+)/ ){
			$prog->{$ord}->{STATUS} = "done";
			$prog->{$ord}->{TIME} = $1;
			my ($h,$m,$s) = $1 =~ /(\d+):(\d+):(\d+)/;
			$tol_running_sec += $h*3600+$m*60+$s;
		}
		elsif( /Running/ ){
			$prog->{$ord}->{STATUS} = "running";
			$proj_status="running";
		}
		elsif( /failed/ ){
			$prog->{$ord}->{STATUS} = "failed";
			$proj_status="failed";
		}
		elsif( /^Cannot/ ){
			$prog->{$ord}->{STATUS} = "failed";
			$proj_status="failed";
		}
		elsif( /^All Done\./ ){
			$prog->{$ord}->{STATUS} = "done";
			$proj_status="Complete";
		}
		$lastline = $_;
	}
	close LOG;

	#unstarted project
	$proj_status            = "unstarted"   if $lastline =~ /EDGE_UI.*unstarted/;
	$proj_status            = "interrupted" if $lastline =~ /EDGE_UI.*interrupted/;
	$proj_status            = "archived" if $lastline =~ /EDGE_UI.*archived/;
	$proj_start             = $1            if $lastline =~ /\[(\S+ \S+)\] EDGE_UI/;
	$prog->{$ord}->{STATUS} = "unfinished"  if $proj_status eq "interrupted"; #turn last step to unfinished
	
	my $proj_runtime = sprintf("%02d:%02d:%02d", int($tol_running_sec / 3600), int(($tol_running_sec % 3600) / 60), int($tol_running_sec % 60));
	#change unfinished auto steps to "skip"
	my $flag=0;
	foreach my $ord ( sort {$b<=>$a} keys %$prog ){
		if( $prog->{$ord}->{STATUS} =~ /(finished|done|running|failed)/ ){
			$flag=1;
		}
		if( $flag && $prog->{$ord}->{DO} eq "auto" && $prog->{$ord}->{STATUS} eq "unfinished" ){
			$prog->{$ord}->{STATUS}="skip";
		}
	}

	return ($proj_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id,$proj_runtime,$rnaPipeline);
}

sub scanNewProjToList {
	my $cnt = 1;
	
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	my @dirfiles = readdir(BIN);
	foreach my $file (@dirfiles)  {
		next if ($file eq '.' || $file eq '..' || ! -d "$out_dir/$file");
		my $config = "$out_dir/$file/config.txt";
		my $processLog = "$out_dir/$file/process_current.log";
		$cnt++;
		if (-r "$config"){
			$list->{$cnt}->{NAME} = $file ;
			$list->{$cnt}->{TIME} =  strftime "%F %X",localtime((stat("$config"))[9]); 
			$list->{$cnt}->{STATUS} = "unknown";
			if ( -r "$processLog"){
				open (my $fh, "<",$processLog) or die "$!";
				while(<$fh>){
					if (/queued/){
						$list->{$cnt}->{STATUS} = "unstarted";
					}
					if (/All Done/){
						$list->{$cnt}->{STATUS} = "Complete";
					}
					if (/failed/i){
						$list->{$cnt}->{STATUS} = "failed";
					}
				}
				close $fh;
			}
			$list->{$cnt}->{STATUS} = "running" if $name2pid->{$file};
			my $projname = $file;
			# if the system change from User management on to off. will need parse the project name from config file
			#  using grep is slow than open file and regrex
		#unlink glob("$edge_input/$user_dir/.sid*");
			#if ( length($projname)==32 ){
			#	open (my $fh, $config);
			#	while(<$fh>){if (/projname=(.*)/){$projname=$1;chomp $projname;last;};}
			#	close $fh;
			#}
			chomp $projname;
			$list->{$cnt}->{PROJNAME} = $projname;

			## sample metadata
			$list->{$cnt}->{ISOWNER} = 1;
			if($sys->{edge_sample_metadata}) {
				$list->{$cnt}->{SHOWMETA} = 1;
			}
			my $metaFile = "$out_dir/$file/metadata_sample.txt";
			my $runFile = "$out_dir/$file/metadata_run.txt";
			my $pathogensFile = "$out_dir/$file/pathogens.txt";
			if(-r $metaFile) {
				$list->{$cnt}->{HASMETA} = 1;
			} 

			if($sys->{edge_sample_metadata_share2bsve} && hasPathogens($pathogensFile)) {
				$list->{$cnt}->{SHAREBSVE} = 1;
			}
			if(-r $runFile) {
				my $bsveId;
				open (my $fh, "<", $runFile) or die "$!";
				while(<$fh>){
					chomp;
					if (/bsve_id=(.*)/){ $bsveId = $1; last;}
				}
				close $fh;
				chomp $bsveId;
				$list->{$cnt}->{METABSVE} = $bsveId;
			}
			## END sample metadata
		}
	}
	closedir(BIN);
}

sub availableToRun {
	my ($num_cpu, $cpu_been_used) = @_;
	return 0 if $cpu_been_used + $num_cpu >= $sys->{edgeui_tol_cpu};

	if( $sys->{edgeui_auto_queue} && $sys->{edgeui_tol_cpu} ){
		foreach my $i ( keys %$list ){
			if( $list->{$i}->{STATUS} eq "running" ){
				$cpu_been_used += $list->{$i}->{CPU};
				return 0 if $cpu_been_used + $num_cpu >= $sys->{edgeui_tol_cpu};
			}
		}
		return 1;
	}
	else{
		return 0;
	}
}

sub getSystemUsage {
	my $mem = `vmstat -s | awk  '\$0 ~/total memory/ {total=\$1 } \$0 ~/free memory/ {free=\$1} \$0 ~/buffer memory/ {buffer=\$1} \$0 ~/cache/ {cache=\$1} END{print (total-free-buffer-cache)/total*100}'`;
	my $cpu = `top -bn1 | grep load | awk '{printf "%.1f", \$(NF-2)}'`;
	my $disk;
	open (my $fh, "-|")
	  or exec ("df","-h","$out_dir");
	while(<$fh>){ my @array= split/\s+/; $disk=($array[4] =~ /\%/)?$array[4]:$array[3];}
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

sub checkProjVital {
	my $ps = `ps auxww | grep run[P]`;
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

sub addInfo {
	my ($cate, $type, $note) = @_;
	$info->{$cate}->{$type}=$note;
}


sub updateDBProjectStatus{
	my $project = shift;
	my $status = shift;
	my $submitted = shift;
	my $running_time = shift;
	my %data = (
                email => $username,
                password => $password,
		project_id => $project,
		new_project_status => $status,
		run_submitted => $submitted,
        );
	if($running_time) {
		$data{running_time} = $running_time;
	}
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
        my $result =  decode_json($result_json);
	if (! $result->{status})
	{
		$info->{INFO}->{ERROR}="Update Project status in database failed.";
	}
}

sub getUserProjFromDB{
    my %data = (
            email => $username,
            password => $password
    );
    # Encode the data structure to JSON
    #w Set the request parameters
	my $service;
	if ($username && $password){ 
		$service = "WS/user/getProjects";
		#$data{project_display} = "yes";
	}else{
		$service = "WS/user/publishedProjects";
	}

    # Encode the data structure to JSON
	my $data = to_json(\%data);
    my $url = $um_url .$service;
    #w Set the request parameters
	my $browser = LWP::UserAgent->new;
	my $req = PUT $url;
    $req->header('Content-Type' => 'application/json');
    $req->header('Accept' => 'application/json');
    #must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
    $req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	#print $result_json,"\n" if ($ARGV[1]);
	if ($result_json =~ /\"error_msg\":"(.*)"/)
	{
		$info->{INFO}->{ERROR}=$1;
		return;
	}
	my $array_ref =  decode_json($result_json);
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $projCode = $hash_ref->{code};
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		my $created = $hash_ref->{created};
		my $proj_dir = (-d "$out_dir/$projCode")?"$out_dir/$projCode":"$out_dir/$id";
		next if (! -r "$proj_dir/process.log" && !$cluster);
		next if ($cluster && ! -r "$proj_dir/clusterSubmit.sh");
		next if ( $status =~ /delete/i);
		$list->{$id}->{NAME} = $id;
		$list->{$id}->{PROJNAME} = $project_name;
		$list->{$id}->{PROJCODE} = $projCode;
		$list->{$id}->{DBSTATUS} = $status;
		$list->{$id}->{STATUS} = $status;
		$list->{$id}->{TIME} = $created;
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{OWNER_FisrtN} = $hash_ref->{owner_firstname};
		$list->{$id}->{OWNER_LastN} = $hash_ref->{owner_lastname};
		$list->{$id}->{PROJTYPE} = $hash_ref->{type} if ($username && $password);

		## sample metadata
		if($sys->{edge_sample_metadata}) {
			$list->{$id}->{SHOWMETA} = 1;
		}
		if($username eq  $hash_ref->{owner_email}) {
			$list->{$id}->{ISOWNER} = 1;
		}
		my $metaFile = "$proj_dir/metadata_sample.txt";
		my $runFile = "$proj_dir/metadata_run.txt";
		my $pathogensFile = "$proj_dir/pathogens.txt";
		if(-r $metaFile) {
			$list->{$id}->{HASMETA} = 1;
		}
		if($sys->{edge_sample_metadata_share2bsve} && hasPathogens($pathogensFile)) {
			$list->{$id}->{SHAREBSVE} = 1;
		}
		if(-r $runFile) {
			my $bsveId; 
			open (my $fh, "<", $runFile) or die "$!";
			while(<$fh>){
				chomp;
				if (/bsve_id=(.*)/){ $bsveId = $1; last;}
			}
			close $fh;
			chomp $bsveId;
			$list->{$id}->{METABSVE} = $bsveId;
		} 
		## END sample metadata
		
	}
}

sub hasPathogens {
	my $file = shift;
	my $top = 0;

	if(-e $file) {
	    	open (my $fh , "<", $file) or die "No config file $!\n";
	   	 while (<$fh>) {
	       	 	chomp;
			if(/pathogen\thost\(s\)\tdisease/) {
				next;
			}
			$top ++;
	    	}
	    	close $fh;
	}

 	return $top;
}

sub getProjInfoFromDB{
    my $project=shift;
    $project = &getProjID($project);
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
	my $hash_ref = decode_json($result_json);

	my $id = $hash_ref->{id};
	my $project_name = $hash_ref->{name};
	my $projCode = $hash_ref->{code};
	my $status = $hash_ref->{status};
	my $created = $hash_ref->{created};
	my $projtype = ($hash_ref->{isPublished})?"publish":"false";
	my $proj_dir = (-d "$out_dir/$projCode")?"$out_dir/$projCode":"$out_dir/$id";
        return if (! -r "$proj_dir/process.log");
        return if ( $status =~ /delete/i);
	$list->{$id}->{NAME} = $id;
	$list->{$id}->{PROJNAME} = $project_name;
	$list->{$id}->{PROJCODE} = $projCode;
	$list->{$id}->{DBSTATUS} = $status;
	$list->{$id}->{STATUS} = $status;
	$list->{$id}->{PROJTYPE} = $projtype;
	$list->{$id}->{TIME} = $created;
}

sub getProjID {
  my $project=shift;
  my $projID = $project;
  if ( -d "$out_dir/$project"){ # use ProjCode as dir
    open (my $fh, "<", "$out_dir/$project/config.txt") or die "Cannot open $out_dir/$project/config.txt\n";
    while(<$fh>){
      if (/^projid=(\S+)/){
        $projID = $1;
        last;
      }
    }
  }
  return $projID;
}

sub loadInitSetup{
	# module on/off
	$info->{INFO}->{UMSYSTEM}= ( $sys->{user_management} )? "true":"false";
	$info->{INFO}->{UPLOAD}  = ( $sys->{user_upload} )?"true":"false";
	$info->{INFO}->{ARCHIVE} = ( -w $sys->{edgeui_archive} ) ? "true":"false";
	$info->{INFO}->{MQC}     = ( $sys->{m_qc} )?"true":"false";
	$info->{INFO}->{MAA}     = ( $sys->{m_assembly_annotation} )?"true":"false";
	$info->{INFO}->{MRBA}    = ( $sys->{m_reference_based_analysis} )?"true":"false";
	$info->{INFO}->{MTC}     = ( $sys->{m_taxonomy_classfication} )?"true":"false";
	$info->{INFO}->{MPA}     = ( $sys->{m_phylogenetic_analysis} )?"true":"false";
	$info->{INFO}->{MSGP}    = ( $sys->{m_specialty_genes_profiling} )?"true":"false";
	$info->{INFO}->{MPPA}    = ( $sys->{m_pcr_primer_analysis} )?"true":"false";
	$info->{INFO}->{MQIIME}  = ( $sys->{m_qiime} )?"true":"false";
	$info->{INFO}->{MTARGETEDNGS}  = ( $sys->{m_targetedngs} )?"true":"false";
	$info->{INFO}->{MPIRET}  = ( $sys->{m_piret} )?"true":"false";
	#parameters
	$info->{INFO}->{UPLOADEXPIRE}  = ( $sys->{edgeui_proj_store_days} )?$sys->{edgeui_proj_store_days}:"1095";
	$info->{INFO}->{UPLOADFILEEXT}  = ( $sys->{user_upload_fileext} )?$sys->{user_upload_fileext}:"fastq,fq,fa,fasta,fna,contigs,gbk,gbff,genbank,gb,txt,text,config,ini,xls,xlsx";

}

sub returnStatus {
	my $json;
	$json = to_json( $info ) if $info;
	$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[1];
	print $cgi->header('application/json'), $json;
	exit;
}
sub stringSanitization{
	my $opt=shift;
	my $dirtybit=0;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		next if $key eq "password";
		next if $key eq "keywords";

		if ($key eq "username"){
			$dirtybit =1  if ! &emailValidate($str);
		}else{
			$opt->{$key} =~ s/[`;'"]/ /g;
			$dirtybit=1 if ($opt->{$key} =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/);
		}
		if ($dirtybit){
			$info->{INFO} = "Invalid characters detected \'$str\'.";
			&returnStatus();
		}
	}
}

sub emailValidate{
	my $email=shift;
	$email = lc($email);
	return Email::Valid->address($email);
}


