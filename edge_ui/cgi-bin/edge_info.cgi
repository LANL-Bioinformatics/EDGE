#!/usr/bin/env perl
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-07
#

use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser);
use POSIX qw(strftime);
use Data::Dumper;
use LWP::UserAgent;
use HTTP::Request::Common;
use Digest::MD5 qw(md5_hex);
require "edge_user_session.cgi";
require "../cluster/clusterWrapper.pl";

######################################################################################
# DATA STRUCTURE:
#
#     $info->{LIST}->{1}->{NAME}    // project name
#                       ->{TIME}    // submission time
#                       ->{PID}     // runpipeline pid
#                       ->{STATUS}  // status [finished|running]
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
$ENV{REMOTE_ADDR} = $ip if $ip;

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $um_url      = $sys->{edge_user_management_url};
my $out_dir     = $sys->{edgeui_output};
my $www_root    = $sys->{edgeui_wwwroot};
my $edge_total_cpu = $sys->{"edgeui_tol_cpu"};
my $edge_projlist_num = $sys->{"edgeui_project_list_num"};
my $domain      = $ENV{'HTTP_HOST'};
my $hideProjFromList = 0;
$domain ||= "edgeset.lanl.gov";
$um_url ||= "$protocol//$domain/userManagement";

#cluster
my $cluster 	= $sys->{cluster};
my $cluster_job_prefix = $sys->{cluster_job_prefix};
my $cluster_job_max_cpu= $sys->{cluster_job_max_cpu};
my $list; # ref for project list
my @projlist; # project list index
my $prog; # progress for latest job
my $info; # info to return

my ($memUsage, $cpuUsage, $diskUsage) = &getSystemUsage();
$info->{INFO}->{CPUU} = $cpuUsage;
$info->{INFO}->{MEMU} = $memUsage;
$info->{INFO}->{DISKU} = $diskUsage;

my $runcpu = ($cluster)? int($cluster_job_max_cpu/8): int($edge_total_cpu/8);
$info->{INFO}->{RUNCPU} = ($runcpu>1)? $runcpu :1;

$info->{INFO}->{PROJLISTNUM} = $edge_projlist_num;

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

&returnStatus() if ($init);
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
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid){
		($username,$password,$viewType) = getCredentialsFromSession($sid);
		my $user_config = $sys->{edgeui_input}."/". md5_hex($username)."/user.properties";
		&getUserProjFromDB();
		&getProjInfoFromDB($pname) if ($pname and ! defined $list->{$pname});
		$info->{INFO}->{SESSION_STATUS} = "valid";
	}
	else{
		&getUserProjFromDB();
		$info->{INFO}->{SESSION_STATUS} = "invalid";
	}
}
else{
	&scanNewProjToList();
}

my $time = strftime "%F %X", localtime;

@projlist = sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
# selected project on index 0
@projlist=((grep $list->{$_}->{NAME} eq $pname, @projlist), (grep $list->{$_}->{NAME} ne $pname, @projlist)) if ($pname);
if( scalar @projlist ){
	my $idx;
	my $progs;
	my $count=0;
	
	# retrive progress info of a project that is selected by the following priorities:
	#  1. assigned project
	#  2. latest running project
	#  3. lastest project
	my @running_idxs = grep { $list->{$_}->{STATUS} eq "running" or $list->{$_}->{NAME} eq $pname or $list->{$_}->{STATUS} eq "unstarted" or $list->{$_}->{STATUS} eq "in process" } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @not_running_idxs;
	$idx = $projlist[0];
	if(scalar(@running_idxs)){
		$idx = $running_idxs[0] if (!$pname);
		@not_running_idxs = grep { $list->{$_}->{STATUS} ne "running" and $list->{$_}->{NAME} ne $pname and $list->{$_}->{STATUS} ne "unstarted" and $list->{$_}->{STATUS} ne "in process"} sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
		@projlist = (@running_idxs,@not_running_idxs);
	}
	
	foreach my $i ( @projlist ) {
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
		my $current_log      = "$proj_dir/process_current.log";
		my $config   = "$proj_dir/config.txt";

		#remove project from list if output directory has been removed
		unless( -e $log ){
			delete $list->{$i};
			next;
		}	
		
		# update current project status
		if( -r $log ){
			my ($p_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id) = &parseProcessLog($log);
			$list->{$i}->{TIME} = $proj_start;
			$list->{$i}->{TIME} ||= strftime "%F %X", localtime;
			$list->{$i}->{PID} = $realpid;
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{DESC} = $proj_desc;
			($list->{$i}->{PROJLOG} = $current_log) =~ s/$www_root//;

			#for unstarted project, read steps from config file
			(my $tmp,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id) = &parseProcessLog($config) if $p_status eq "unstarted";
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{PROJNAME} = $proj_name;

			if( defined $name2pid->{$lproj} || defined $name2pid->{$lprojc} ){ #running
				$list->{$i}->{STATUS} = "running";
			}
			elsif( $p_status =~ /running/ ){
				# the process log reports it's running, but can't find vital
				# Unexpected exit detected
				$list->{$i}->{STATUS} = "failed";
				`echo "\n*** [$time] EDGE_UI: Pipeline failed (PID:$realpid). Unexpected exit detected! ***" |tee -a $log >> $proj_dir/process_current.log`;
			}
			else{
				$list->{$i}->{STATUS} = $p_status;
			}

			$progs->{$i} = $prog;
			if ( $list->{$i}->{DBSTATUS} && ($list->{$i}->{STATUS} ne $list->{$i}->{DBSTATUS})){
				&updateDBProjectStatus($i, $list->{$i}->{STATUS});
			}
		}
	}

	$info->{PROG} = $progs->{$idx};
	# with user management, NAME becomes unique project id
	$info->{INFO}->{NAME}   = $list->{$idx}->{NAME};
	$info->{INFO}->{PROJNAME}   = $list->{$idx}->{PROJNAME}; 
	$info->{INFO}->{PROJCODE}   = $list->{$idx}->{PROJCODE};
	$info->{INFO}->{PROJLOG} = $list->{$idx}->{PROJLOG};;
	$info->{INFO}->{STATUS} = $list->{$idx}->{STATUS};
	$info->{INFO}->{TIME}   = strftime "%F %X", localtime;
	$info->{INFO}->{PROJTYPE} = $list->{$idx}->{PROJTYPE} if ($list->{$idx}->{PROJTYPE});
}

#autorun
if( scalar keys %$list && $sys->{edgeui_auto_run} ){
	my ( $progs, $proj, $projCode, $p_status, $proj_start, $proj_dir, $log, $config );
	my $num_cpu_used = 0;
	foreach my $i ( sort {$list->{$a}->{TIME} cmp $list->{$b}->{TIME}} keys %$list ) {
		$proj     = $list->{$i}->{NAME};
		$projCode = $list->{$i}->{PROJCODE};
		$proj_dir = "$out_dir/$projCode";
		my $run=0;
		$run = &availableToRun($list->{$i}->{CPU}, $num_cpu_used ) if $list->{$i}->{STATUS} eq "unstarted";
		if($run && ($umSystemStatus eq 'ture'?$sid:1) ){
			my $json = `$RealBin/edge_action.cgi $proj rerun "" "" "" $sid $domain 2>> $proj_dir/error.log`;
			#print STDERR "$json";
			my $info = decode_json($json);
			$list->{$i}->{STATUS} = "running" if $info->{STATUS} eq "SUCCESS";
			$num_cpu_used += $list->{$i}->{CPU};
		}
	}
}

$info->{LIST} = $list if $list;

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
	my ($step,$ord,$do,$status);
	my %map;

	open LOG, $log or die "Can't open $log.";
	foreach(<LOG>){
		chomp;
		next if /^$/;
		next if /^#/;
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
			$proj_status="finished";
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

	return ($proj_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name,$proj_id);
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
			$list->{$cnt}->{TIME} = (stat("$out_dir/$file"))[10]; 
			$list->{$cnt}->{STATUS} eq "running" if $name2pid->{$file};
			if ( -r "$processLog"){
				open (my $fh, $processLog);
				while(<$fh>){
					if (/queued/){
						$list->{$cnt}->{STATUS} eq "unstarted";
						last;
					}
				}
				close $fh;
			}
			#my $projname = `grep -a "projname=" $out_dir/$file/config.txt | awk -F'=' '{print \$2}'`;
			my $projname = $file;
			chomp $projname;
			$list->{$cnt}->{PROJNAME} = $projname;
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
	my $mem = `free -m | awk 'NR==2{printf "%.1f", \$3*100/(\$4+\$3)}'`;
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

sub addInfo {
	my ($cate, $type, $note) = @_;
	$info->{$cate}->{$type}=$note;
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
	#print $result_json,"\n";
	if ($result_json =~ /\"error_msg\":"(.*)"/)
	{
		$info->{INFO}->{ERROR}=$1;
		return;
	}
	my $array_ref =  from_json($result_json);
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $projCode = $hash_ref->{code};
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		my $created = $hash_ref->{created};
		next if (! -r "$out_dir/$id/process.log" && ! -r "$out_dir/$projCode/process.log");
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
	}
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
	print $result_json if (@ARGV);
	my $hash_ref = from_json($result_json);

	my $id = $hash_ref->{id};
	my $project_name = $hash_ref->{name};
	my $projCode = $hash_ref->{code};
	my $status = $hash_ref->{status};
	my $created = $hash_ref->{created};
	my $projtype = ($hash_ref->{isPublished})?"publish":"false";
	#next if (! -r "$out_dir/$id/process.log");
	#next if ( $status =~ /delete/i);
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
    open (my $fh, "$out_dir/$project/config.txt") or die "Cannot open $out_dir/$project/config.txt\n";
    while(<$fh>){
      if (/^projid=(\S+)/){
        $projID = $1;
        last;
      }
    }
  }
  return $projID;
}

sub returnStatus {
	my $json;
	$json = to_json( $info ) if $info;
	$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[1];
	print $cgi->header('application/json'), $json;
	exit;
}

