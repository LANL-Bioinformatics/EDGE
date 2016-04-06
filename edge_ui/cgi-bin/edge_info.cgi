#!/usr/bin/env perl
#
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
require "edge_user_session.cgi";

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
######################################################################################

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
$pname ||= $ARGV[0];
my $username    = $opt{'username'}|| $ARGV[1];
my $password    = $opt{'password'}|| $ARGV[2];
my $viewType    = "user";
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[3];
my $protocol = $opt{protocol} || 'http:';
my $sid         = $opt{'sid'}|| $ARGV[4];
my $ip          = $ARGV[5];
$ENV{REMOTE_ADDR} = $ip if $ip;

# read system params from config template
my $config_tmpl = "$RealBin/edge_config.tmpl";
my $sys         = &getSysParamFromConfig($config_tmpl);
my $um_url      = $sys->{edge_user_management_url};
my $out_dir     = $sys->{edgeui_output};
my $domain      = $ENV{'HTTP_HOST'};
my $hideProjFromList = 0;
$domain ||= "edgeset.lanl.gov";
$um_url ||= "$protocol//$domain/userManagement";

my $list; # ref for project list
my $prog; # progress for latest job
my $info; # info to return

my ($memUsage, $cpuUsage, $diskUsage) = &getSystemUsage();
$info->{INFO}->{CPUU} = $cpuUsage;
$info->{INFO}->{MEMU} = $memUsage;
$info->{INFO}->{DISKU} = $diskUsage;

$info->{INFO}->{UPLOAD} = "true" if  ( $sys->{user_upload} );

#($umSystemStatus =~ /true/i)? &getUserProjFromDB():&scanNewProjToList();

# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid){
		($username,$password,$viewType) = getCredentialsFromSession($sid);
		&getUserProjFromDB();
		&getProjInfoFromDB($pname) if ! defined $list->{$pname};
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


#check projects vital
my ($vital, $name2pid) = &checkProjVital();

my $time = strftime "%F %X", localtime;

if( scalar keys %$list ){
	my $idx;
	my $progs;

	foreach my $i ( keys %$list ) {
		my $lproj    = $list->{$i}->{NAME};
		my $lcpu     = $list->{$i}->{CPU};
		my $lstatus  = $list->{$i}->{STATUS};
		my $lpid     = $list->{$i}->{PID};
		my $realpid  = $name2pid->{$lproj}; 
		my $proj_dir = "$out_dir/$lproj";
		my $log      = "$proj_dir/process.log";
		my $config   = "$proj_dir/config.txt";
		$idx         = $i if $lproj eq $pname;

		#remove project from list if output directory has been removed
		unless( -e $log ){
			delete $list->{$i};
			next;
		}	
		
		# update current project status
		if( -r $log ){
			my ($p_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name) = &parseProcessLog($log);
			$list->{$i}->{TIME} = $proj_start;
			$list->{$i}->{TIME} ||= strftime "%F %X", localtime;
			$list->{$i}->{PID} = $realpid;
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{DESC} = $proj_desc;

			#for unstarted project, read steps from config file
			(my $tmp,$prog,$proj_start,$numcpu,$proj_desc,$proj_name) = &parseProcessLog($config) if $p_status eq "unstarted";
			$list->{$i}->{CPU} = $numcpu;
			$list->{$i}->{PROJNAME} = $proj_name;

			if( defined $name2pid->{$lproj} ){ #running
				$list->{$i}->{STATUS} = "running";
			}
			elsif( $p_status =~ /running/ ){
				# the process log reports it's running, but can't find vital
				# Unexpected exit detected
				$list->{$i}->{STATUS} = "failed";
				`echo "\n*** [$time] EDGE_UI: Pipeline failed (PID:$lpid). Unexpected exit detected! ***" |tee -a $log >> $proj_dir/process_current.log`;
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

	# retrive progress info of a project that is selected by the following priorities:
	#  1. assigned project
	#  2. latest running project
	#  3. lastest project

	unless( $idx ){
		my @idxs1 = grep { $list->{$_}->{STATUS} eq "running" } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
		my @idxs2 = sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
		my @idxs = (@idxs1,@idxs2);
		$idx = shift @idxs;
	}
	$info->{PROG} = $progs->{$idx};
	# with user management, NAME becomes unique project id
	$info->{INFO}->{NAME}   = $list->{$idx}->{NAME};
	$info->{INFO}->{PROJNAME}   = $list->{$idx}->{PROJNAME}; 
	$info->{INFO}->{STATUS} = $list->{$idx}->{STATUS};
	$info->{INFO}->{TIME}   = strftime "%F %X", localtime;
	$info->{INFO}->{PROJTYPE} = $list->{$idx}->{PROJTYPE} if ($list->{$idx}->{PROJTYPE});
}

#autorun
if( scalar keys %$list && $sys->{edgeui_auto_run} ){
	my ( $progs, $proj, $p_status, $proj_start, $proj_dir, $log, $config );
	my $num_cpu_used = 0;
	foreach my $i ( sort {$list->{$a}->{TIME} cmp $list->{$b}->{TIME}} keys %$list ) {
		$proj     = $list->{$i}->{NAME};
		$proj_dir = "$out_dir/$proj";
		my $run=0;
		$run = &availableToRun($list->{$i}->{CPU}, $num_cpu_used ) if $list->{$i}->{STATUS} eq "unstarted";
		if($run){
			my $json = `$RealBin/edge_action.cgi $proj rerun "" "" "" $sid $domain 2>> $proj_dir/error.log`;
			#print STDERR "$json";
			my $info = decode_json($json);
			$list->{$i}->{STATUS} = "running" if $info->{STATUS} == "SUCCESS";
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
		elsif( /^\[(.*)\]/ ){
			my $step = $1;
			next if $step eq "system";

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
			$prog->{$cnt}->{DO}=$do;
			$prog->{$cnt}->{STATUS}="unfinished";
			$prog->{$cnt}->{STATUS}="skip" if $do eq 0;
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

	return ($proj_status,$prog,$proj_start,$numcpu,$proj_desc,$proj_name);
}

sub scanNewProjToList {
	my $list_idx;
	my $cnt = (sort {$b<=>$a} keys %$list)[0];
	
	foreach my $i ( keys %$list ) {
		my $n = $list->{$i}->{NAME};
		$list_idx->{$n}=$i;
	}
	
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	while( defined (my $file = readdir BIN) ) {
		next if $file eq '.' or $file eq '..';
		if ( -d "$out_dir/$file" ) {
			next if defined $list_idx->{$file};
			$list->{++$cnt}->{NAME} = $file if -r "$out_dir/$file/process.log";
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
	my $mem = `free -m | awk 'NR==3{printf "%.1f", \$3*100/(\$4+\$3)}'`;
	my $cpu = `top -bn1 | grep load | awk '{printf "%.1f", \$(NF-2)}'`;
	my $disk = `df -h $out_dir | tail -1 | awk '{printf "%.1f", \$5}'`;
	$cpu = $cpu/$sys->{edgeui_tol_cpu}*100;
	$disk =~ s/\%//;
	if( $mem || $cpu || $disk ){
		$mem = sprintf "%.1f", $mem;
		$cpu = sprintf "%.1f", $cpu;
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
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		next if (! -r "$out_dir/$id/process.log");
		next if ( $status =~ /delete/i);
		$list->{$id}->{NAME} = $id;
		$list->{$id}->{PROJNAME} = $project_name;
		$list->{$id}->{DBSTATUS} = $status;
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{OWNER_FisrtN} = $hash_ref->{owner_firstname};
		$list->{$id}->{OWNER_LastN} = $hash_ref->{owner_lastname};
		$list->{$id}->{PROJTYPE} = $hash_ref->{type} if ($username && $password);
	}
}

sub getProjInfoFromDB{
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
	print $result_json if (@ARGV);
	my $hash_ref = from_json($result_json);

	my $id = $hash_ref->{id};
	my $project_name = $hash_ref->{name};
	my $status = $hash_ref->{status};
	my $projtype = ($hash_ref->{isPublished})?"publish":"false";
	#next if (! -r "$out_dir/$id/process.log");
	#next if ( $status =~ /delete/i);
	$list->{$id}->{NAME} = $id;
	$list->{$id}->{PROJNAME} = $project_name;
	$list->{$id}->{DBSTATUS} = $status;
	$list->{$id}->{PROJTYPE} = $projtype;
}

sub returnStatus {
	my $json;
	$json = to_json( $info ) if $info;
	$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
	print $cgi->header('application/json'), $json;
	exit;
}

