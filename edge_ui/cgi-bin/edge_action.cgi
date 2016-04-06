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
require "edge_user_session.cgi";

my $cgi    = CGI->new;
my %opt    = $cgi->Vars();
my $pname  = $opt{proj};
my $username = $opt{username};
my $password = $opt{password};
my $action = lc($opt{action});
my $shareEmail = $opt{shareEmail};
my $userType = $opt{userType}||"user";
my $protocol = $opt{protocol}||"http:";
my $sid = $opt{sid};
my $domain	= $ENV{'HTTP_HOST'};

$pname      ||= $ARGV[0];
$action     ||= $ARGV[1];
$username   ||= $ARGV[2];
$password   ||= $ARGV[3];
$shareEmail ||= $ARGV[4];
$sid        ||= $ARGV[5];
$domain     ||= $ARGV[6];

# read system params from config template
my $config_tmpl = "$RealBin/edge_config.tmpl";
my $sys         = &getSysParamFromConfig($config_tmpl);
my $out_dir     = $sys->{edgeui_output};
my $input_dir   = $sys->{edgeui_input};
my $um_url      = $sys->{edge_user_management_url};
$domain       ||= "edgeset.lanl.gov";
$um_url	      ||= "$protocol//$domain/userManagement";
$out_dir      ||= "/tmp"; #for security
my $info;
my $proj_dir    = abs_path("$out_dir/$pname");
my $list;
my $permission;

#check projects vital
my ($vital, $name2pid) = &checkProjVital();
my $time = strftime "%F %X", localtime;

$info->{STATUS} = "FAILURE";
#$info->{INFO}   = "Project $pname not found.";


#session check
my $real_name = $pname;
my $user_proj_dir;
if ( $sys->{user_management} )
{
	my $valid = verifySession($sid);
	unless($valid){
		$info->{INFO} = "ERROR: Invalid session found.";
		&returnStatus();
	}
	else{
		($username,$password) = getCredentialsFromSession($sid);
	}
	
	$list = &getUserProjFromDB("owner");

	$real_name=getProjNameFromDB($pname);

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
	}
	#print STDERR "User: $username; Sid: $sid; Valid: $valid; Pname: $pname; Realname: $real_name; List:",Dumper($list),"\n";
}

if( $action eq 'empty' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	if( $name2pid->{$pname} ){
		$info->{INFO} = "ERROR: Project $real_name is running.";
		&returnStatus();
	}

	if( -d $proj_dir ){
		opendir(BIN, $proj_dir) or die "Can't open $proj_dir: $!";
		while( defined (my $file = readdir BIN) ) {
			next if $file eq '.' or $file eq '..';
			`rm -rf $proj_dir/$file` if -d "$proj_dir/$file";
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

		my $pid = $name2pid->{$pname};
		if( $pid ){
			my $invalid = &killProcess($pid);
			if( $invalid ){
				$info->{INFO} = "Failed to kill the running process. (PID: $pid)";
			}
		}
	
		`rm -rf $proj_dir`;
		`rm -rf $out_dir/$pname`;
		if( !-e $proj_dir && !-e "$out_dir/$pname" ){
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "Project $real_name has been deleted.";
		}
		if ($username && $password){
			&updateDBProjectStatus($pname,"delete");
			`rm -f $user_proj_dir`;
			`rm -f $input_dir/public/projects/${real_name}_$pname`;
			`rm -f $input_dir/*/SharedProjects/${real_name}_$pname`;
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
	
	my $pid = $name2pid->{$pname};

	if( $pid ){
		my $invalid = &killProcess($pid);
		if( !$invalid ){
			`echo "\n*** [$time] EDGE_UI: This project has been interrupted. ***" |tee -a $proj_dir/process.log >> $proj_dir/process_current.log`;
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "The process (PID: $pid) has been stopped.";
		}
	}
	else{
		$info->{INFO} = "Project $real_name is not running.";
	}
}
elsif( $action eq 'rerun' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}

	$info->{STATUS} = "FAILURE";
	$info->{INFO}   = "Failed to rerun project $real_name.";

	my $pid = $name2pid->{$pname};

	if( ! defined $pid ){
		my $cmd = "";
		open LOG, "$proj_dir/process.log" or die "Can't open process log:$!.";
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
			 &returnStatus();
		}
		if( $cmd ){
			chdir($proj_dir);
			my $newpid = open RUNPIPLINE, "-|", "$cmd > $proj_dir/process_current.log 2>&1 &" or die $!;
			close RUNPIPLINE;
			if( $newpid ){
				$newpid++;
				$info->{STATUS} = "SUCCESS";
				$info->{INFO}   = "Project $real_name has been restarted (PID: $newpid).";
				$info->{PID}    = $newpid;
			}
			else{
				$info->{INFO} = "Failed to restart this project.";
			}
		}
		else{
			$info->{INFO} = "Failed to restart this project. No runPipeline command found.";
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
			&updateDBProjectStatus($pname,"archive");
			`rm -f $user_proj_dir`;
			`rm -f $input_dir/public/projects/${real_name}_$pname`;
			`rm -f $input_dir/*/SharedProjects/${real_name}_$pname`;
		}
	}
}
elsif( $action eq 'share' || $action eq 'unshare' ){
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	&shareProject($pname,$shareEmail,$action);
}
elsif( $action eq 'publish' || $action eq 'unpublish'){
	print STDERR "USERMANAGMENT: $sys->{user_management}; $action: $permission->{$action}";
	if( $sys->{user_management} && !$permission->{$action} ){
		$info->{INFO} = "ERROR: Permission denied. Only project owner can perform this action.";
		&returnStatus();
	}
	&publishProject($pname,$action);
	my $public_proj_dir = "$input_dir/public/projects/${real_name}_$pname";
	`ln -sf $out_dir/$pname $public_proj_dir` if ($action eq 'publish' && ! -e "$public_proj_dir");
	`rm -f $public_proj_dir` if ($action eq 'unpublish');
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
	return 0 if ($num_cpu > $sys->{edgeui_tol_cpu});
	if( $sys->{edgeui_auto_queue} && $sys->{edgeui_tol_cpu} ){
		foreach my $pid ( keys %$vital ){
			$cpu_been_used += $vital->{$pid}->{CPU};
			return 0 if (($cpu_been_used + $num_cpu) > $sys->{edgeui_tol_cpu});
		}       
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
		return $result->{name};
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

sub shareProject{
	my $project=shift;
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
				`ln -sf $out_dir/$project $shared_proj_dir` if (!-e $shared_proj_dir);
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
		my $status = $hash_ref->{status};
		next if ($status =~ /delete/i);
		next if (! -r "$out_dir/$id/process.log");
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{REAL_PROJNAME} = $project_name;
		$list->{$id}->{OWNER} = $hash_ref->{owner_firstname};
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};
	}
	return $list;
}


