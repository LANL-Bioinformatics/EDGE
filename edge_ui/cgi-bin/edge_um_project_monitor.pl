#!/usr/bin/env perl
use strict;
use JSON;
use POSIX qw(strftime);
use FindBin qw($RealBin);
use LWP::UserAgent;
use HTTP::Request::Common;
use DBI;

my $in_status= $ARGV[0];
# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $username = $sys->{edgeui_admin};
my $password = $sys->{edgeui_admin_password};
my $out_dir = $sys->{edgeui_output};
my $umURL = $sys->{edge_user_management_url};
exit if ( $ENV{"REQUEST_METHOD"} );
my $numProj=0;
my $list={};
my %data = (
   email => $username,
   password => $password,
);

if($in_status) {
	$data{project_status} = $in_status;
}
# Encode the data structure to JSON
my $json = JSON->new;
my $data = $json->encode(\%data);
#print "$data\n";


# Set the request parameters
my $url = $umURL.'WS/user/admin/getProjects';

my $browser = LWP::UserAgent->new;
my $req = PUT $url;
$req->header('Content-Type' => 'application/json');
$req->header('Accept' => 'application/json');
#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
$req->header( "Content-Length" => length($data) );
$req->content($data);

my $response = $browser->request($req);
my $result_json = $response->decoded_content;
	
if ($result_json =~ /\"error_msg\":"(.*)"/){
      	$list->{INFO}->{ERROR}=$1;
	print "Failed to get projects from UM\n";
     	return;
}

my $array_ref =  decode_json($result_json);
my @projectlist=@$array_ref;
foreach my $hash_ref ( @projectlist) {
	my $id = $hash_ref->{id};
	my $projCode = $hash_ref->{code};
	my $project_name = $hash_ref->{name};
	my $status = $hash_ref->{status};
	my $created_time = $hash_ref->{created};

	next if (! -r "$out_dir/$id/process.log" && ! -r "$out_dir/$projCode/process.log" );
	my $proj_dir=(-d "$out_dir/$projCode")?"$out_dir/$projCode":"$out_dir/$id";
	my $processlog = (-r "$proj_dir/process.log")? "$proj_dir/process.log":"$proj_dir/config.txt";
print "$processlog\n";
		
	if ( $status =~ /finished|archived/i && -e "$proj_dir/.AllDone"){
		$list=&get_start_run_time("$proj_dir/.AllDone",$id,$list);
		$status=($status =~ /finished/i)?"Complete":"Archived";
	}elsif( $status =~ /unstarted/i){
		$list->{$id}->{PROJSUBTIME}=$created_time;
		$list->{$id}->{RUNTIME}="";
		$status="Unstarted";
	}elsif( $status =~ /in process/i){
       		$list->{$id}->{PROJSUBTIME}=$created_time;
      		$list=&pull_summary($processlog,$id,$list,$proj_dir) if (  -e $processlog);
                        #$list->{$id}->{RUNTIME}="";
             	$status="In process";
	}else{
		$list=&pull_summary($processlog,$id,$list,$proj_dir) if (  -e $processlog);
	}

	$list->{$id}->{PROJNAME} = $id;
	if (!$list->{$id}->{PROJSTATUS}){
		$list->{$id}->{PROJSTATUS} = $status;
		$list->{$id}->{PROJSTATUS} = "Running" if ($status =~ /running/);
		$list->{$id}->{PROJSTATUS} = "Failure" if ($status =~ /failure/);
	}
	$list->{$id}->{PROJSUBTIME}=$created_time if (!$list->{$id}->{PROJSUBTIME});
	$list->{$id}->{REAL_PROJNAME} = $project_name if (!$list->{$id}->{REAL_PROJNAME});
	$list->{$id}->{PROJCODE} = $projCode;
	$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}";
	$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
	$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};

	$numProj++;
}

##update project in UM db
foreach (keys %$list) {
	my $projStatus = $list->{$_}->{PROJSTATUS};
	my $projID = $list->{$_}->{PROJNAME};
	my $projTime = $list->{$_}->{TIME};
	my $projRunTime = $list->{$_}->{RUNTIME};

	if($projStatus eq "Complete") {
		&updateProject($projID, $projStatus, $projTime, $projRunTime);
	} else {
		&updateProject($projID, $projStatus, $projTime);
	}
}

print "updated projects: $numProj\n";
###################
sub updateProject{
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
	my $data = to_json(\%data);

        my $url = $umURL ."WS/project/update";
        my $browser = LWP::UserAgent->new;
        my $req = PUT $url;
        $req->header('Content-Type' => 'application/json');
        $req->header('Accept' => 'application/json');
       
        $req->header( "Content-Length" => length($data) );
        $req->content($data);

        my $response = $browser->request($req);
        my $result_json = $response->decoded_content;
        my $result =  from_json($result_json);
        if (! $result->{status})
        {
                print " Update Project status in database failed.";
        }
}   
     
sub get_start_run_time{
	my $log = shift;
	my $id = shift;
	my $list = shift;
	my $tol_running_sec=0;
	my ($start_time, $run_time);
	open (my $log_fh,"<", $log) or die "Cannot read $log\n";
	while(<$log_fh>){
		chomp;
		($start_time, $run_time) = split /\t/,$_;
		$list->{$id}->{PROJSUBTIME} = $start_time if ($start_time);
		if ($run_time =~ /(\d+):(\d+):(\d+)/ ){
			$list->{$id}->{LASTRUNTIME} = "$1h $2m $3s";
                        $tol_running_sec += $1*3600+$2*60+$3;
		}
	}
	my ($yyyy,$mm,$dd,$hms) = $list->{$id}->{PROJSUBTIME} =~ /(\d{4}) (\w{3})\s+(\d+)\s+(.*)/;
	my %mon2num = qw(jan 1  feb 2  mar 3  apr 4  may 5  jun 6  jul 7  aug 8  sep 9  oct 10 nov 11 dec 12);
	$mm = $mon2num{ lc substr($mm, 0, 3) };
	$mm = sprintf "%02d", $mm;
	$dd = sprintf "%02d", $dd;
	$list->{$id}->{TIME} = "$yyyy-$mm-$dd $hms";
	$list->{$id}->{RUNTIME} = sprintf("%02d:%02d:%02d", int($tol_running_sec / 3600), int(($tol_running_sec % 3600) / 60), int($tol_running_sec % 60));
	close $log_fh;
	return $list;
}
sub pull_summary {
	my $log = shift;
	my $cnt= shift;
	my $list = shift;
	my $proj_dir =shift;
	my @INFILES;
	
	my ($step,$lastline);
	my $tol_running_sec=0;
	my $run_time;
	open(my $sumfh, "<", "$log") or die "unable to open $log: $!";
	while(<$sumfh>) {
		chomp;

		if( /Total Running time: (\d+):(\d+):(\d+)/){
			$run_time = $_;
			$list->{$cnt}->{LASTRUNTIME} = "$1h $2m $3s";
			next;
		}
		if( /^\[(.*)\]/ ){
			$step = $1;
			if( $step eq "project" or $step eq "system"){
				while(<$sumfh>){
					chomp;
					if ( /^([^=]+)=([^=]+)/ ){
						$list->{$cnt}->{uc($1)}=$2;
						$list->{$cnt}->{REAL_PROJNAME}=$2 if ($1 eq "projname");
					}
					elsif ( /^\[(.*)\]/ ){
						$step = $1;
						last;
					}
				}
			}
			
		}
		elsif( /Project Start: (.*)/ ){
			$list->{$cnt}->{PROJSUBTIME} = $1;
			my ($yyyy,$mm,$dd,$hms) = $list->{$cnt}->{PROJSUBTIME} =~ /(\d{4}) (\w{3})\s+(\d+)\s+(.*)/;
			my %mon2num = qw(jan 1  feb 2  mar 3  apr 4  may 5  jun 6  jul 7  aug 8  sep 9  oct 10 nov 11 dec 12);
			$mm = $mon2num{ lc substr($mm, 0, 3) };
			$mm = sprintf "%02d", $mm;
			$dd = sprintf "%02d", $dd;
			my $proj_start  = "$yyyy-$mm-$dd $hms";
			$list->{$cnt}->{TIME} = $proj_start;
			$list->{$cnt}->{PROJSTATUS} = "Unfinished";
		}
		elsif( /^Do.*=(.*)$/ ){
			my $do = $1;
			$list->{$cnt}->{$step}->{GNLRUN}= "Auto";
			$list->{$cnt}->{$step}->{GNLRUN}= "On" if $do eq 1;
			$list->{$cnt}->{$step}->{GNLRUN}= "Off" if $do eq 0;
			
			$list->{$cnt}->{$step}->{GNLSTATUS}="Skipped";
			$list->{$cnt}->{$step}->{GNLSTATUS}="Incomplete" if $do eq 1;
		}
		elsif( /Finished/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "Skipped (result exists)";
		}
		elsif( /Running time: (\d+:\d+:\d+)/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "Complete";
			$list->{$cnt}->{$step}->{GNLTIME} = $1;
        	        my ($h,$m,$s) = $1 =~ /(\d+):(\d+):(\d+)/;
                	$tol_running_sec += $h*3600+$m*60+$s;
			$list->{$cnt}->{PROJSTATUS} = "Unfinished";
		}
		elsif( / Running/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "<span class='edge-fg-orange'>Running</span>";
			$list->{$cnt}->{PROJSTATUS} = "Running";
		}
		elsif( /failed/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "<span class='edge-fg-red'>Failed</span>";
			$list->{$cnt}->{PROJSTATUS} = "Failure";
		}
		elsif( /All Done/){
			$list->{$cnt}->{PROJSTATUS} = "Complete";
		}
		$lastline = $_;
	}

        #$list->{$cnt}->{RUNTIME} = strftime("\%H:\%M:\%S", gmtime($tol_running_sec));
	$list->{$cnt}->{RUNTIME} = sprintf("%02d:%02d:%02d", int($tol_running_sec / 3600), int(($tol_running_sec % 3600) / 60), int($tol_running_sec % 60));

	$list->{$cnt}->{PROJSTATUS}        = "Unstarted"   if $lastline =~ /EDGE_UI.*unstarted/;
	$list->{$cnt}->{PROJSTATUS}        = "Interrupted" if $lastline =~ /EDGE_UI.*interrupted/;
	$list->{$cnt}->{PROJSTATUS}        = "Archived"    if $lastline =~ /EDGE_UI.*archived/;
	$list->{$cnt}->{TIME}              = $1            if $lastline =~ /\[(\S+ \S+)\] EDGE_UI/;
	$list->{$cnt}->{$step}->{GNLSTATUS} = "Interrupted" if $list->{$cnt}->{$step}->{PROJSTATUS} eq "Interrupted"; #turn last step to unfinished
	
	$list->{$cnt}->{INFILES} = join ", ", @INFILES;
	$list->{$cnt}->{TIME} ||= strftime "%F %X", localtime;
	
	close ($sumfh);
	if ($list->{$cnt}->{PROJSTATUS} eq "Complete" && ! -e "$proj_dir/.AllDone"){
		if  ($list->{$cnt}->{RUNTIME}){
			open (my $fh, ">","$proj_dir/.AllDone") or die "$!"; 
			print $fh "$list->{$cnt}->{PROJSUBTIME}\t$run_time";
			close $fh;
		}
	}
	return $list;
}

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open (CONF, "<", $config )or die "Can't open $config: $!";
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




