#!/usr/bin/env perl 
#
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-14
#

use strict;
use JSON;
use FindBin qw($RealBin);
use lib "../../lib";
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use Data::Dumper;
use Email::Valid;
use LWP::UserAgent;
use HTTP::Request::Common;
use POSIX qw(strftime);
use Digest::MD5 qw(md5_hex);
require "./edge_user_session.cgi";

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

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
my $domain      = $ENV{'HTTP_HOST'}|| 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;
&stringSanitization(\%opt);

exit if (!$pname);
# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};

my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;
my $um_url      = $sys->{edge_user_management_url};
my $out_dir     = $sys->{edgeui_output};
my $projDir;
my $proj_code;
my $proj_status;
# Generates the project list (pname = encoded name) Scans output dir
if( $sys->{user_management} && $pname !~ /\D/ && $sid ){
	($username,$password,$viewType) = getCredentialsFromSession($sid);
	($proj_code,$proj_status)=&getProjCodeFromDB($pname, $username, $password);
	if (!$proj_code){
		my $html = "<p>The project does not exist</p>";
		print "Content-Type: text/html\n\n",
		$html;
		exit 0;
	}
	$projDir = (-d "$edgeui_output/$pname")? "$relpath/$pname": $relpath ."/". $proj_code; 
}
if( !$sys->{user_management} || !$username ){
	if ( -e "$edgeui_output/$pname/config.txt"){
		$projDir = $relpath."/".$pname;
	}else{
		my $proj_list=&scanProjToList($edgeui_output);
		if( $proj_list->{$pname} ){
			$projDir = $relpath . "/". $proj_list->{$pname};
		}else{
			my $html = "<p>The project does not exist</p>";
			print "Content-Type: text/html\n\n",
			$html;
			exit 0;
		}
	}
}
generateReport($projDir);


exit;

######################################################
sub generateReport {
	my $projDir = shift;
	my $complete_report_exist = 0;
	chdir $edgeui_wwwroot;

	if( ! -e "$projDir/HTML_Report/.complete_report_web" ){
		system("$EDGE_HOME/scripts/munger/outputMunger_w_temp.pl","$projDir","$projDir/HTML_Report/report_web.html");
		#print STDERR "outputMunger_w_temp.pl ran!\n";
	}
	else{
		$complete_report_exist = 1;
		#print STDERR "outputMunger_w_temp.pl NOT ran!\n";
	}

	open REP, "<", "$projDir/HTML_Report/report_web.html" or die "Can't open report_web.html: $!";
	my $pr=0;
	my @htmls;
	foreach(<REP>){
		last if /<!-- \/content -->/;
		push @htmls, $_ if $pr;
		$pr=1 if /id='edge-content-report'/;
		
		if( $_ =~ /edge-output-projstatus.*(complete|archived)/i && !$complete_report_exist ){
			open (my $fh,">","$projDir/HTML_Report/.complete_report_web") or die "$!"; close $fh;
		}
	}

	close REP;

	my $html = join "", @htmls;

	print "Content-Type: text/html\n\n",
		  $html;
}

sub scanProjToList {
	my $out_dir = shift;
	my $list;
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	while( defined (my $file = readdir BIN) ) 
	{
		next if $file eq '.' or $file eq '..';
		my $projid;
		my $projCode;
		if ( -d "$out_dir/$file" && -r "$out_dir/$file/config.txt"  ) 
		{
			open ( CONFIG, "<","$out_dir/$file/config.txt") or die "Cannot open $out_dir/$file/config.txt\n";
			while(<CONFIG>)
			{
				last if (/^\[Down/);
				if (/^projid=(\S+)/)
				{
					$projid=$1;
				}
				if (/^projcode=(\S+)/)
				{
					$projCode=$1;
				}
			}
			close CONFIG;
			$projid ||= $file;
			$list->{$projid} = $file;
			$list->{$file} = $file;
			$list->{$projCode} = $file;
		}
	}
	closedir(BIN);
	return $list;
}

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open CONF, $config or die "Can't open $config: $!";
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


sub getProjCodeFromDB{
	my $projectID=shift;
	$projectID = &getProjID($projectID);
	my $username = shift;
	my $password = shift;

	my %data = (
		email => $username,
		password => $password,
		project_id => $projectID 
	);
    
	$um_url ||= "$protocol//$domain/userManagement";
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
	# print $result_json if (@ARGV);
	my $hash_ref = from_json($result_json);

	my $id = $hash_ref->{id};
	my $project_name = $hash_ref->{name};
	my $projCode = $hash_ref->{code};
	my $projStatus = $hash_ref->{status};
	
	return ($projCode,$projStatus);
}
sub getProjID {
  my $project=shift;
  my $projID = $project;
  if ( -d "$out_dir/$project"){ # use ProjCode as dir
    open (my $fh, "<","$out_dir/$project/config.txt") or die "Cannot open $out_dir/$project/config.txt\n";
    while(<$fh>){
      if (/^projid=(\S+)/){
        $projID = $1;
        last;
      }
    }
  }
  return $projID;
}

sub stringSanitization{
	my $opt=shift;
	my $dirtybit=0;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		next if $key eq "password";
		next if $key eq "keywords";

		if ($key eq "username" || $key eq "shareEmail"){
			my @email = split(',',$str);
			map { $dirtybit =1  if ! &emailValidate($_); } @email;
		}else{
			$opt->{$key} =~ s/[`;'"]/ /g;
			$dirtybit=1 if ($opt->{$key} =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/);
		}
		if ($dirtybit){
			print "Content-Type: text/html\n\n", "Invalid characters detected \'$str\'.\n\n";
			exit;
		}
	}
}

sub emailValidate{
	my $email=shift;
	$email = lc($email);
	return Email::Valid->address($email);
}

