#!/usr/bin/env perl
#
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-14
#

use strict;
use JSON;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use Data::Dumper;

use LWP::UserAgent;
use HTTP::Request::Common;
use POSIX qw(strftime);
use Digest::MD5 qw(md5_hex);
require "edge_user_session.cgi";



my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
$pname ||= $ARGV[0];
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

my $username    = $opt{'username'}|| $ARGV[1];
my $password    = $opt{'password'}|| $ARGV[2];
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[3];
my $protocol = $opt{protocol} || 'http:';
my $sid         = $opt{'sid'}|| $ARGV[4];
my $viewType;
my $domain      = $ENV{'HTTP_HOST'}|| 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");

my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;
my $um_url      = $sys->{edge_user_management_url};
my $out_dir     = $sys->{edgeui_output};
my $projDir;
# Generates the project list (pname = encoded name) Scans output dir
if( $sys->{user_management} ){
	($username,$password,$viewType) = getCredentialsFromSession($sid);
	my $proj_code=&getProjCodeFromDB($pname, $username, $password);
	$projDir = $relpath ."/". $proj_code;
} else {
	if ( -e "$edgeui_output/$pname/config.txt"){
		$projDir = $relpath."/".$pname;
	}else{
		my $proj_list=&scanProjToList($edgeui_output);
		if( $proj_list->{$pname} ){
			$projDir = $relpath . "/". $proj_list->{$pname};
		}
	}
}


generateReport($projDir);


exit;

######################################################
sub generateReport {
	my $projDir = shift;
	chdir $edgeui_wwwroot;
	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/munger/outputMunger_w_temp.pl $projDir $projDir/HTML_Report/report_web.html";
	`$cmd`;
	#print STDERR "$cmd";

	open REP, "$projDir/HTML_Report/report_web.html" or die "Can't open report_web.html: $!";
	my $pr=0;
	my @htmls;
	while(<REP>){
		last if /<!-- \/content -->/;
		push @htmls, $_ if $pr;
		$pr=1 if /id='edge-content-report'/;
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
			open ( CONFIG, "$out_dir/$file/config.txt") or die "Cannot open $out_dir/$file/config.txt\n";
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
	
	return $projCode;
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

