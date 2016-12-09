#! /usr/bin/env perl

use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use JSON;
use CGI qw(:standard);
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
require "edge_user_session.cgi";

my $cgi    = CGI->new;
my %opt    = $cgi->Vars();
my $username = $opt{username};
my $password = $opt{password};
my $action = lc($opt{action});
my @pname = split /,/,$opt{proj}; 
my $protocol = $opt{protocol};
my $sid    = $opt{sid};
$action   ||= $ARGV[0];
$username ||= $ARGV[1];
$password ||= $ARGV[2];
$pname[0] ||= $ARGV[3];
$sid      ||= $ARGV[4];

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $edgeui_admin = $sys->{edgeui_admin};
my $edgeui_adminpw = $sys->{edgeui_admin_password}; 
my $um_switch	 = $sys->{user_management};
my $um_url       = $sys->{edge_user_management_url};
my $edge_input   = $sys->{edgeui_input};
my $domain       = $ENV{'HTTP_HOST'};
my $sessionValid = 0;
$domain ||= "edgeset.lanl.gov";
$protocol ||= "http:";
$um_url ||= "$protocol//$domain/userManagement";
#print Dumper (\%ENV);
my $info;
my %data;

#check session
if( $sys->{user_management} ){
	$sessionValid = verifySession($sid);
	($username,$password) = getCredentialsFromSession($sid) if $sessionValid;
	#print STDERR "$0\nSID: $sid\nVALID: $valid\nUSER: $username\nPASS: $password\n";
}

if ($sys->{user_upload_maxFileSize}){
	$info->{maxFileSize} = $sys->{user_upload_maxFileSize};
}
if ($sys->{user_social_login}){
	$info->{socialLogin} = $sys->{user_social_login};
}

if ($action eq "check"){
	if ($um_switch){
		if (&check_um_service("$um_url")){
			$info->{SUCCESS} = "User management system ($um_url) is live ";
			$info->{url} = "$um_url";
		}else{
			$info->{error}  = "User management system ($um_url) is down. Please contact your system administrator.";
		}
	}else{
			$info->{error}  = "User management system ($um_url) is off ";
	}
	$info->{sid} = $sessionValid ? $sid : "";
}
elsif ($action eq "login"){
	%data = (
		email => $username, 
		password => $password,
	);
	my $data = to_json(\%data);
	&um_service($um_url,$data,"WS/user/login");
	&um_service($um_url,$data,"WS/user/getInfo");

	unless( $info->{error} ){
		my $sid = createSession($username,$password,$info->{type});
		$info->{SESSION} = $sid;
		my $user_dir=md5_hex(lc($username));
		$info->{UserDir} = $user_dir; 
		my $cronjobs = `crontab -l 2>/dev/null`;
		$info->{CleanData} = $sys->{edgeui_proj_store_days} if ($sys->{edgeui_proj_store_days} > 0 && $cronjobs =~ /edge_data_cleanup/);
		`mkdir -p $edge_input/$user_dir/MyProjects/`;
		`ln -sf $edge_input/public/data $edge_input/$user_dir/PublicData` if (! -e "$edge_input/$user_dir/PublicData");a
		`ln -sf $edge_input/public/projects $edge_input/$user_dir/PublicProjects` if (! -e "$edge_input/$user_dir/PublicProjects");
	}
}
elsif ($action eq "logout"){
	my $valid = closeSession($sid);

	$info->{error} = "failed to close the session:$sid\n" unless $valid;
}
elsif ($action eq "update"){
	%data = (
		email => $username, 
		password => $password
	);
	$data{new_firstname}=$opt{newFn};
	$data{new_lastname}=$opt{newLn};
	$data{new_password}=$opt{newPass} if ($opt{newPass});
	my $data = to_json(\%data);
	&um_service($um_url,$data,"WS/user/update");
}
elsif ($action eq "sociallogin"){
	
	my %admin_data = (
                admin_email => "$edgeui_admin",
                admin_password => "$edgeui_adminpw"
        );
	my $admin_data = to_json(\%admin_data);
	&um_service($um_url,$admin_data,"WS/user/admin/getUsers");
	my $social_email = $opt{email};
	my $social_id = $opt{id};
	my $social_fn = $opt{first_name};
	my $social_ln = $opt{last_name};
	my $auto_pass = md5_hex("$social_id$social_email");
	%data = (
                email => $social_email,
		password => $auto_pass
	);
	if ($info->{$social_email}){
		$info={};
		#try to login 
		my $data = to_json(\%data);
		&um_service($um_url,$data,"WS/user/getInfo");
		
		if ( $info->{error}){ # login fail, update the spassword and login
			$admin_data{email} = $social_email;
			$admin_data{new_spassword} = $auto_pass;
			$admin_data{new_active} = 'yes';
			$admin_data = to_json(\%admin_data);
 			&um_service($um_url,$admin_data,"WS/user/admin/update");
			&um_service($um_url,$data,"WS/user/getInfo");
		}
		$info->{password} = $auto_pass;
	}else{	# no existing account, auto register user by user social login info
		$data{firstname} = $social_fn;
		$data{lastname} = $social_ln;
		my $data = to_json(\%data);
		&um_service($um_url,$data,"WS/user/register");
		&um_service($um_url,$data,"WS/user/getInfo");
		$info->{password} = $auto_pass;
	}
}
elsif ($action eq "share" || $action eq "unshare"){
	my $valid = verifySession($sid);
	
	if( $valid ){
		my $service = ($action eq "share")? "WS/project/getNonguests":"WS/project/getGuests"; 
        %data = (
			email => $username,
            password => $password,
			project_id => $pname[0]
			);
		if (scalar(@pname)>1 || $pname[0] !~ /\d/ ){
			$service = "WS/user/admin/getUsers"; 
			%data = (
            	admin_email => $edgeui_admin,
            	admin_password => $edgeui_adminpw,
    	    	);
		}
    	my $data = to_json(\%data);
    	&um_service($um_url,$data,$service);
	}
	else{
		$info->{error} = "Invalid session.\n";
	}
}

#$info->{SUCCESS} = "login is successful";
#$info->{error}  = "username or password is wrong";

&returnStatus();

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

sub check_um_service {
	my $url=shift;
	$url .= "WS/user/hello";
#       $url = 'http://bioedge.lanl.gov/userManagementWS/user/hello';
        #print $url,"\n";
	my $browser = LWP::UserAgent->new;
 	my $req = GET $url;
	my $response = $browser->request($req);
	if ($response -> is_success) {
            return 1;
	}else{
        #       warn "The User managment Service is DOWN!!!! Will pull all projects from EDGE output direcotry";
            return 0;
    	}
}

sub um_service {
	my $url = shift;
	my $data = shift;
	my $service = shift;
	$url .= $service;
	my $browser = LWP::UserAgent->new;
	my $req = ($service =~ /register/)? POST $url : PUT $url;
	my $response = $browser->request($req);
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	print $result_json,"\n" if (@ARGV);
	if ($result_json =~ /\"error_msg\":"(.*)"/)
	{
		$info->{error}=$1;
		return;
    }
	if ($service =~ /login|getInfo/)
	{
		my $tmp_r = from_json($result_json);
		$info->{$_} = $tmp_r->{$_} foreach (keys %$tmp_r);
		
	}else{
		my $array_ref =  from_json($result_json);
		if ($action =~ /sociallogin/){
			foreach (@$array_ref){
				my $email=$_->{email};
				$info->{"$email"} = 1;
			}
		}else{
			$info =  $array_ref;
		}
	}
}

sub returnStatus {
    my $json;
    $json = to_json( $info ) if $info;
    $json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
    print $cgi->header('application/json'), $json;
    exit;
}
