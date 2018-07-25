#!/usr/bin/env perl
use strict;
use LWP::Simple qw();
use FindBin qw($RealBin);
use lib "../../lib";
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use Digest::MD5 qw(md5_hex);
use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser);
use POSIX qw(strftime);
use Data::Dumper;
use Email::Valid;
#use CGI::Pretty;
require "./edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $username    = $opt{'username'}|| $ARGV[0];
my $password    = $opt{'password'}|| $ARGV[1];
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[2];
my $userType    = $ARGV[3]||"user";
my $viewType    = $opt{'view'}|| $ARGV[4];
my $protocol    = $opt{protocol}||'http:';
my $sid         = $opt{'sid'}|| $ARGV[5];
my $loadNum    = $opt{'loadNum'} || 9999999;
my $forceupdate = $opt{'forceupdate'} || $ARGV[6];
my $domain      = $ENV{'HTTP_HOST'} || 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
my $www_root    = $sys->{edgeui_wwwroot};
my $out_dir     = $sys->{edgeui_output};
my $um_config	= $sys->{user_management};
my $um_url      = $sys->{edge_user_management_url};
$um_url ||= "$protocol//$domain/userManagement";

my $cluster     = $sys->{cluster};
# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid || scalar(@ARGV) > 1 ) {
		($username,$password,$userType) = getCredentialsFromSession($sid) unless (scalar(@ARGV) > 1);
	}
}

#print Dumper ($list);
my $relt;
print $cgi->header('application/json');
my $html = "<h2>Project List</h2>";

if ( $username && $password || $um_config == 0){
	#Action buttons
	$html .= "<div id='edge-projpage-action' class='flex-container'>";
	if ($userType =~ /admin/i){
		$html .= "<a href='' title='See All Projects List (admin)' class='tooltip ui-btn ui-btn-d ui-icon-bars ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>show-all</a>";
	}
	$html .= "<a href='' title='Force Selected Projects to rerun' class='tooltip ui-btn ui-btn-d ui-shadow-icon ui-icon-refresh ui-btn-icon-notext ui-corner-all' data-role='button' >rerun</a>";
	$html .= "<a href='' title='Interrupt running Projects' class='tooltip ui-btn ui-btn-d ui-icon-forbidden ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>interrupt</a>";
	$html .= "<a href='' title='Delete Selected Projects' class='tooltip ui-btn ui-btn-d ui-icon-delete ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>delete</a>";
	$html .= "<a href='' title='Empty Selected Projects Output' class='tooltip ui-btn ui-btn-d ui-icon-recycle ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>empty</a>";
	if ($sys->{edgeui_archive}){
		$html .= "<a href='' title='Archive Selected Projects' class='tooltip ui-btn ui-btn-d ui-icon-arrow-u-r ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>archive</a>";
 	}
	 if ($um_config != 0){
		$html .= "<a href='' title='Share Selected Projects' class='tooltip ui-btn ui-btn-d ui-icon-forward ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>share</a>";
		$html .= '<a href="" title="Make Selected Projects Public" class="tooltip ui-btn ui-btn-d ui-icon-eye ui-btn-icon-notext ui-corner-all" data-role="button" role="button">publish</a>';
 	}
	$html .= "<a href='' title='Compare Selected Projects Taxonomy Classification (HeatMap)' class='tooltip ui-btn ui-btn-d ui-icon-bullets ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>compare</a>";
 	if($sys->{edge_sample_metadata}) {
 		$html .= "<a href='' title='Export Selected Projects Metadata' class='tooltip ui-btn ui-btn-d ui-icon-arrow-d ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>metadata-export</a>";
 	}
 	if($sys->{edge_sample_metadata} && $sys->{edge_sample_metadata_share2bsve}) {
 		$html .= "<a href='' title='Share Selected Projects Metadata/Pathogens with BSVE' class='tooltip ui-btn ui-btn-d ui-icon-arrow-u ui-btn-icon-notext ui-corner-all' data-role='button' role='button'>metadata-bsveadd</a>";
 	}
	
 	$html .= "</div>";
}

$html .= "<table id='edge-project-page-table' class='output-table ui-responsive ui-table ui-table-reflow'>";
my $head_checkbox="<input type='checkbox' id='edge-projpage-ckall'>";

if ($umSystemStatus=~ /true/i && $username && $password ){
	if($viewType =~ /user/i ) {
		&getUserProjList();
	} else {
		&getAdminProjList();
	}

}elsif ($umSystemStatus=~ /true/i) {
	&getPublicProjList();
}elsif ($um_config == 0) {
	&getProjList();
}

$html .= "</table>";

$relt -> {html} = $html;
my $relt_json = encode_json($relt);
#print STDERR $relt_json;
print $relt_json;


## END MAIN## 
sub getUserProjList {
	$html .= "<thead><tr>";
	$html .= "<th>$head_checkbox</th><th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th><th>Type</th><th>Owner</th>";
	$html .= " </tr></thead>";
	&getProjFromUM("user");
}

sub getAdminProjList {
	$html .= "<thead><tr>";
	$html .= "<th>$head_checkbox</th><th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th><th>Owner</th>";
	$html .= " </tr></thead>";
	&getProjFromUM("admin");
}

sub getPublicProjList {
	$html .= "<thead><tr>";
	$html .= "<th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th><th>Owner</th>";
	$html .= " </tr></thead>";
	&getProjFromUM("public");
}

sub getProjList {	
	# all projects in the EDGE_output
	$html .= "<thead><tr>";
	$html .= "<th>$head_checkbox</th><th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th>";
	$html .= " </tr></thead>";
	
	my $list= &scanProjToList();
	my ($idxs)= &sortList($list);
	&tdList($idxs,$list);
}

sub getProjFromUM {
	my $request_type = shift;
	my @tds_for_list;
	
        my %data = (
                email => $username,
                password => $password
        );
	my $api_path;
	if($request_type eq "user") {
		$api_path = "user/getRuns";
	} elsif($request_type eq "admin") {
		$api_path = "user/admin/getRuns";
	} else {
		$api_path = "user/publicRuns";
	}       

	my $service="WS/$api_path";
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
	#print STDERR "$result_json";

	if ($result_json =~ /\"error_msg\":"(.*)"/)
        {                
                return "ERROR";
        }
        my $array_ref =  decode_json($result_json);

	my @projectlist=@$array_ref;
	my ($id, $code, $name, $status, $submitted, $running_time,$type, $owner);
	foreach my $hash_ref (@projectlist)
	{
		$id = $hash_ref->{id};
		$code = $hash_ref->{code};
		$name = $hash_ref->{name};
		$status = $hash_ref->{status};
		$submitted = $hash_ref->{submitted};
		$running_time = $hash_ref->{running_time};
		$type = $hash_ref->{type};
		$owner = $hash_ref->{owner};
		
		next if ( ! -e "$out_dir/$code" and ! -e "$out_dir/$id");

		my $projname = "<a href='#' class='edge-project-page-link ' title='alt-click to open in a new tab' data-pid='$id'>$name</a>";
		my $checkbox = "<input type='checkbox' class='edge-projpage-ckb' name='edge-projpage-ckb' value='$code'>";

		if($request_type eq "user") {
			@tds_for_list = ( $checkbox,$projname,$status,$submitted,$running_time,$type,$owner );
		} elsif($request_type eq "admin") {
			@tds_for_list = ( $checkbox,$projname,$status,$submitted,$running_time,$owner );
		} else {
			@tds_for_list = ( $projname,$status,$submitted,$running_time,$owner );
		}
		push @{$relt->{data}}, [@tds_for_list];
	
	}
}


sub scanProjToList {
	my $cnt = 0;
	my $list;
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	my @dirs= readdir BIN;
	foreach my $file (@dirs) {
		next if $file eq '.' or $file eq '..' or $file eq 'sample_metadata_export';
		if ( -d "$out_dir/$file"  && -r "$out_dir/$file/config.txt") {
			$cnt ++;
			if ( -e "$out_dir/$file/.AllDone" && -e "$out_dir/$file/HTML_Report/writeHTMLReport.finished"){
				$list=&get_start_run_time("$out_dir/$file/.AllDone",$cnt,$list);
				$list->{$cnt}->{PROJSTATUS} = "Complete";
			}else{
				if (-r "$out_dir/$file/process.log" ){
                                	$list=&pull_summary("$out_dir/$file/process.log",$cnt,$list,"$out_dir/$file")
                               	}else{
                                  	$list->{$cnt}->{PROJSTATUS} = "Unstarted";
                             	}
				$list=&pull_summary("$out_dir/$file/config.txt",$cnt,$list,"$out_dir/$file") if ($list->{$cnt}->{PROJSTATUS} =~ /unstart/i);
				
			}
			$list->{$cnt}->{REAL_PROJNAME} = $list->{$cnt}->{PROJNAME} || $file;
			$list->{$cnt}->{PROJNAME} = $file;
		}
	}
	closedir(BIN);
	return $list;
}

sub sortList {
	my $list = shift;
	
	my @idxs1 = grep { $list->{$_}->{PROJSTATUS} =~ /running/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs2 = grep { $list->{$_}->{PROJSTATUS} =~ /unstarted/i } sort {$list->{$a}->{REAL_PROJNAME} cmp $list->{$b}->{REAL_PROJNAME}} keys %$list;
	my @idxs3 = grep { $list->{$_}->{PROJSTATUS} !~ /running|unstarted/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs = (@idxs1,@idxs2,@idxs3);
	
	return (\@idxs);
}

sub tdList {
	my $idx_ref = shift;
	my $list = shift;
	my @idxs = ($idx_ref)? @{$idx_ref}:"";
	my @tds_for_list;
	#return if (@ARGV);
	if ($list->{INFO}->{ERROR})
	{
		print "<p class='error'>$list->{INFO}->{ERROR}</p>\n";
	}
	foreach (@idxs)
	{
		my $projStatus = $list->{$_}->{PROJSTATUS};
		my $projID = $list->{$_}->{PROJNAME};
		my $projname = "<a href='#' class='edge-project-page-link ' title='$list->{$_}->{PROJDESC} (alt-click to open in a new tab)' data-pid='$projID'>$list->{$_}->{REAL_PROJNAME}</a>";
		my $projSubTime = $list->{$_}->{PROJSUBTIME};
		my $projRunTime = $list->{$_}->{RUNTIME};
		my $projCode = $list->{$_}->{PROJCODE} || $list->{$_}->{REAL_PROJNAME};
		my $checkbox = "<input type='checkbox' class='edge-projpage-ckb' name='edge-projpage-ckb' value='$projCode'>";
		@tds_for_list = ( $checkbox,$projname,$projStatus,$projSubTime,$projRunTime ); 
		push @{$relt->{data}}, [@tds_for_list];
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
			$list->{$cnt}->{PROJSTATUS} = "<span class='edge-fg-orange'>Running</span>";
		}
		elsif( /failed/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "<span class='edge-fg-red'>Failed</span>";
			$list->{$cnt}->{PROJSTATUS} = "<span class='edge-fg-red'>Failure</span>";
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
