#!/usr/bin/env perl
use strict;
use LWP::Simple qw();
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser);
use POSIX qw(strftime);
use Data::Dumper;
use CGI::Pretty;
require "edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $username    = $opt{'username'}|| $ARGV[0];
my $password    = $opt{'password'}|| $ARGV[1];
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[2];
my $userType    = $opt{'userType'}|| $ARGV[3];
my $viewType    = $opt{'view'}|| $ARGV[4];
my $protocol    = $opt{protocol}||'http:';
my $sid         = $opt{'sid'}|| $ARGV[5];

# read system params from config template
my $config_tmpl = "$RealBin/edge_config.tmpl";
my $sys         = &getSysParamFromConfig($config_tmpl);
my $out_dir     = $sys->{edgeui_output};
my $um_config	= $sys->{user_management};
my $um_url      = $sys->{edge_user_management_url};
my $domain      = $ENV{'HTTP_HOST'};
$domain ||= "edgeset.lanl.gov";
$um_url ||= "$protocol//$domain/userManagement";

# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid){
		($username,$password) = getCredentialsFromSession($sid);
	}
}

#print Dumper ($list);
print  $cgi->header( "text/html" );
if ($userType =~ /admin/i && $viewType !~ /admin/i){
	print "<h2>Project List	 <a id=\"edge-project-page-admin\" href=\"#\" class='ui-btn ui-mini ui-btn-inline ui-icon-bars ui-btn-icon-left'>See Project List as Admin</a></h2>\n";
}else{
	print  $cgi->h2("Project List");
}

print "<div data-filter='true' id='edge-project-list-filter' data-filter-placeholder='Search projects ...'> ";

if ($umSystemStatus=~ /true/i && $username && $password && $viewType =~ /user/i ){
	# My Table
	my $list = &getUserProjFromDB("owner");
	#<div data-role='collapsible-set' id='edge-project-list-collapsibleset'> 
	print <<"HTML";  
		<div data-role='collapsible' data-collapsed-icon="carat-d" data-expanded-icon="carat-u" data-collapsed="false" >
		<h4>My Projects</h4>
HTML
	;

	my @theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Submission Time"),th("Total Running Time"),th("Type"));
	#my @theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Total Running Time"),th("Type"),th("Action"));
	my $idxs = &sortList($list);
	my $table_id = "edge-project-page-Mytable";
	&printTable($table_id,$idxs,$list,\@theads);
	print "\t\t</div>\n";
	
	# Projects shared to me Table
	undef $list;
	undef $idxs;
	$list =  &getUserProjFromDB("guest");
	print <<"HTML";  
		<div data-role='collapsible' data-collapsed-icon="carat-d" data-expanded-icon="carat-u" data-collapsed="true" >
		<h4>Projects shared with me</h4>
HTML
;
	@theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Submission Time"),th("Total Running Time"),th("Owner"));
	$idxs = &sortList($list);
	$table_id = "edge-project-paga-Guesttable";
	&printTable($table_id,$idxs,$list,\@theads);
	print "\t\t</div>\n";

	# Projects published by others Table
	undef $list;
	undef $idxs;
	$list =  &getUserProjFromDB("other_published");
	print <<"HTML";  
		<div data-role='collapsible' data-collapsed-icon="carat-d" data-expanded-icon="carat-u" data-collapsed="true" >
		<h4>Projects published by others</h4>
HTML
;
	@theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Submission Time"),th("Total Running Time"),th("Owner"));
	$idxs = &sortList($list);
	$table_id = "edge-project-page-OtherPubtable";
	&printTable($table_id,$idxs,$list,\@theads);
	print "\t\t</div>\n";

}elsif ($umSystemStatus=~ /true/i) {
	my $list =  &getUserProjFromDB();
	my $idxs = &sortList($list);
	my @theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Submission Time"),th("Total Running Time"),th("Owner"));
	my $table_id = "edge-project-page-table";
	&printTable($table_id,$idxs,$list,\@theads);
}elsif ($um_config == 0) {
	# all projects in the EDGE_output
	my $list= &scanProjToList();
	my $idxs = &sortList($list);
	my @theads = (th("Project Name"),th({-style=>"width: 8%"},"Status"),th("Submission Time"),th("Total Running Time"),th("Last Run Time"));
	my $table_id = "edge-project-page-table";
	&printTable($table_id,$idxs,$list,\@theads);
}

print "\t\t\t</div>"; # data-filter

#print "<div role='main' id='edge-content-report' class='ui-content edge-content edge-output'>\n";

#print "<form class=\"ui-filterable\"><input id=\"edge-project-table-filter\" data-type=\"search\"> </form>\n";

#print $cgi->table( {-id=>"edge-project-page-table" , -class=>"tablesorter"},
sub sortList {
	my $list = shift;
	
	my @idxs1 = grep { $list->{$_}->{STATUS} eq "running" } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs2 = grep { $list->{$_}->{STATUS} ne "running" } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs = (@idxs1,@idxs2);
	return \@idxs;
}

sub printTable {
	my $table_id = shift;
	my $idx_ref = shift;
	my $list = shift;
	my $theads = shift;
	my @idxs = @{$idx_ref};
	my @tbodys;
	return if (@ARGV);
	if ($list->{INFO}->{ERROR})
	{
		print "<p class='error'>$list->{INFO}->{ERROR}</p>\n";
	}
	foreach (@idxs)
	{
		my $projOwner = $list->{$_}->{OWNER};
		my $projStatus = $list->{$_}->{PROJSTATUS};
		my $projname = "<a href=\"#\" class=\"edge-project-page-link \" title=\"$list->{$_}->{PROJDESC}\" data-pid=\"$list->{$_}->{PROJNAME}\">$list->{$_}->{REAL_PROJNAME}</a>";
		my $projSubTime = $list->{$_}->{PROJSUBTIME};
		my $projRunTime = $list->{$_}->{RUNTIME};
		my $projLastRunTime = $list->{$_}->{LASTRUNTIME};
		my $projType = $list->{$_}->{PROJ_TYPE};
		my $projAction;
		my $publish_action= ($projType =~ /published/)? "unpublished":"published";
		$projAction = "<a href='#edge_confirm_dialog' id='action-sblehare-btn1' aria-haspopup='true' data-rel='popup' data-position-to='window' datadata-transition='pop' class='edge-icon-bg-grey ui-icon-forward ui-btn ui-overlay-shadow ui-btn-icon-notext ui-corner-all' data='share'>Share project</a>\n"; 
		$projAction .= "<a href='#edge_confirm_dialog' id='action-share-btn2' aria-haspopup='true' data-rel='popup' data-position-to='window' datadata-transition='pop' class='edge-icon-bg-grey ui-icon-back ui-btn ui-overlay-shadow ui-btn-icon-notext ui-corner-all' data='share'>Unshare project</a>\n"; 
		$projAction .= "<a href='#edge_confirm_dialog' id='action-publish-btn' aria-haspopup='true' data-rel='popup' data-position-to='window' datadata-transition='pop' class='edge-icon-bg-grey ui-icon-eye ui-btn ui-overlay-shadow ui-btn-icon-notext ui-corner-all' data=\"$publish_action\">Make project public/private</a>\n"; 
		$projType =~ s/published/public/;
		my @tds;
		if ($umSystemStatus=~ /true/i){
			@tds = ($table_id =~ /Mytable/i)?  ( td($projname),td($projStatus),td($projSubTime),td($projRunTime),td($projType)):
							( td($projname),td($projStatus),td($projSubTime),td($projRunTime),td($projOwner));
		}else{
			@tds = ( td($projname),td($projStatus),td($projSubTime),td($projRunTime),td($projLastRunTime));
		}
		push @tbodys, \@tds;
	}
	if (scalar(@idxs)<1){
		my @tds = (td("No Projects"),td(""),td(""),td(""),td(""));
		push @tbodys, \@tds;
	}
	print $cgi->table( 
			{-id=>"$table_id" , -class=>"output-table ui-responsive ui-table ui-table-reflow" },
			thead(Tr(@{$theads})),
			tbody(
			map { Tr(@{$_}) } @tbodys
			)
	);
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

sub scanProjToList {
	my $cnt = 0;
	my $list;
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	while( defined (my $file = readdir BIN) ) {
		next if $file eq '.' or $file eq '..';
		if ( -d "$out_dir/$file" && -r "$out_dir/$file/process.log"  ) {
			++$cnt;
			$list=&pull_summary($file,$cnt,$list);
			$list->{$cnt}->{REAL_PROJNAME} = $list->{$cnt}->{PROJNAME};
			$list->{$cnt}->{PROJNAME} = $file;
		}
	}
	closedir(BIN);
	return $list;
}

sub getUserProjFromDB{
	my $project_type = shift;
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
	print Dumper ($array_ref) if @ARGV;
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		next if ($status =~ /delete/i);
		next if (! -r "$out_dir/$id/process.log");
		$list=&pull_summary($id,$id,$list);
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{REAL_PROJNAME} = $project_name;
		$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}";
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};
	}
	return $list;
}


sub pull_summary {
	my $pname = shift;
	my $cnt= shift;
	my $list = shift;
	my @INFILES;
	
	my ($step,$lastline);
	my $tol_running_sec=0;
	$list->{$cnt}->{PROJSTATUS} = "Unfinished";

	open(my $sumfh, "<", "$out_dir/$pname/process.log") or die $!;
	while(<$sumfh>) {
		chomp;
		#parse input files
		if( /runPipeline/ ) {
			undef @INFILES;	
		}
		if( /runPipeline .*-p (.*) -\w/ || /runPipeline .*-p (.*) >/ || /runPipeline .*-p (.*)$/) {
			push @INFILES, split /\s+/,$1;
		}
		if(/runPipeline .*-u (.*) -\w/ || /runPipeline .*-u (.*) >/ || /runPipeline .*-u (.*)$/){
			push @INFILES, split /\s+/,$1;
		}

		#parse reference files
		if( /runPipeline .*-ref (\S+)/){
			$list->{$cnt}->{REFFILE} = $1;
		}
		elsif( /^reference=(\S+)/ ){
			 $list->{$cnt}->{REFFILE} = $1;
		}

		if( /Total Running time: (\d+):(\d+):(\d+)/){
			$list->{$cnt}->{LASTRUNTIME} = "$1h $2m $3s";
			next;
		}
		elsif( /^Host=(.*)/ ){
			$list->{$cnt}->{HOSTFILE} = $1;
			next;
		}
		elsif( /^SNPdbName=(.*)/ ){
			$list->{$cnt}->{SPDB} = $1;
			next;
		}

		if( /^\[(.*)\]/ ){
			$step = $1;
			if( $step eq "system" ){
				while(<$sumfh>){
					chomp;
					if ( /^([^=]+)=([^=]+)/ ){
						$list->{$cnt}->{uc($1)}=$2;
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
			#$list->{$cnt}->{PROJSTATUS} = "Finished";
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
		}
		elsif( / Running/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "<span class='edge-fg-orange'>Running</span>";
			$list->{$cnt}->{PROJSTATUS} = "<span class='edge-fg-orange'>Running</span>";
		}
		elsif( /failed/ ){
			$list->{$cnt}->{$step}->{GNLSTATUS} = "<span class='edge-fg-red'>Failed</span>";
			$list->{$cnt}->{PROJSTATUS} = "<span class='edge-fg-red'>Failure</span>";
		}
		$lastline = $_;
	}

        #$list->{$cnt}->{RUNTIME} = strftime("\%H:\%M:\%S", gmtime($tol_running_sec));
	$list->{$cnt}->{RUNTIME} = sprintf("%02d:%02d:%02d", int($tol_running_sec / 3600), int(($tol_running_sec % 3600) / 60), int($tol_running_sec % 60));

	$list->{$cnt}->{PROJSTATUS}        = "Unstarted"   if $lastline =~ /EDGE_UI.*unstarted/;
	$list->{$cnt}->{PROJSTATUS}        = "Interrupted" if $lastline =~ /EDGE_UI.*interrupted/;
	$list->{$cnt}->{TIME}              = $1            if $lastline =~ /\[(\S+ \S+)\] EDGE_UI/;
	$list->{$cnt}->{$step}->{GNLSTATUS} = "Interrupted" if $list->{$cnt}->{$step}->{PROJSTATUS} eq "Interrupted"; #turn last step to unfinished
	
	#Reads Taxonomy Classification
	my %toolmap;
	open PROC_CUR, "<", "$out_dir/$pname/process_current.log" or die $!;
	while(<PROC_CUR>) {
		chomp;
		#parse input files
		if( /^\[RUN_TOOL\] \[(.*)\] COMMAND/ ){
			$step = $1;
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLANALYSIS} = "<span style='margin-left:3em'>$step</span>";
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLRUN}      = "On";
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLSTATUS}   = "<span class='edge-fg-orange'>Running</span>";
		}
		elsif( /^\[RUN_TOOL\] \[(.*)\] Error occured/ ){
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLSTATUS}   = "Error";
		}
		elsif( /^\[RUN_TOOL\] \[(.*)\] Running time: (.*)/ ){
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLSTATUS}   = "Complete";
			$list->{$cnt}->{"Reads Taxonomy Classification"}->{$step}->{GNLTIME}     = $2;
		}
		elsif( /All Done/){
			$list->{$cnt}->{PROJSTATUS} = "Complete";
		}
	}
	close PROC_CUR;


	$list->{$cnt}->{INFILES} = join ", ", @INFILES;
	$list->{$cnt}->{TIME} ||= strftime "%F %X", localtime;
	
	close ($sumfh);
	return $list;
}

sub check_um_service {
	my $url=shift;
	if (! LWP::Simple::head($url)) {
  	#	warn "The User managment Service is DOWN!!!! Will pull all projects from EDGE output direcotry";
  		return 0; 
	}else{
		return 1;
	}
}
