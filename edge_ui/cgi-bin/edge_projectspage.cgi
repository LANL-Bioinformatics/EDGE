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
use Email::Valid;
#use CGI::Pretty;
require "edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $username    = $opt{'username'}|| $ARGV[0];
my $password    = $opt{'password'}|| $ARGV[1];
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[2];
my $userType;
my $viewType    = $opt{'view'}|| $ARGV[4];
my $protocol    = $opt{protocol}||'http:';
my $sid         = $opt{'sid'}|| $ARGV[5];
my $loadNum    = $opt{'loadNum'} || 100;
my $domain      = $ENV{'HTTP_HOST'} || 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
my $out_dir     = $sys->{edgeui_output};
my $um_config	= $sys->{user_management};
my $um_url      = $sys->{edge_user_management_url};
$um_url ||= "$protocol//$domain/userManagement";

my $cluster     = $sys->{cluster};
# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid){
		($username,$password) = getCredentialsFromSession($sid);
	}
	my $user_info=&getUserInfo();
        $userType=$user_info->{type};
}

#print Dumper ($list);
print  $cgi->header( "text/html" );
print  $cgi->h2("Project List");
if ( $username && $password || $um_config == 0){
	#Action buttons
	print "<div id='edge-projpage-action' class='flex-container'>\n";
	if ($userType =~ /admin/i){
		print '<a href="" title="See All Projects List (admin)" class="tooltip ui-btn ui-btn-d ui-icon-bars ui-btn-icon-notext ui-corner-all" data-role="button" role="button">show-all</a>';
	}
	print '<a href="" title="Force Selected Projects to rerun" class="tooltip ui-btn ui-btn-d ui-shadow-icon ui-icon-refresh ui-btn-icon-notext ui-corner-all" data-role="button" >rerun</a>';
	print '<a href="" title="Interrupt running Projects" class="tooltip ui-btn ui-btn-d ui-icon-forbidden ui-btn-icon-notext ui-corner-all" data-role="button" role="button">interrupt</a>';
	print '<a href="" title="Delete Selected Projects" class="tooltip ui-btn ui-btn-d ui-icon-delete ui-btn-icon-notext ui-corner-all" data-role="button" role="button">delete</a>';
	print '<a href="" title="Empty Selected Projects Output" class="tooltip ui-btn ui-btn-d ui-icon-recycle ui-btn-icon-notext ui-corner-all" data-role="button" role="button">empty</a>';
	if ($sys->{edgeui_archive}){
		print '<a href="" title="Archive Selected Projects" class="tooltip ui-btn ui-btn-d ui-icon-arrow-u-r ui-btn-icon-notext ui-corner-all" data-role="button" role="button">archive</a>';
 	}
	 if ($um_config != 0){
		print '<a href="" title="Share Selected Projects" class="tooltip ui-btn ui-btn-d ui-icon-forward ui-btn-icon-notext ui-corner-all" data-role="button" role="button">share</a>';
		print '<a href="" title="Make Selected Projects Public" class="tooltip ui-btn ui-btn-d ui-icon-eye ui-btn-icon-notext ui-corner-all" data-role="button" role="button">publish</a>';
		print '<a href="" title="Disable Selected Projects Display" class="tooltip ui-btn ui-btn-d ui-icon-minus ui-btn-icon-notext ui-corner-all" data-role="button" role="button">disable-project-display</a>';
		print '<a href="" title="Enable Selected Projects Display" class="tooltip ui-btn ui-btn-d ui-icon-plus ui-btn-icon-notext ui-corner-all" data-role="button" role="button">enable-project-display</a>';
 	}
	print '<a href="" title="Compare Selected Projects Taxonomy Classification (HeatMap)" class="tooltip ui-btn ui-btn-d ui-icon-bullets ui-btn-icon-notext ui-corner-all" data-role="button" role="button">compare</a>';
 	if($sys->{edge_sample_metadata}) {
 		print '<a href="" title="Export Selected Projects Metadata" class="tooltip ui-btn ui-btn-d ui-icon-arrow-d ui-btn-icon-notext ui-corner-all" data-role="button" role="button">metadata-export</a>';
 	}
 	if($sys->{edge_sample_metadata} && $sys->{edge_sample_metadata_share2bsve}) {
 		print '<a href="" title="Share Selected Projects Metadata/Pathogens with BSVE" class="tooltip ui-btn ui-btn-d ui-icon-arrow-u ui-btn-icon-notext ui-corner-all" data-role="button" role="button">metadata-bsveadd</a>';
 	}
 	print '</div>';
}

#print "<div data-filter='true' id='edge-project-list-filter' data-filter-placeholder='Search projects ...'> \n";
#print "<form id='edge-projpage-form'>\n";
my $head_checkbox="<input type='checkbox' id='edge-projpage-ckall'>";
if ($umSystemStatus=~ /true/i && $username && $password && $viewType =~ /user/i ){
	# My Table
	my $list = &getUserProjFromDB("owner",$loadNum);
	my $list_g = &getUserProjFromDB("guest",$loadNum);
	my $list_p = &getUserProjFromDB("other_published",$loadNum);
	$list = &ref_merger($list, $list_p) if $list_p;
	$list = &ref_merger($list, $list_g) if $list_g;
	#<div data-role='collapsible-set' id='edge-project-list-collapsibleset'> 

	my @theads = (th("$head_checkbox"),th("Project Name"),th("Status"),th("Display"),th("Submission Time"),th("Total Running Time"),th("Type"),th("Owner"));
	my ($idxs,$return_idxs)= &sortList($list,$loadNum);
	my $table_id = "edge-project-page-table";
	&printTable($table_id,$return_idxs,$list,\@theads);

}elsif ($umSystemStatus=~ /true/i) {
	# show admin list or published project
	$head_checkbox="" if (!$username && !$password);
	my $list =  &getUserProjFromDB("",$loadNum);
	my ($idxs,$return_idxs)= &sortList($list,$loadNum);
	my @theads = (th("$head_checkbox"),th("Project Name"),th("Status"),th("Submission Time"),th("Total Running Time"),th("Owner"));
	my $table_id = "edge-project-page-table";
	&printTable($table_id,$return_idxs,$list,\@theads);
}elsif ($um_config == 0) {
	# all projects in the EDGE_output
	my $list= &scanProjToList($loadNum);
	my ($idxs,$return_idxs)= &sortList($list,$loadNum);
	my @theads = (th("$head_checkbox"),th("Project Name"),th("Status"),th("Submission Time"),th("Total Running Time"),th("Last Run Time"));
	my $table_id = "edge-project-page-table";
	&printTable($table_id,$return_idxs,$list,\@theads);
}



## END MAIN## 

sub sortList {
	my $list = shift;
	my $loadNum = shift;
	
	my @idxs1 = grep { $list->{$_}->{PROJSTATUS} =~ /running/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs2 = grep { $list->{$_}->{PROJSTATUS} =~ /unstarted/i } sort {$list->{$a}->{REAL_PROJNAME} cmp $list->{$b}->{REAL_PROJNAME}} keys %$list;
	my @idxs3 = grep { $list->{$_}->{PROJSTATUS} !~ /running|unstarted/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs = (@idxs1,@idxs2,@idxs3);
	my $totalNum=scalar(@idxs);
	my @return_idxs = @idxs;
	my $load_num_button="<div>Load recent <input type='number' id='edge-projpage-loadnum' data-role='none' min='100' max=\"$totalNum\" value=\'$loadNum\' step='10'> /$totalNum projects. <input type='button' data-role='none' id='edge-projpage-loadnum-submit' value='Go'></div>";
	print $load_num_button if ($totalNum>100);
	if ( $totalNum > $loadNum) {
		@return_idxs=@idxs[0..$loadNum];
	}
	#print join ("<br/>",@return_idxs,$loadNum);
	return (\@idxs,\@return_idxs);
}

sub printTable {
	my $table_id = shift;
	my $idx_ref = shift;
	my $list = shift;
	my $theads = shift;
	my @idxs = @{$idx_ref};
	my @tbodys;
	#return if (@ARGV);
	if ($list->{INFO}->{ERROR})
	{
		print "<p class='error'>$list->{INFO}->{ERROR}</p>\n";
	}
	foreach (@idxs)
	{
		my $projOwner = $list->{$_}->{OWNER};
		my $projStatus = $list->{$_}->{PROJSTATUS};
		my $projDisplay = $list->{$_}->{PROJDISPLAY};
		my $projID = $list->{$_}->{PROJNAME};
		my $projname = "<a href=\"#\" class=\"edge-project-page-link \" title=\"$list->{$_}->{PROJDESC} (alt-click to open in a new tab)\" data-pid=\"$projID\">$list->{$_}->{REAL_PROJNAME}</a>";
		my $projSubTime = $list->{$_}->{PROJSUBTIME};
		my $projRunTime = $list->{$_}->{RUNTIME};
		my $projLastRunTime = $list->{$_}->{LASTRUNTIME};
		my $projType = $list->{$_}->{PROJ_TYPE};
		my $projCode = $list->{$_}->{PROJCODE} || $list->{$_}->{REAL_PROJNAME};
		my $checkbox = "<input type='checkbox' class='edge-projpage-ckb' name='edge-projpage-ckb' value=\'$projCode\'>";
		my $publish_action= ($projType =~ /published/)? "unpublished":"published";
		$projType =~ s/published/public/;
		my @tds;
		if ($umSystemStatus=~ /true/i){
			$checkbox="" if (!$username && !$password);
			if( scalar @$theads == 8 ){
				@tds = ( td($checkbox),td($projname),td($projStatus),td($projDisplay),td($projSubTime),td($projRunTime),td($projType),td($projOwner) );
			}
			else{
				@tds = ( td($checkbox), td($projname),td($projStatus),td($projSubTime),td($projRunTime),td($projOwner) );
			}
		}else{
			@tds = ( td($checkbox),td($projname),td($projStatus),td($projSubTime),td($projRunTime),td($projLastRunTime));
		}
		push @tbodys, \@tds;
	}

	if (scalar(@idxs)<1){
		my @tds = (td(""),td("No Projects"),td(""),td(""),td(""),td(""));
		if( scalar @$theads == 8 ){
			@tds = (td(""),td("No Projects"),td(""),td(""),td(""),td(""),td(""),td(""));
		}

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

sub scanProjToList {
	my $loadNum=shift;
	my $cnt = 0;
	my $list;
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	my @dirs= sort { -M $a <=> -M $b } readdir BIN;
	foreach my $file (@dirs) {
		next if $file eq '.' or $file eq '..' or $file eq 'sample_metadata_export';
		if ( -d "$out_dir/$file"  && -r "$out_dir/$file/config.txt") {
			++$cnt;
			if ($cnt <= $loadNum){
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
			}
			$list->{$cnt}->{REAL_PROJNAME} = $list->{$cnt}->{PROJNAME} || $file;
			$list->{$cnt}->{PROJNAME} = $file;
		}
	}
	closedir(BIN);
	return $list;
}
sub getUserInfo {
	my $service = "WS/user/getInfo";
	my $url = $um_url .$service;
	my %user_info;
	my %data = (
		email => $username,
		password => $password
 	);
	my $data = to_json(\%data);
	my $browser = LWP::UserAgent->new;
	my $req = PUT $url;
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

 	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	print $result_json,"\n" if (@ARGV);
	if ($result_json =~ /\"error_msg\":"(.*)"/){
		print "Content-Type: text/html\n\n", "Fail to get user info\n\n";
		exit;
 	}
	my $tmp_r = from_json($result_json);
	$user_info{$_} = $tmp_r->{$_} foreach (keys %$tmp_r);
	return \%user_info;
}

sub getUserProjFromDB{
	my $project_type = shift;
	my $loadNum = shift;
	my $numProj=0;
	my $list={};
        my %data = (
                email => $username,
                password => $password
        );
        # Encode the data structure to JSON
        #w Set the request parameters
	my $service;
	if ($username && $password){ 
		$service= ($viewType =~ /admin/i)? "WS/user/admin/getProjects" :"WS/user/getProjects";
		$data{project_type} = $project_type if ($viewType =~ /user/i && $project_type);
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
        my $array_ref =  decode_json($result_json);
	#print Dumper ($array_ref) if @ARGV;
	my @projectlist=@$array_ref;
	foreach my $hash_ref ( sort {$b->{created} cmp $a->{created}} @projectlist)
	{
		my $id = $hash_ref->{id};
		my $projCode = $hash_ref->{code};
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		my $created_time = $hash_ref->{created};
		my $display = $hash_ref->{display}||"yes";
		if($project_type eq "other_published") {
			$display = "no";
		}
		next if ($status =~ /delete/i);
		next if (! -r "$out_dir/$id/process.log" && ! -r "$out_dir/$projCode/process.log" && ! $cluster);
		my $proj_dir=(-d "$out_dir/$projCode")?"$out_dir/$projCode":"$out_dir/$id";
		my $processlog = (-r "$proj_dir/process.log")? "$proj_dir/process.log":"$proj_dir/config.txt";
		if ($numProj < $loadNum){
		#		$list=&pull_summary($processlog,$id,$list,$proj_dir) if (  -e $processlog);
		#		&updateDBProjectStatus("$id","finished") if ( $list->{$id}->{PROJSTATUS} =~ /Complete/i);
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
                        	$list->{$id}->{RUNTIME}="";
                        	$status="In process";
			}else{
				$list=&pull_summary($processlog,$id,$list,$proj_dir) if (  -e $processlog);
			}
		}
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{PROJSTATUS} = $status if (!$list->{$id}->{PROJSTATUS});
		$list->{$id}->{PROJDISPLAY} = $display;
		$list->{$id}->{REAL_PROJNAME} = $project_name if (!$list->{$id}->{REAL_PROJNAME});
		$list->{$id}->{PROJCODE} = $projCode;
		$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}";
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};

		$numProj++;
	}
	return $list;
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
	my $data = to_json(\%data);

        my $url = $um_url ."WS/project/update";
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
                #$info->{INFO} .= " Update Project status in database failed.";
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

sub check_um_service {
	my $url=shift;
	if (! LWP::Simple::head($url)) {
  	#	warn "The User managment Service is DOWN!!!! Will pull all projects from EDGE output direcotry";
  		return 0; 
	}else{
		return 1;
	}
}

sub ref_merger {
	my ($r1, $r2) = @_;
	foreach my $key (keys %$r2){
		if($r1->{$key}) {
			$r2->{$key}->{PROJ_TYPE} .= ",public";
		}
		$r1->{$key} = $r2->{$key};
	}
	return $r1;
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
		}
	}
}

sub emailValidate{
	my $email=shift;
	$email = lc($email);
	return Email::Valid->address($email);
}

