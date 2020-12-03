#!/usr/bin/env perl
use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use LWP::UserAgent;
use HTTP::Request::Common;
use JSON;
use File::Basename;
use POSIX qw(strftime);
use Data::Dumper;
use Email::Valid;
require "./edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $username    = $opt{'username'}|| $ARGV[0];
my $password    = $opt{'password'}|| $ARGV[1];
my $umSystemStatus    = $opt{'umSystem'}|| $ARGV[2];
my $userType;
my $viewType    = $opt{'view'}|| $ARGV[4];
my $protocol    = $opt{protocol}||'http:';
my $sid         = $opt{'sid'}|| $ARGV[5];
my $action = lc($opt{action});
my $domain      = $ENV{'HTTP_HOST'} || 'edge-prod.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;
my $info;
&stringSanitization(\%opt);
# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");
my $out_dir     = $sys->{edgeui_output};
my $report_dir     = $sys->{edgeui_report};
my $www_root	= $sys->{edgeui_wwwroot};
my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
(my $report_rel_dir = $report_dir) =~ s/$www_root//;
my $um_config	= $sys->{user_management};
my $um_url      = $sys->{edge_user_management_url};
$um_url ||= "$protocol//$domain/userManagement";
my $maintenance= ($sys->{"maintenance"})? $sys->{"maintenance"}:"0";

my $cluster     = $sys->{cluster};

# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if($valid){
		($username,$password) = getCredentialsFromSession($sid);
		my $user_info=&getUserInfo($username,$password);
		$userType=$user_info->{type};
	}
}

my $relt;
if($action eq "form") {
	print $cgi->header('application/json');
	my $reportFormTmpl    = "$RealBin/../../scripts/projects_report/edge_projects_report_form.tmpl";
	my $head_checkbox="<input type='checkbox' id='edge-reportform-ckall'>";
	my $projectListHtml;
	my $tableHtml = "<table id='edge-report-form-table' class='output-table ui-responsive ui-table ui-table-reflow'>";
	if ($username && $password){
		# My Table
		$tableHtml .= "<thead><tr>";
		$tableHtml .= "<th>$head_checkbox</th><th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th><th>Type</th><th>Owner</th>";
		$tableHtml .= " </tr></thead>";
		&getProjFromUM("user");

	} else {
		# all projects in the EDGE_output
		$tableHtml .= "<thead><tr>";
		$tableHtml .= "<th>$head_checkbox</th><th>Project Name</th><th>Status</th><th>Submission Time</th><th>Total Running Time</th>";
		$tableHtml .= " </tr></thead>";
		my $list= &scanProjToList();
		my $idxs = &sortList($list);
		&dtList($idxs,$list);
	}
	$tableHtml .= "</table>";

	my $html;
	open(my $fh, '<', $reportFormTmpl) or die "cannot open file $reportFormTmpl";
	{
		local $/;
		$html = <$fh>;
	}
	close($fh);

	$html =~ s/<PROJECT-LIST PLACEHOLDER>/$tableHtml/;
	$relt -> {html} = $html;
	my $relt_json = encode_json($relt);
	#print STDERR $relt_json;
	print $relt_json;
} elsif ($action eq "create") {
	my $html;
	my $report_name = $opt{'report'};
	if ($maintenance){
		print "Content-Type: text/html\n\n",
                       "<h3 class='error'>System is under maintenance. Please submit later or contact system administrator.</h3>";
		exit;
	}
	if ($username && $password){
		$report_name = &addReport2DB($opt{'edge-report-name'},$opt{'edge-report-desc'});
	}
	if(!$report_name || $info->{error}) {
		$html = $info->{error};
	} else {
		my $report_out_dir = "$report_dir/". $report_name;
		my $relative_outdir = "$report_rel_dir/". $report_name;
		system("mkdir -p $report_out_dir");

		my $reportSettings = "$report_out_dir/reportSettings.txt";
		open OUT,  ">$reportSettings"; 

		print OUT "report-name=".$opt{'edge-report-name'}."\n";
		print OUT "report-desc=".$opt{'edge-report-desc'}."\n";
		print OUT "run-name=".$opt{'checkbox-report-1-1'}."\n";
		print OUT "run-desc=".$opt{'checkbox-report-1-2'}."\n";
		print OUT "run-files=".$opt{'checkbox-report-1-3'}."\n";
		print OUT "sample-metadata=".$opt{'checkbox-report-1-4'}."\n";
		print OUT "sample-metadata-selected=".join(',',split(/[\x0]/,$opt{'metadata-list'}))."\n";
		print OUT "preprocess-stats=".$opt{'checkbox-report-2-1'}."\n";
		print OUT "preprocess-stats-selected=".join(',',split(/[\x0]/,$opt{'preprocess-stats-list'}))."\n";
		print OUT "preprocess-figures=".$opt{'checkbox-report-2-2'}."\n";
		print OUT "preprocess-figures-selected=".join(',',split(/[\x0]/,$opt{'preprocess-figure-list'}))."\n";
		print OUT "assembly-annotation-stats=".$opt{'checkbox-report-3-1'}."\n";
		print OUT "assembly-annotation-stats-selected=".join(',',split(/[\x0]/,$opt{'aa-stats-list'}))."\n";
		print OUT "assembly-annotation-figures=".$opt{'checkbox-report-3-2'}."\n";
		print OUT "assembly-annotation-figures-selected=".join(',',split(/[\x0]/,$opt{'aa-figure-list'}))."\n";
		print OUT "ref-stats=".$opt{'checkbox-report-4-1'}."\n";
		print OUT "ref-stats-selected=".join(',',split(/[\x0]/,$opt{'ref-stats-list'}))."\n";
		print OUT "ref-figures=".$opt{'checkbox-report-4-2'}."\n";
		print OUT "ref-figures-selected=".join(',',split(/[\x0]/,$opt{'ref-figure-list'}))."\n";
		print OUT "tax-stats=".$opt{'checkbox-report-5-1'}."\n";
		print OUT "tax-stats-selected=".join(',',split(/[\x0]/,$opt{'tax-stats-list'}))."\n";
		print OUT "tax-figures=".$opt{'checkbox-report-5-2'}."\n";
		print OUT "tax-figures-selected=".join(',',split(/[\x0]/,$opt{'tax-figure-list'}))."\n";
		print OUT "tax-tools=".$opt{'checkbox-report-5-3'}."\n";
		print OUT "tax-tools-selected=".join(',',split(/[\x0]/,$opt{'tax-tool-list'}))."\n";
		print OUT "tax-tool-pangia-score=".$opt{'tax-tool-pangia-score'}."\n";
		print OUT "projects-selected=".join(',',split(/[\x0]/,$opt{'edge-reportform-ckb'}))."\n";
		close OUT;
	
		my $reportTmpl    = "$RealBin/../../scripts/projects_report/edge_projects_report_html.tmpl";
		my $reportHtml =  "$report_out_dir/index.html";
		my $cmd = "cd $edgeui_wwwroot; $edgeui_wwwroot/../scripts/projects_report/create_report_w_temp.pl $report_out_dir $reportHtml $reportSettings";
		#print STDERR "$cmd\n";
		`$cmd`;

		my $pr=0;
		my @htmls;
		open REP, "$reportHtml" or die "Can't open $reportHtml: $!";
		foreach(<REP>){
			last if /<!-- \/content -->/;
			push @htmls, $_ if $pr;
			$pr=1 if /id='edge-projects-report-page'/;
		}

		close REP;

		$html = join "", @htmls;
	}
	print "Content-Type: text/html\n\n",
				  $html;
} elsif ($action eq "list") {
	my $list;
	if ($username && $password){
		$list = &getReportsFromDB();
	} else {
		$list = &scanReportsToList();
	}
	my $idxs = &sortReportList($list);
	my $head_checkbox="<input type='checkbox' id='edge-reportlistpage-ckall'>";
	my @theads = (th("$head_checkbox"),th("Report Name"),th("Description"),th("Owner"),th("Created Time"));
	my $table_id = "edge-reports-page-table";
	my $html = "<h2>Report List</h2><div id='edge-reportlist-action' class='flex-container'>\n";
	$html .= '<a href="" title="Delete Selected Reports" class="tooltip ui-btn ui-btn-d ui-icon-delete ui-btn-icon-notext ui-corner-all" data-role="button" role="button">delete</a>';
	if ($um_config != 0){
		$html .= '<a href="" title="Share Selected Reports" class="tooltip ui-btn ui-btn-d ui-icon-forward ui-btn-icon-notext ui-corner-all" data-role="button" role="button">share</a>';
		$html .= '<a href="" title="Unshare Selected Reports" class="tooltip ui-btn ui-btn-d ui-icon-forbidden ui-btn-icon-notext ui-corner-all" data-role="button" role="button">unshare</a>';
 	}
	$html .= "</div>\n";
	$html .= &printReportTable($table_id,$idxs,$list,\@theads);
	print "Content-Type: text/html\n\n",
				  $html;
} elsif ($action eq "delete") {
	my $report = $opt{'report'};
	$report =~ s/<li>//;
	$report =~ s/<\/li>//;
	$report =~ s/^\s+|\s+$//g;
	my $report_id = $opt{'report_id'};
	my $report_db_id;
	if ($username && $password){
		$report_db_id= $report_id;		
		my $list = &getReportsFromDB("owner");
		if( ! $list->{$report_db_id}) {
			$info->{STATUS} = "FAILURE";
			$info->{INFO} = "ERROR: Permission denied. Only report owner can delete report '$report'.<br>";
			&returnStatus();
		} else {
			$report_id = $list->{$report_db_id}{REPORTNAME};
		}
	}
	my $report_out_dir = "$report_dir/". $report_id;
	if( -d $report_out_dir ){
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Failed to delete the report $report output directory.";
		`rm -rf $report_out_dir`;
		if( !-e $report_out_dir){
			$info->{STATUS} = "SUCCESS";
			$info->{INFO}   = "Report $report has been deleted.";
			if ($username && $password){
				&deleteReportFromDB($report_db_id);
			}
		}
	} else{
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Report $report output directory not found.";
	}
	&returnStatus();
} elsif ($action eq "share" || $action eq 'unshare') {
	my $report = $opt{'report'};
	$report =~ s/<li>//;
	$report =~ s/<\/li>//;
	$report =~ s/^\s+|\s+$//g;
	my $report_id = $opt{'report_id'};
	if ($username && $password){		
		my $list = &getReportsFromDB("owner");
		if( ! $list->{$report_id}) {
			$info->{STATUS} = "FAILURE";
			$info->{INFO} = "ERROR: Permission denied. Only report owner can $action report '$report'.<br>";
			&returnStatus();
		} else {
			my $shareEmail = $opt{shareEmail};
			&shareReport($report,$report_id,$shareEmail,$action);
		}
	} else{
		$info->{STATUS} = "FAILURE";
		$info->{INFO}   = "Please login to do this action.";
	}
	&returnStatus();
}
## END MAIN## 

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
		
		my $projname = "<a href='#' class='edge-report-form-link ' title='alt-click to open in a new tab' data-pid='$id'>$name</a>";
		my $checkbox = "<input type='checkbox' class='edge-reportform-ckb' name='edge-reportform-ckb' value='$code'>";

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

sub sortReportList {
	my $list = shift;
	
	my @idxs = sort {$list->{$b}->{REPORTTIME} cmp $list->{$a}->{REPORTTIME}} keys %$list;
	return \@idxs;
}

sub printReportTable {
	my $reportListHtml = '';
	my $table_id = shift;
	my $idx_ref = shift;
	my $list = shift;
	my $theads = shift;
	my @idxs = @{$idx_ref};
	my @tbodys;
	#return if (@ARGV);
	if ($list->{INFO}->{ERROR})
	{
		$reportListHtml = "<p class='error'>$list->{INFO}->{ERROR}</p>\n";
	}
	foreach (@idxs)
	{
		my $reportOwner = $list->{$_}->{OWNER};
		my $reportStatus = $list->{$_}->{REPORTSTATUS};
		my $reportDesc = $list->{$_}->{REPORTDESC};
		my $reportID = $list->{$_}->{REPORTID};
		my $reportname = "<a href='$report_rel_dir/$list->{$_}->{REPORTNAME}' class='edge-reportpage-link' data-rid='$reportID' target='_blank'>$list->{$_}->{REAL_REPORTNAME}</a>";
		my $reportTime = $list->{$_}->{REPORTTIME};
		my $reportCode = $list->{$_}->{REPORTCODE} || $list->{$_}->{REAL_REPORTNAME};
		my $checkbox = "<input type='checkbox' class='edge-reportlist-ckb' name='edge-reportlist-ckb' value=\'$reportCode\'>";
		
		my @tds;
		@tds = ( td($checkbox),td($reportname),td($reportDesc),td($reportOwner),td($reportTime));
		
		push @tbodys, \@tds;
	}

	if (scalar(@idxs)<1){
		my @tds = (td(""),td("No Reports"),td(""),td(""),td(""));
		
		push @tbodys, \@tds;
	}
	$reportListHtml = $cgi->table( 
			{-id=>"$table_id" , -class=>"output-table ui-responsive ui-table ui-table-reflow" },
			thead(Tr(@{$theads})),
			tbody(
			map { Tr(@{$_}) } @tbodys
			)
	);

	return $reportListHtml;
}

sub sortList {
	my $list = shift;
	
	my @idxs1 = grep { $list->{$_}->{PROJSTATUS} =~ /running/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs2 = grep { $list->{$_}->{PROJSTATUS} =~ /unstarted/i } sort {$list->{$a}->{REAL_PROJNAME} cmp $list->{$b}->{REAL_PROJNAME}} keys %$list;
	my @idxs3 = grep { $list->{$_}->{PROJSTATUS} !~ /running|unstarted/i } sort {$list->{$b}->{TIME} cmp $list->{$a}->{TIME}} keys %$list;
	my @idxs = (@idxs1,@idxs2,@idxs3);
	return \@idxs;
}

sub dtList {
	my $idx_ref = shift;
	my $list = shift;
	my @idxs = @{$idx_ref};
	
	if ($list->{INFO}->{ERROR})
	{
		return  "<p class='error'>$list->{INFO}->{ERROR}</p>\n";
	}
	foreach (@idxs)
	{
		my $projOwner = $list->{$_}->{OWNER};
		my $projStatus = $list->{$_}->{PROJSTATUS};
		my $projID = $list->{$_}->{PROJNAME};
		my $projname = "<a href='#' class='edge-report-form-link ' title='$list->{$_}->{PROJDESC} (alt-click to open in a new tab)' data-pid='$projID'>$list->{$_}->{REAL_PROJNAME}</a>";
		my $projSubTime = $list->{$_}->{PROJSUBTIME};
		my $projRunTime = $list->{$_}->{RUNTIME};
		my $projType = $list->{$_}->{PROJ_TYPE};
		my $projCode = $list->{$_}->{PROJCODE} || $list->{$_}->{REAL_PROJNAME};
		my $checkbox = "<input type='checkbox' class='edge-reportform-ckb' name='edge-reportform-ckb' value='$projCode'>";
		my $publish_action= ($projType =~ /published/)? "unpublished":"published";
		$projType =~ s/published/public/;
		my @tds;
		if ($umSystemStatus=~ /true/i){
			@tds = ($checkbox,$projname,$projStatus,$projSubTime,$projRunTime,$projType,$projOwner);
			
		}else{
			@tds = ($checkbox,$projname,$projStatus,$projSubTime,$projRunTime);
		}
		push @{$relt->{data}}, [@tds];
	}
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

sub scanProjToList {
	my $cnt = 0;
	my $list;
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	while( defined (my $file = readdir BIN) ) {
		next if $file eq '.' or $file eq '..' or $file eq 'sample_metadata_export';
		if ( -d "$out_dir/$file"  && -r "$out_dir/$file/config.txt") {
			++$cnt;
			if ( -e "$out_dir/$file/.AllDone" && -e "$out_dir/$file/HTML_Report/writeHTMLReport.finished"){
				$list=&get_start_run_time("$out_dir/$file/.AllDone",$cnt,$list);
				$list->{$cnt}->{PROJSTATUS} = "Complete";
			}else
			{
				if (-r "$out_dir/$file/process.log" ){
                                        $list=&pull_summary("$out_dir/$file/process.log",$cnt,$list,"$out_dir/$file")
                                }else{
                                        $list->{$cnt}->{PROJSTATUS} = "Unstarted";
                                }
				$list=&pull_summary("$out_dir/$file/config.txt",$cnt,$list,"$out_dir/$file") if ($list->{$cnt}->{PROJSTATUS} =~ /unstart/i);
			}
			$list->{$cnt}->{REAL_PROJNAME} = $file unless $list->{$cnt}->{REAL_PROJNAME};
			$list->{$cnt}->{PROJNAME} = $file;
		}
	}
	closedir(BIN);
	return $list;
}

sub getUserProjFromDB{
	my $project_type = shift;
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
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $projCode = $hash_ref->{code};
		my $project_name = $hash_ref->{name};
		my $status = $hash_ref->{status};
		my $created_time = $hash_ref->{created};
		next if ($status =~ /delete/i);
		next if (! -r "$out_dir/$id/process.log" && ! -r "$out_dir/$projCode/process.log" && ! $cluster);
		my $proj_dir=(-d "$out_dir/$projCode")?"$out_dir/$projCode":"$out_dir/$id";
		my $processlog = (-r "$proj_dir/process.log")? "$proj_dir/process.log":"$proj_dir/config.txt";
		next if(! -e $processlog);
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
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{PROJSTATUS} = $status if (!$list->{$id}->{PROJSTATUS});
		$list->{$id}->{REAL_PROJNAME} = $project_name if (!$list->{$id}->{REAL_PROJNAME});
		$list->{$id}->{PROJCODE} = $projCode;
		$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}";
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};
	}
	return $list;
}

sub get_start_run_time{
	my $log = shift;
	my $id = shift;
	my $list = shift;
	my $tol_running_sec=0;
	my ($start_time, $run_time);
	open (my $log_fh, $log) or die "Cannot read $log\n";
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
        if(! -e $log) {return;}
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
		`echo "$list->{$cnt}->{PROJSUBTIME}\t$run_time" > $proj_dir/.AllDone` if ($list->{$cnt}->{RUNTIME});
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

sub getSettings {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		chomp;
		if ( /^([^=]+)=([^=]+)/ ){
			$sys->{$1}=$2;
		}
	}
	close CONF;
	return $sys;
}

sub scanReportsToList {
	my $cnt = 0;
	my $list;
	my $settings_file;
	opendir(BIN, $report_dir) or die "Can't open $report_dir: $!";
	while( defined (my $file = readdir BIN) ) {
		next if $file eq '.' or $file eq '..';
		$settings_file = "$report_dir/$file/reportSettings.txt";
		if ( -d "$report_dir/$file"  && -r $settings_file) {
			++$cnt;
			
			my $settings = &getSettings($settings_file);
			if ( $settings -> {'Created Time'}){
				$list->{$cnt}->{REPORTSTATUS} = "Complete";
				$list->{$cnt}->{REPORTTIME} = $settings -> {'Created Time'};
			} else {
				$list->{$cnt}->{REPORTSTATUS} = "In Process";
			}
			$list->{$cnt}->{REPORTDESC} = $settings -> {'report-desc'};
			$list->{$cnt}->{REAL_REPORTNAME} = $settings -> {'report-name'} || $file;
			$list->{$cnt}->{REPORTNAME} = $file;
			$list->{$cnt}->{REPORTID} = $file;
			$list->{$cnt}->{OWNER} = "All";
		}
	}
	closedir(BIN);
	return $list;
}

sub returnStatus {
	my $json = "{}";
	$json = to_json($info) if $info;
	$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
	print $cgi->header('application/json') unless $ARGV[0];
	print $json;
	exit;
}

sub addReport2DB{
	my $name = shift; 
	my $desc = shift;
	$desc =~ s/(['"])/\\$1/g;

	my %data = (
		email => $username,
   		password => $password,
   		report_name => $name,
   		description => $desc,
	);
	# Encode the data structure to JSON
	my $data = to_json(\%data);
	#w Set the request parameters
	my $url = $um_url ."WS/report/add";
	my $browser = LWP::UserAgent->new;
	my $req = POST $url;
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	if ($result_json =~ /\"error_msg\":"(.*)"/){
		$info->{error}=$1;
		return;
    	}
	my $info =  from_json($result_json);
	my $new_name =  $info->{"code"};
	return $new_name;
}
sub getUserInfo {
	my $username=shift;
	my $password=shift;
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
sub getReportsFromDB{
	my $type = shift;
	my $list;
	my $settings_file;
        my %data = (
                email => $username,
                password => $password
        );
	if($type) {
		$data{report_type} = $type;
	}
        # Encode the data structure to JSON
        #w Set the request parameters
	my $service = "WS/user/getReports";
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
                $info->{error}=$1;
                return;
        }
        my $array_ref =  decode_json($result_json);
	#print Dumper ($array_ref) if @ARGV;
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $reportCode = $hash_ref->{code};
		$settings_file = "$report_dir/$reportCode/reportSettings.txt";
		if ( -d "$report_dir/$reportCode"  && -r $settings_file) {
			$list->{$id}->{REPORTTIME} = $hash_ref->{created};
			$list->{$id}->{REPORTDESC} = $hash_ref->{description};
			$list->{$id}->{REAL_REPORTNAME} = $hash_ref->{name};
			$list->{$id}->{REPORTNAME} = $reportCode;
			$list->{$id}->{REPORTID} = $id;
			$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}";
			$list->{$id}->{REPORTSTATUS} = "Complete";
		}
	}
	return $list;
}

sub deleteReportFromDB {
	my $id = shift;
	my %data = (
		email => $username,
		password => $password,
		report_id => $id
	);
	# Encode the data structure to JSON		
	my $data =  encode_json(\%data);

	# Set the request parameters
	my $url = $um_url. 'WS/report/delete';
	my $browser = LWP::UserAgent->new;
	my $req = PUT $url;
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
	my $info =  from_json($result_json);
	return  $info->{"code"};
}

sub shareReport{
	my $report_name=shift;
	my $report_id=shift;
	my $email=shift;
	my $action =shift;
	$email =~ s/ //g;
	# avoid share to owner self.
	$email = join(',', grep (!/$username/, split(',',$email)));
	my %data = (
                email => $username,
                password => $password,
                report_id => $report_id,
        );
	$data{"${action}_to"} = $email;
	 # Encode the data structure to JSON
        my $data = to_json(\%data);
        #w Set the request parameters
        my $url = $um_url ."WS/report/${action}_user";
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
                $info->{INFO} .= " Report $report_name has ${action}d to $email.";
        }else{
		$action =~ s/^([a-z])/\u$1/;
                $info->{STATUS} = "FAILURE";
                $info->{INFO} .= " $action Report $report_name to $email failed. $result->{error_msg}";
	}
}

sub stringSanitization{
	my $opt=shift;
	my $dirtybit=0;
	foreach my $key (keys %$opt){
		my $str = $opt->{$key};
		next if $key eq "password";
		next if $key eq "keywords";

		if ($key eq "username" or $key eq "shareEmail"){
			my @email = split(',',$str);
			map { $dirtybit =1  if ! &emailValidate($_); } @email;
		}else{
			$opt->{$key} =~ s/[`;'"]/ /g;
			$str =~ s/[`;'"]/ /g;
			$str =~ s/[\x0]//g;
			next if $key eq "edge-report-desc";
			$dirtybit=1 if ($str =~ /[^0-9a-zA-Z\_\@\-\=\:\\\.\/\+ ]/);
		}
		if ($dirtybit){
			$info->{INFO} = "Invalid characters detected \'$str\'.";
			&returnStatus();
		}
	}
}

sub emailValidate{
	my $email=shift;
	$email = lc($email);
	return Email::Valid->address($email);
}


