#!/usr/bin/env perl

use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use File::Basename;
use Digest::MD5 qw(md5_hex);
use POSIX qw(strftime);
require "./edge_user_session.cgi";
require "../metadata_scripts/metadata_api.pl";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
my $sid         = $opt{'sid'};
my $userType;
my $prealname = $opt{projname};
my $action = $opt{'action'};
my $userDir = $opt{'userDir'};
my $msg;

######################################################################################
# DATA STRUCTURE:
#	$msg->{"PARAMS"}->{0}->{metadata-virus-name}  =
#			     ->{metadata-sample-host} =
#			     ...
# 	                ->{1}...

#	$msg->{SUBMISSION_STATUS}
#	$msg->{PATH}   # file for action download 
     
#######################################################################################

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $um_url      = $sys->{edge_user_management_url};
my ($username,$password);
my $domain      = $ENV{'HTTP_HOST'} || 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");

my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};
my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;
my $ownProjlist;
my ($projCode,$projStatus,$projOwnerEmail);
# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if(!$valid){
		my $html = "Invalid session. Please login again!";
		print "Content-Type: text/html\n\n",
		$html;
		exit;
	}else{
		($username,$password) = getCredentialsFromSession($sid);
	}
	$ownProjlist = &getUserProjFromDB("owner");
	($prealname,$projCode,$projStatus,$projOwnerEmail) = &getProjNameFromDB($pname) if ($action !~ /batch-/);
}else{
	($prealname,$projCode,$projStatus)= &scanProjToList($edgeui_output,$pname) if ($action !~ /batch-/);
}


if($action eq "create-form") {
	my $proj_abs_dir = $edgeui_output . "/". $projCode;
	if ( ! -d "$proj_abs_dir"){
		my $html = "Cannot find Project $prealname direcotry!";
		print "Content-Type: text/html\n\n",
		$html;
		exit;
	}
	my $projDir = $relpath . "/". $projCode;
	chdir $edgeui_wwwroot;
	my $outHtml = "$projDir/GISAID/upload.html";
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/metadata/outputGisaidMetadata_w_temp.pl $projDir $outHtml $prealname $profileDir";
	`$cmd`;
	#print STDERR "$cmd";

	open REP, "$outHtml" or die "Can't open $outHtml: $!";
	my $pr=0;
	my @htmls;
	while(<REP>){
		last if /<!-- \/content -->/;
		push @htmls, $_ if $pr;
		$pr=1 if /id='edge-content-report'/;
	}

	close REP;

	my $html = join "", @htmls;
	#$html =~ s/>\s+</></mg;
	#$html =~ s/\t/ /mg;
	#$html =~ s/ {2,}//mg;
	#$html =~ s/\n/ /mg;

	print "Content-Type: text/html\n\n",
		  $html;

		  #open LOG, ">/tmp/edge_report$$.log";
		  #print LOG $html;
		  #close LOG;
}elsif( $action eq 'create-batch-form'){
	my @projCodes = split /,/,$opt{proj};
	my $short_username = $1 if $username =~ /(.*)\@/;
        $short_username =~ s/\W/_/;
	my $metadata_out_dir = "$edgeui_output/sample_metadata_export/$short_username/". md5_hex(join ('',@projCodes));
	chdir $edgeui_wwwroot;
	my $outHtml = "$metadata_out_dir/batch_metadata.html";
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $projects = join(",",map { "$edgeui_output/$_" } @projCodes);

	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/metadata/outputGisaidMetadata_w_temp_table.pl -projects $projects -out $outHtml -udir $profileDir";
	`$cmd`;
	#print STDERR "$cmd";

	open REP, "$outHtml" or die "Can't open $outHtml: $!";
	my $html_content;
	while(<REP>){
		$html_content .= $_;
	}
	print "Content-Type: text/html\n\n",
		$html_content;
}elsif($action eq "upload2gisaid" or $action eq "download" or $action eq "update") {
	$msg->{SUBMISSION_STATUS}="success";

	# validating parameters
	&checkParams();
	if($msg->{SUBMISSION_STATUS} eq "failure") {
		returnParamsStatus();
	}

	my $projDir = $edgeui_output . "/". $projCode;
	my $projDir_rel = $relpath . "/". $projCode;
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	if(-e $projDir) {
		#sample metadata
		my $gisaid_done = "$projDir/gisaid_submission.done";
		my $metadata_out= "$projDir/gisaid_submission.txt";
		my $selected_consensus = $opt{'metadata-sample-consensus'};
		$selected_consensus =~ s/::/ /g;
		my ($consensus_ref,$consensus_cov,$consensus_depth) = split(/\s+/,$selected_consensus);
		$consensus_ref =~ s/_?\d?$//g;
		$selected_consensus = "Consensus to ". $selected_consensus;
		writeMetaData($projDir,$selected_consensus);
		#gisaid profile
		writeProfile($profileDir);
		#create an input file for submission pipeline
		writeSubmissionFile($projDir,$selected_consensus);

		unlink "$projDir/HTML_Report/.complete_report_web";
		if ($action eq "update"){
			&returnParamsStatus();
		}

		my $cmd;
		if ($action eq "download"){
			# download
			my $zip_file = "$projDir/GISAID/${prealname}_gisaid_data.zip";
			my $zip_file_rel = "$projDir_rel/GISAID/${prealname}_gisaid_data.zip";
			$cmd = "zip -j $zip_file $metadata_out $projDir/ReadsBasedAnalysis/readsMappingToRef/${consensus_ref}*_consensus.fasta 2>/dev/null";
			`$cmd`;
			$msg->{PATH} = $zip_file_rel;
			addMessage("DOWNLOAD","failure","failed to zip metadata for downloading") unless (-e $zip_file);
			&returnParamsStatus();
		}
		#
		#call gisaid upload 
		my $fasta = `ls $projDir/ReadsBasedAnalysis/readsMappingToRef/${consensus_ref}*_consensus.fasta`;
		chomp $fasta;
		$cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/gisaid_EpiCoV_uploader.py -m $metadata_out -p " . $opt{'metadata-gisaid-pw'} . " -u ". $opt{'metadata-gisaid-id'};
		$cmd .= " -f $fasta --headless | tee $projDir/GISAID/submit.log";
		open (my $fh,">","$projDir/GISAID/submit.sh");
		print $fh "source $EDGE_HOME/thirdParty/Anaconda3/bin/activate base\n";
		print $fh $cmd,"\n";
		close $fh;
		#my $return=`/bin/bash $projDir/GISAID/submit.sh`;
		my $script_fh;
		my $pid = open ($script_fh, "-|")
			or exec ('/bin/bash', "$projDir/GISAID/submit.sh");
		#create .done file after successful submission
		while(<$script_fh>){
			if ($_ =~ /Complete/i){
				&touchFile($gisaid_done);
			}
		}
		if ( ! -e "$gisaid_done"){
				$msg->{SUBMISSION_STATUS}="failure";
		}
		unlink "$projDir/GISAID/submit.sh";
		
	} else {
		addMessage("$action","failure","Cannot find Project $prealname direcotry!");
	}

	#delete HTML_Report/.complete_report_web
	unlink "$projDir/HTML_Report/.complete_report_web";

}elsif( $action eq "batch-download" or $action eq "batch-update" or $action eq "batch-upload2gisaid"){
	# init status
	$msg->{SUBMISSION_STATUS}="success";
	my @selectedProjCodes = split /,/,$opt{proj};
	my %selectedProjCode = map { $_ => 1 } @selectedProjCodes;
	my $short_username = $1 if $username =~ /(.*)\@/; 
	$short_username =~ s/\W/_/;
	my $metadata_out_dir = "$edgeui_output/sample_metadata_export/$short_username/". md5_hex(join ('',@selectedProjCodes));
	my $relative_outdir = "$relpath/sample_metadata_export/$short_username/". md5_hex(join ('',@selectedProjCodes));
	unlink $metadata_out_dir;
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $projects = join(",",map { "$edgeui_output/$_" } @selectedProjCodes);
	my @projCodes = split /[\x0]/, $opt{'metadata-projcodes'};
	my @projNames =  split /[\x0]/, $opt{'metadata-projnames'};
	my @virusNames = split /[\x0]/, $opt{'metadata-virus-names'};
	my @virusPassages = split /[\x0]/, $opt{'metadata-virus-passages'};
	my @sampleCollectionDates = split /[\x0]/, $opt{'metadata-sample-collection-dates'};
	my @sampleLocations = split /[\x0]/, $opt{'metadata-sample-locations'};
	my @sampleHosts = split /[\x0]/, $opt{'metadata-sample-hosts'};
	my @sampleGenders = split /[\x0]/,$opt{'metadata-sample-genders'};
	my @sampleAges = split /[\x0]/,$opt{'metadata-sample-ages'};
	my @sampleStatus = split /[\x0]/,$opt{'metadata-sample-status'};
	my @sampleSequencingTech = split /[\x0]/,$opt{'metadata-sample-sequencing-techs'};
	my @sampleConsensus = split /[\x0]/,$opt{'metadata-sample-consensus'};
	my @gisaidDoneFiles;
	foreach my $i (0..$#projCodes){
		next if ! $selectedProjCode{$projCodes[$i]};
		$opt{'metadata-virus-name'} = $virusNames[$i];
		$opt{'metadata-virus-passage'} = $virusPassages[$i];
		$opt{'metadata-sample-collection-date'} = $sampleCollectionDates[$i];
		$opt{'metadata-sample-location'} = $sampleLocations[$i];
		$opt{'metadata-sample-host'} = $sampleHosts[$i];
		$opt{'metadata-sample-gender'} = $sampleGenders[$i];
		$opt{'metadata-sample-age'} = $sampleAges[$i];
		$opt{'metadata-sample-status'} = $sampleStatus[$i];
		$opt{'metadata-sample-sequencing-tech'} = $sampleSequencingTech[$i];
		my $selected_consensus = $sampleConsensus[$i];
		$selected_consensus =~ s/::/ /g;
		my ($consensus_ref,$consensus_cov,$consensus_depth) = split(/\s+/,$selected_consensus);
		$consensus_ref =~ s/_?\d?$//g;
		$selected_consensus = "Consensus to ". $selected_consensus;
		#	print STDERR join("\t",$i, $virusNames[$i], $opt{'metadata-virus-name'}, $projNames[$i], $projCodes[$i],"\n");
		my $projDir = $edgeui_output . "/". $projCodes[$i];
		my $gisaidDone = "$projDir/gisaid_submission.done"; 
		push @gisaidDoneFiles, $gisaidDone;
		my $ownProjectFlag = 1 if ($ownProjlist->{$projCodes[$i]});
		addMessage("BATCH-SUBMIT","failure","Not the owner of project $projNames[$i]") if ( $action eq "batch-upload2gisaid" && !$ownProjectFlag);
		addMessage("BATCH-SUBMIT","failure","$projNames[$i] had submmited") if ( $action eq "batch-upload2gisaid" && $gisaidDone);
		checkParams($projNames[$i],$i);
		writeMetaData($projDir,$selected_consensus);
		#gisaid profile
		writeProfile($profileDir);
		#create an input file for submission pipeline
		writeSubmissionFile($projDir,$selected_consensus);
		unlink "$projDir/HTML_Report/.complete_report_web";
	}
	if ($action eq "batch-update"){
		&returnParamsStatus();
	}
	my $date_str = strftime "%Y%m%d", localtime;
	my $metadata_out = "$metadata_out_dir/${date_str}_edge_covid19_metadata.xls";
	my $all_sequences = "$metadata_out_dir/all_sequences.fasta";
	my $metadata_out_tsv = "$metadata_out_dir/${date_str}_edge_covid19_metadata.tsv";
	system("$EDGE_HOME/scripts/metadata/export_metadata_xls.pl", "-out", "$metadata_out", "-projects", "$projects", "-udir",$profileDir);
	if ($action eq "batch-download"){
		my $zip_file = "$metadata_out_dir/${date_str}_edge_covid19_consensus_fasta_metadata.zip";
		my $zip_file_rel = "$relative_outdir/${date_str}_edge_covid19_consensus_fasta_metadata.zip";
		my $cmd = "zip -j $zip_file $metadata_out $all_sequences $metadata_out_tsv 2>/dev/null";
		`$cmd`;
		$msg->{PATH} = $zip_file_rel;

		addMessage("BATCH-DOWNLOAD","failure","failed to export sample metadata to .xls file") unless (-e "$metadata_out" );
		addMessage("BATCH-DOWNLOAD","failure","failed to export sample consensus fasta files") unless (-e "$metadata_out_dir/all_sequences.fasta" );
		addMessage("BATCH-DOWNLOAD","failure","failed to zip exported file") unless (-e $zip_file );
		&returnParamsStatus();
	}
	if ($msg->{SUBMISSION_STATUS}="success" && $action eq "batch-upload2gisaid"){
		my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/gisaid_EpiCoV_batch_uploader.py -m $metadata_out -f $all_sequences -p " . $opt{'metadata-gisaid-pw'} . " -u ". $opt{'metadata-gisaid-id'};
		$cmd .= " --headless | tee $metadata_out_dir/GISAID_submit.log";
		open (my $fh,">","$metadata_out_dir/GISAID_submit.sh");
		print $fh "source $EDGE_HOME/thirdParty/Anaconda3/bin/activate base\n";
		print $fh $cmd,"\n";
		close $fh;
		#my $return=`/bin/bash $projDir/GISAID/submit.sh`;
		my $script_fh;
		my $pid = open ($script_fh, "-|")
			or exec ('/bin/bash', "$metadata_out_dir/GISAID_submit.sh");
		#create .done file after successful submission
		my $upload_success_flag=0;
		my $error_msg;
		while(<$script_fh>){
			if ($_ =~ /ERROR.*:(.*)/i){
				$error_msg = $1;
				addMessage("BATCH-SUBMISSION","failure","$error_msg"); 
			}
			if ($_ =~ /Complete/i){
				map {&touchFile($_)} @gisaidDoneFiles;
				$upload_success_flag=1 if (!$error_msg);
			}
		}
		if ( ! $upload_success_flag ){
			addMessage("BATCH-SUBMISSION","failure","Please see $relative_outdir/GISAID_submit.log");
		}
		unlink "$metadata_out_dir/GISAID_submit.sh";

	}
}else {
	addMessage("$action","failure","$action not recognized.");
}

if ($action !~ /create.*form/){
	&returnParamsStatus();
}

######################################################
sub touchFile{
    my $file=shift;
    open (my $fh,">",$file) or die "$!";
    close $fh;
}
sub writeSubmissionFile{
	my $outdir = shift;
	my $selected_consensus = shift;
	my $metadata_out = "$outdir/gisaid_submission.txt";
	open my $ofh,  ">$metadata_out";
	print $ofh <<"OUTMSG";
virus_name=$opt{'metadata-virus-name'}
virus_passage=$opt{'metadata-virus-passage'}
collection_date=$opt{'metadata-sample-collection-date'}
location=$opt{'metadata-sample-location'}
host=$opt{'metadata-sample-host'}
gender=$opt{'metadata-sample-gender'}
age=$opt{'metadata-sample-age'}
status=$opt{'metadata-sample-status'}
sequencing_technology=$opt{'metadata-sample-sequencing-tech'}
assembly_method=$opt{'metadata-sample-assembly-method'}
originating_lab=$opt{'metadata-orig-lab'}
originating_address=$opt{'metadata-orig-address'} 
submitting_lab=$opt{'metadata-sub-lab'}
submitting_address=$opt{'metadata-sub-address'}
authors=$opt{'metadata-authors'}
submitter=$opt{'metadata-submitter'}
gisaid_id=$opt{'metadata-gisaid-id'}
coverage=$selected_consensus
OUTMSG

	#print $ofh "gisaid_pw=".$opt{'metadata-gisaid-pw'}."\n";#
	close $ofh;
}

sub writeMetaData{
	my $outdir = shift;
	my $selected_consensus = shift;
	my $metadata_out = "$outdir/metadata_gisaid.txt";
	my @virusName = split /[\x0]/, $opt{'metadata-virus-names'};
	open my $ofh,  ">$metadata_out";
	print $ofh <<"OUTMSG";
virus_name=$opt{'metadata-virus-name'}
virus_passage=$opt{'metadata-virus-passage'}
collection_date=$opt{'metadata-sample-collection-date'}
location=$opt{'metadata-sample-location'}
host=$opt{'metadata-sample-host'}
gender=$opt{'metadata-sample-gender'}
age=$opt{'metadata-sample-age'}
status=$opt{'metadata-sample-status'}
sequencing_technology=$opt{'metadata-sample-sequencing-tech'}
coverage=$selected_consensus
OUTMSG
	close $ofh;
}

sub writeProfile{
	my $dir = shift;
	my $profile_out =  "$dir/gisaid_submission_profile.txt";
	if(-e $profile_out) {
		open my $fh, $profile_out or die "Can't open $profile_out $!";
		while(<$fh>){
			chomp;
			next if(/^#/);
			if ( /(.*)=(.*)/ ){
				$opt{'metadata-submitter'} =$2 if ($1 eq "submitter" && !$opt{'metadata-submitter'});
				$opt{'metadata-gisaid-id'} =$2 if ($1 eq "gisaid_id" && !$opt{'metadata-gisaid-id'});
			}
		}
		close $fh;
	}

	open my $ofh,  ">$profile_out";
	print $ofh <<"OUTMSG";
originating_lab=$opt{'metadata-orig-lab'}
originating_address=$opt{'metadata-orig-address'}
submitting_lab=$opt{'metadata-sub-lab'}
submitting_address=$opt{'metadata-sub-address'}
authors=$opt{'metadata-authors'}
submitter=$opt{'metadata-submitter'}
gisaid_id=$opt{'metadata-gisaid-id'}
OUTMSG
	close $ofh;
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
            addMessage("getUserProjFromDB","failure",$1);
            return;
    }
    my $array_ref =  from_json($result_json);
	#print Dumper($array_ref) if @ARGV;
	foreach my $hash_ref (@$array_ref)
	{
		my $id = $hash_ref->{id};
		my $project_name = $hash_ref->{name};
		my $projCode = $hash_ref->{code};
		my $status = $hash_ref->{status};
		next if ($status =~ /delete/i);
		next if (! -r "$edgeui_output/$id/process.log" && ! -r "$edgeui_output/$projCode/process.log");
		$list->{$id}->{PROJNAME} = $id;
		$list->{$id}->{REAL_PROJNAME} = $project_name;
		$list->{$id}->{PROJCODE} = $projCode;
		$list->{$id}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}"; 
		$list->{$id}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$id}->{PROJ_TYPE} = $hash_ref->{type};
		$list->{$projCode}->{PROJNAME} = $id;
		$list->{$projCode}->{REAL_PROJNAME} = $project_name;
		$list->{$projCode}->{PROJID} = $id;
		$list->{$projCode}->{OWNER} = "$hash_ref->{owner_firstname} $hash_ref->{owner_lastname}"; 
		$list->{$projCode}->{OWNER_EMAIL} = $hash_ref->{owner_email};
		$list->{$projCode}->{PROJ_TYPE} = $hash_ref->{type};
	}
	return $list;
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
        my $service= ($userType =~ /admin/i)? "WS/user/admin/getProjects" :"WS/project/getInfo";
        #w Set the request parameters
        my $url = $um_url . $service;
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
        $result = shift @$result if (ref($result) eq "ARRAY");
        if ($result->{error_msg})
        {
			addMessage("getProjNameFromDB","failure",$result->{error_msg});
        }
        else{
		return ($result->{name} , $result->{code}, $result->{status} , $result->{owner_email});
        }
}
sub getUserInfo {
	my $username=shift;
	my $password=shift;
	my $msg=shift;
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
		addMessage("getUserInfo","failure", $1);
		#$msg->{ERROR}=$1;
		#return;
        }
        my $tmp_r = from_json($result_json);
        $user_info{$_} = $tmp_r->{$_} foreach (keys %$tmp_r);
        return \%user_info;
}

sub scanProjToList{
        my $out_dir = shift;
        my $pname = shift;
        my $config_file;
        my $processLog;
        my ($projid,$projCode,$projName,$projStatus);
        if ($pname && -d "$out_dir/$pname"){
                $config_file = "$out_dir/$pname/config.txt";
                $processLog = "$out_dir/$pname/process_current.log";
        }else{
                $config_file = `grep -a "projid=$pname" $out_dir/*/config.txt | awk -F':' '{print \$1}'`;
        }
        chomp $config_file;
        return ($projName,$projCode,$projStatus) if ( ! -e $config_file);
        if ( -r "$processLog"){
                open (my $fh, $processLog);
                while(<$fh>){
                        if (/queued/){
                                $projStatus="unstarted";
                                last;
                        }
                        if (/^All Done/){
                                $projStatus="finished";
                        }
                }
                close $fh;
        }
        open (my $fh, $config_file) or die "Cannot read $config_file\n";
        while(<$fh>){
                last if (/^\[Down/);
                $projid=$1 if (/^projid=(\S+)/);
                $projCode=$1 if (/^projcode=(\S+)/);
                $projName=$1 if (/^projname=(\S+)/);

        }
        close $fh;
        return ($projName,$projid,$projStatus);
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


sub stringSanitization{
	my $opt=shift;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		if ($key =~ /metadata|edgesite|locality|administrative|country|lat|lng/){
			$opt->{$key} =~ s/[`";'&|]/ /g;
			next;
		}
		if($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/){
			addMessage("CHECK_INPUT","failure","Invalid characters detected \'$str\'.");
			&returnParamsStatus();
		}
	}
}

sub returnParamsStatus {
	# return json
	my $json = encode_json($msg);
	print $cgi->header('application/json'), $json;
	exit 0;
}

sub addMessage {
	my ($job, $status, $note, $index) = @_;
	$index = 0 if (!$index);
	if( $job eq "PARAMS" ){
		#e.g.: $msg->{"PARAMS"}->{"edge-ref-file"}="File not found";
		$msg->{$job}->{$index}->{$status}=$note;
		$msg->{SUBMISSION_STATUS}="failure";
	}
	else{
		#&addMessage("PROJECT_NAME","info","$info");
		$msg->{$job}->{STATUS}=$status;
		$msg->{$job}->{NOTE}=$note if defined $note;
		$msg->{SUBMISSION_STATUS}=$status if $status eq "failure";
	}
}

sub checkParams {
	my $projname = shift;
	$projname .= ": " if ($projname);
	my $index = shift;
	&addMessage("PARAMS", "metadata-virus-name", "$projname Virus name is required." , $index) unless ( $opt{'metadata-virus-name'}); 
	&addMessage("PARAMS", "metadata-virus-passage", "$projname Passage details/history is required." , $index) unless ( $opt{'metadata-virus-passage'}); 
	&addMessage("PARAMS", "metadata-sample-collection-date", "$projname Collection date is required." , $index) unless ( $opt{'metadata-sample-collection-date'}); 
	&addMessage("PARAMS", "metadata-sample-location", "$projname Location is required.", $index) unless ( $opt{'metadata-sample-location'}); 
	&addMessage("PARAMS", "metadata-sample-host", "$projname Host is required.", $index) unless ( $opt{'metadata-sample-host'});
	&addMessage("PARAMS", "metadata-sample-gender", "$projname Gender is required.", $index) unless ( $opt{'metadata-sample-gender'});
	&addMessage("PARAMS", "metadata-sample-age", "$projname Age is required.", $index) unless ( $opt{'metadata-sample-age'});
	&addMessage("PARAMS", "metadata-sample-status", "$projname Status is required.", $index) unless ( $opt{'metadata-sample-status'});
	&addMessage("PARAMS", "metadata-sample-sequencing-tech", "$projname Sequencing Technology is required.", $index) unless ( $opt{'metadata-sample-sequencing-tech'});
	&addMessage("PARAMS", "metadata-sample-consensus", "$projname Consensus FASTA is required.", $index) unless ( $opt{'metadata-sample-consensus'} );


	&addMessage("PARAMS", "metadata-orig-lab", "Originating lab is required.") unless ( $opt{'metadata-orig-lab'});
	&addMessage("PARAMS", "metadata-orig-address", "Originating address is required.") unless ( $opt{'metadata-orig-address'});
	&addMessage("PARAMS", "metadata-org-lab", "Submitting lab is required.") unless ( $opt{'metadata-sub-lab'});
	&addMessage("PARAMS", "metadata-sub-address", "Submitting address is required.") unless ( $opt{'metadata-sub-address'});
	&addMessage("PARAMS", "metadata-authors", "Authors is required.") unless ( $opt{'metadata-authors'});
	if ($action eq 'upload2gisaid' or $action eq 'batch-upload2gisaid'){
		&addMessage("PARAMS", "metadata-submitter", "Submitter is required.") unless ( $opt{'metadata-submitter'});
		&addMessage("PARAMS", "metadata-gisaid-id", "GISAID id is required.") unless ( $opt{'metadata-gisaid-id'});
		&addMessage("PARAMS", "metadata-gisaid-pw", "GISAID password is required.") unless ( $opt{'metadata-gisaid-pw'});
	}
}
