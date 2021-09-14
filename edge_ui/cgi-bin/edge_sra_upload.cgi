#!/usr/bin/env perl

use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use lib "../../lib";
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use File::Basename;
use File::Path qw(make_path remove_tree);
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
my $debug= 0;

######################################################################################
# DATA STRUCTURE:
#	$msg->{"PARAMS"}->{0}->{metadata-virus-name}  =
#			     ->{metadata-sample-host} =
#			     ...
# 	                ->{1}...

#	$msg->{SUBMISSION_STATUS}
#	$msg->{PATH}   # file for action download or log to display
#	$msg->{PID}   # process id of child program exectue in the background
     
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
	my $outHtml = "$projDir/UPLOAD/upload.html";
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/metadata/outputSRA_Metadata_w_temp.pl $projDir $outHtml $prealname $profileDir";
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

}elsif( $action eq 'create-batch-form'){
	my @projCodes = split /,/,$opt{proj};
	my $short_username = $1 if $username =~ /(.*)\@/;
        $short_username =~ s/\W/_/;
	my $metadata_out_dir = "$edgeui_output/sample_metadata_export/$short_username/". md5_hex(join ('',@projCodes));
	chdir $edgeui_wwwroot;
	my $outHtml = "$metadata_out_dir/batch_metadata.html";
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $projects = join(",",map { "$edgeui_output/$_" } @projCodes);

	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/metadata/outputSRA_Metadata_w_temp_table.pl -projects $projects -out $outHtml -udir $profileDir";
	`$cmd`;
	#print STDERR "$cmd";

	open REP, "$outHtml" or die "Can't open $outHtml: $!";
	my $html_content;
	while(<REP>){
		$html_content .= $_;
	}
	print "Content-Type: text/html\n\n",
		$html_content;
}elsif($action eq "upload2sra" or $action eq "download" or $action eq "update") {
	$msg->{SUBMISSION_STATUS}="success";

	# validating parameters
	&checkParams() if ($action eq 'upload2sra');
	if($msg->{SUBMISSION_STATUS} eq "failure") {
		returnParamsStatus();
	}

	my $projDir = $edgeui_output . "/". $projCode;
	my $upload_content_dir = "$projDir/UPLOAD";
	my $projDir_rel = $relpath . "/". $projCode;
	my $upload_content_reldir = "$projDir_rel/UPLOAD";
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $projCompleteReport_cache = "$projDir/HTML_Report/.complete_report_web";
	make_dir($upload_content_dir);
	if(-e $projDir) {
		#sample metadata
		my $sra_done = "$upload_content_dir/ncbi_sra_submission.done";
		my $metadata_project  = "$upload_content_dir/sra_project.txt";
		my $metadata_sample   = "$upload_content_dir/sra_samples.txt";
		my $metadata_exp      = "$upload_content_dir/sra_experiments.txt";
		my $metadata_filelist    = "$upload_content_dir/sra_filelist.txt";
		my $metadata_other_info  = "$upload_content_dir/sra_additional_info.txt";
		
		my $date_str = strftime "%Y%m%d%H%M%S", localtime;
		my $sra_submit_data_dir = "$upload_content_dir/sra_transfer_".$date_str.$pname;
		my ($input_pe_fastq, $input_se_fastq) = get_input_fastq($projDir,$sra_submit_data_dir);
		
		writeBioProject($metadata_project);
		writeBioSamples($projDir,$metadata_sample);
		writeExperiment($projDir,$metadata_exp,$metadata_filelist, $input_pe_fastq, $input_se_fastq);
		wrtieAdditional($projDir,$metadata_other_info);

		unlink "$projCompleteReport_cache";
		if ($action eq "update"){
			remove_tree($sra_submit_data_dir);
			&returnParamsStatus();
		}

		if ($action eq "download"){
			# download
			my $zip_file = "$upload_content_dir/${prealname}_sra_metadata.zip";
			my $zip_file_rel = "$upload_content_reldir/${prealname}_sra_metadata.zip";
			my $cmd = "zip -j $zip_file $metadata_project $metadata_sample $metadata_exp $metadata_filelist $metadata_other_info 2>/dev/null";
			`$cmd`;
			$msg->{PATH} = $zip_file_rel;
			addMessage("DOWNLOAD","failure","failed to zip metadata for downloading") unless (-e $zip_file);
			remove_tree($sra_submit_data_dir);
			&returnParamsStatus();
		}
		#
		#call gisaid upload 
		my $submit_log = "$upload_content_dir/sra_submit.log";
		my $submit_current_log = "$upload_content_dir/sra_submit_current.log";
		my $submit_script = "$upload_content_dir/sra_submit.sh";
		## prepare submission.xml file
		my $ncbi_cmd .= "$EDGE_HOME/scripts/sra_submit/sra_xml.py --input-dir $sra_submit_data_dir --listfile $metadata_filelist --projectfile $metadata_project --metadatafile $metadata_sample ";
		$ncbi_cmd .= "--libselection \'$opt{'metadata-sra-meta-libselection'}\' --libstrategy \'$opt{'metadata-sra-meta-libstrategy'}\' --datatype 'Raw sequence reads' --instrument \'$opt{'metadata-sra-meta-libmodel'}\' ";
		$ncbi_cmd .= "--libsource \'$opt{'metadata-sra-meta-libsource'}\'  --platform \'$opt{'metadata-sra-meta-platform'}\' --samplescope eMultiisolate --liblayout \'$opt{'metadata-sra-meta-liblayout'}\' ";
		$ncbi_cmd .= "--libdesign \'$opt{'metadata-sra-meta-design-description'}\'  --libtitle \'$opt{'metadata-sra-meta-title'}\'  | tee -a $submit_log ";
		## file transfer commands;
		my $ncbi_sra_dir = ($debug)? "submit/Test/": "submit/Production/" ;
		$ncbi_cmd .= "\n cd $upload_content_dir \n";
		$ncbi_cmd .= "\n echo 'Transfer files to NCBI' \n";
		#$ncbi_cmd .= "\n sleep 30 \n";
		$ncbi_cmd .= "$EDGE_HOME/scripts/sra_submit/sra_ascp.py --ncbi-sra-dir \'$ncbi_sra_dir\' --input-dir \'". basename($sra_submit_data_dir) . "\/\' ";
		$ncbi_cmd .= " --ncbi-username ". $sys->{'sra_submission_account'} ;
		$ncbi_cmd .= " --ncbi-private-key ". $sys->{'sra_acsp_keyfile'} .  " | tee -a $submit_log \n";
		
		SetSubmitScript($submit_script, $submit_log, $projCompleteReport_cache, $sra_done, $ncbi_cmd);
		my $script_fh;

		my $pid;
		my $pid = open ($script_fh, "-|")
			or exec ("/bin/bash $submit_script > $submit_current_log 2>\&1 &");
		if (not defined $pid ){
			$msg->{SUBMISSION_STATUS}="failure";
		}else{
			$msg->{PID} = ++$pid;
			$msg->{PATH} = "$upload_content_reldir/sra_submit_current.log";
			addMessage("SUBMIT",'success',"See the log window for status.");
		}
		#or exec ('/bin/bash', "$upload_content_dir/submit.sh");
		#create .done file after successful submission
		#while(<$script_fh>){
		#	if ($_ =~ /Complete/i){
		#		&touchFile($gisaid_done);
		#	}
		#}
		#if ( ! -e "$gisaid_done"){
		#		$msg->{SUBMISSION_STATUS}="failure";
		#}
		#unlink $submit_script if (! $debug );
		#unlink $projCompleteReport_cache;
		
	} else {
		addMessage("$action","failure","Cannot find Project $prealname direcotry!");
	}

	#delete HTML_Report/.complete_report_web
	unlink $projCompleteReport_cache;

}elsif( $action eq "batch-download" or $action eq "batch-update" or $action eq "batch-upload2sra" or $action eq "batch-biosample-template-download" or $action eq "batch-bioexperiment-template-download"){
	# init status
	$msg->{SUBMISSION_STATUS}="success";
	my @selectedProjCodes = split /,/,$opt{proj};
	my %selectedProjCode = map { $_ => 1 } @selectedProjCodes;
	my $short_username = $1 if $username =~ /(.*)\@/; 
	$short_username =~ s/\W/_/;
	my $metadata_out_dir = "$edgeui_output/sample_metadata_export/$short_username/". md5_hex(join ('',@selectedProjCodes));
	my $relative_outdir = "$relpath/sample_metadata_export/$short_username/". md5_hex(join ('',@selectedProjCodes));
	unlink $metadata_out_dir;
	make_dir($metadata_out_dir);
	my $date_str = strftime "%Y%m%d%H%M%S", localtime;
	my $sra_submit_data_dir = "$metadata_out_dir/sra_transfer_".$date_str.md5_hex(join ('',@selectedProjCodes));
	my $batch_sra_bioproject = "$metadata_out_dir/sra_project.txt";
	my $batch_sra_biosample = "$metadata_out_dir/sra_samples.txt";
	my $batch_sra_experiment = "$metadata_out_dir/sra_experiments.txt";
	my $batch_sra_filelist = "$metadata_out_dir/sra_filelist.txt";
	my $batch_sra_other_info = "$metadata_out_dir/sra_additional_info.txt";
	unlink $batch_sra_experiment, $batch_sra_biosample, $batch_sra_filelist;
	my $profileDir = $sys->{edgeui_input}."/$userDir";
	my $projects = join(",",map { "$edgeui_output/$_" } @selectedProjCodes);
	## Data from UI input
	my @projCodes = split /[\x0]/, $opt{'metadata-projcodes'};
	my @projNames =  split /[\x0]/, $opt{'metadata-projnames'};
	my @sampleNames = split /[\x0]/, $opt{'metadata-sra-biosample-names'};
	my @isolates = split /[\x0]/, $opt{'metadata-sra-biosample-isolates'};
	my @isolateSources = split /[\x0]/, $opt{'metadata-sra-biosample-isolate-sources'};
	my @sampleLocations = split /[\x0]/, $opt{'metadata-sra-biosample-locations'};
	my @samplePassages = split /[\x0]/, $opt{'metadata-sra-biosample-passages'};
	my @sampleCollectionDates = split /[\x0]/, $opt{'metadata-sra-biosample-collection-dates'};
	my @sampleCollectionBys = split /[\x0]/, $opt{'metadata-sra-biosample-collect-bys'};
	my @sampleLatlons = split /[\x0]/, $opt{'metadata-sra-biosample-latlons'};
	my @sampleHosts = split /[\x0]/, $opt{'metadata-sra-biosample-hosts'};
	my @sampleGenders = split /[\x0]/,$opt{'metadata-sra-biosample-genders'};
	my @sampleAges = split /[\x0]/,$opt{'metadata-sra-biosample-ages'};
	my @sampleStatus = split /[\x0]/,$opt{'metadata-sra-biosample-statuss'};
	my @samplePurposes = split /[\x0]/,$opt{'metadata-sra-biosample-purposes'};
	my @sampleGISAIDaccs = split /[\x0]/,$opt{'metadata-sra-biosample-gisaid-accs'};
	my @expTitles = split /[\x0]/,$opt{'metadata-sra-meta-titles'};
	my @expDesigns = split /[\x0]/,$opt{'metadata-sra-meta-design-descriptions'};
	my @expLibSelections = split /[\x0]/,$opt{'metadata-sra-meta-library-selections'};
	my @expLibStrategys = split /[\x0]/,$opt{'metadata-sra-meta-library-strategys'};
	my @expLibLayouts = split /[\x0]/,$opt{'metadata-sra-meta-library-layouts'};
	my @expLibSources = split /[\x0]/,$opt{'metadata-sra-meta-library-sources'};
	my @expPlatforms = split /[\x0]/,$opt{'metadata-sra-meta-platforms'};
	my @expLibModels = split /[\x0]/,$opt{'metadata-sra-meta-library-models'};
	my @expSeqPurposes = split /[\x0]/,$opt{'metadata-sra-meta-sequencing-purposes'};
	my $bioprojectID = $opt{'metadata-sra-bioproject-id'};
	my $bioprojectTitle = $opt{'metadata-sra-bioproject-title'};
	my $bioprojectDesc = $opt{'metadata-sra-bioproject-desc'};
	my $bioprojectLinkDesc = $opt{'metadata-sra-bioproject-link-desc'};
	my $bioprojectLinkUrl = $opt{'metadata-sra-bioproject-link-url'};

	writeBioProject($batch_sra_bioproject);
	wrtieAdditional($batch_sra_other_info);
	my @SRADoneFiles;
	my @projCompleteReport_cacheFiles;
	foreach my $i (0..$#projCodes){
		next if ! $selectedProjCode{$projCodes[$i]};
		# save for each project checkParams
		$opt{'metadata-sra-biosample-name'} = $sampleNames[$i];
		$opt{'metadata-sra-biosample-isolate'} = $isolates[$i];
		$opt{'metadata-sra-biosample-isolate-source'} = $isolateSources[$i];
		$opt{'metadata-sra-biosample-location'} = $sampleLocations[$i];
		$opt{'metadata-sra-biosample-passage'} = $samplePassages[$i];
		$opt{'metadata-sra-biosample-collection-date'} = $sampleCollectionDates[$i];
		$opt{'metadata-sra-biosample-collect-by'} = $sampleCollectionBys[$i];
		$opt{'metadata-sra-biosample-latlon'} = $sampleLatlons[$i];
		$opt{'metadata-sra-biosample-host'} = $sampleHosts[$i];
		$opt{'metadata-sra-biosample-gender'} = $sampleGenders[$i];
		$opt{'metadata-sra-biosample-status'} = $sampleStatus[$i];
		$opt{'metadata-sra-biosample-age'} = $sampleAges[$i];
		$opt{'metadata-sra-biosample-purpose'} = $samplePurposes[$i];
		$opt{'metadata-sra-biosample-gisaid'} = $sampleGISAIDaccs[$i];
		$opt{'metadata-sra-meta-title'} = $expTitles[$i];
		$opt{'metadata-sra-meta-design-description'} = $expDesigns[$i];
		$opt{'metadata-sra-meta-libselection'} = $expLibSelections[$i];
		$opt{'metadata-sra-meta-libstrategy'} = $expLibStrategys[$i];
		$opt{'metadata-sra-meta-liblayout'} = $expLibLayouts[$i];
		$opt{'metadata-sra-meta-libsource'} = $expLibSources[$i];
		$opt{'metadata-sra-meta-platform'} = $expPlatforms[$i];
		$opt{'metadata-sra-meta-libmodel'} = $expLibModels[$i];
		$opt{'metadata-sra-meta-sequencing-purpose'} = $expSeqPurposes[$i];
		
		my $projDir = $edgeui_output . "/". $projCodes[$i];
		my $upload_content_dir = "$projDir/UPLOAD";
		make_dir($upload_content_dir);
		my $sra_done = "$upload_content_dir/ncbi_sra_submission.done";
		my $metadata_project  = "$upload_content_dir/sra_project.txt";
		my $metadata_sample   = "$upload_content_dir/sra_samples.txt";
		my $metadata_exp      = "$upload_content_dir/sra_experiments.txt";
		my $metadata_filelist    = "$upload_content_dir/sra_filelist.txt";
		my $metadata_other_info  = "$upload_content_dir/sra_additional_info.txt";
		my ($input_pe_fastq, $input_se_fastq) = get_input_fastq($projDir,$sra_submit_data_dir);
		my $projCompleteReport_cache = "$projDir/HTML_Report/.complete_report_web";
		push @SRADoneFiles, $sra_done;
		push @projCompleteReport_cacheFiles, $projCompleteReport_cache;
		my $ownProjectFlag = 1 if ($ownProjlist->{$projCodes[$i]});
		addMessage("BATCH-SUBMIT","failure","Not the owner of project $projNames[$i]") if ( $action eq "batch-upload2sra" && !$ownProjectFlag);
		addMessage("BATCH-SUBMIT","failure","$projNames[$i] had submmited") if ( $action eq "batch-upload2sra" && -e $sra_done);
		checkParams($projNames[$i],$i) if $action eq "batch-upload2sra";
		
		writeBioProject($metadata_project);
		writeBioSamples($projDir,$metadata_sample,0,0);  
		writeBioSamples($metadata_out_dir,$batch_sra_biosample,$i,1);  ## append data table
		writeExperiment($projDir,$metadata_exp,$metadata_filelist, $input_pe_fastq, $input_se_fastq,0,0);
		writeExperiment($metadata_out_dir,$batch_sra_experiment,$batch_sra_filelist, $input_pe_fastq, $input_se_fastq,$i,1); ## append data table
		wrtieAdditional($metadata_other_info);
		
		unlink $projCompleteReport_cache;
	}
	if ($action eq "batch-update"){
		remove_tree($sra_submit_data_dir);
		&returnParamsStatus();
	}
	my $date_str2 = strftime "%Y%m%d", localtime;
	my $download_biosample_template_txt = "$metadata_out_dir/metadata_biosample_template.tsv";
	my $download_biosample_template_txt_rel = "$relative_outdir/metadata_biosample_template.tsv";
	my $download_bioexperiment_template_txt = "$metadata_out_dir/metadata_sra_experiment_template.tsv";
	my $download_bioexperiment_template_txt_rel = "$relative_outdir/metadata_sra_experiment_template.tsv";
	if ($action eq "batch-biosample-template-download" or $action eq "batch-bioexperiment-template-download"){
		my $zip_file = "$metadata_out_dir/${date_str2}_edge_covid19_sra_biosample_template.tsv.zip";
		my $zip_file_rel = "$relative_outdir/${date_str2}_edge_covid19_sra_biosample_template.tsv.zip";
		my $zip_file2 = "$metadata_out_dir/${date_str2}_edge_covid19_sra_experiment_template.tsv.zip";
		my $zip_file_rel2 = "$relative_outdir/${date_str2}_edge_covid19_sra_experiment_template.tsv.zip";
		
		my @download_biosample_template_tsv_header = ("project-name",'sample-name','sample-isolate','sample-isolate-source','sample-location','sample-passage','sample-collection-date','sample-collect-by','sample-latlon','sample-host','sample-gender','sample-age','sample-purpose','sample-gisaid-acc');
		my @download_biosample_template_tsv_content;
		my @download_bioexperiment_template_tsv_header = ("project-name",'title','design-description','library-selection','library-strategy','library-layout','library-source','platform','library-model','sequencing-purpose');
		my @download_bioexperiment_template_tsv_content;
		foreach my $i (0..$#projCodes){
			my $download_biosample_template_tsv_string = join("\t", $projNames[$i], $sampleNames[$i],$isolates[$i],$isolateSources[$i],$sampleLocations[$i],$samplePassages[$i],$sampleCollectionDates[$i],$sampleCollectionBys[$i],$sampleLatlons[$i],$sampleHosts[$i],$sampleGenders[$i],$sampleAges[$i],$samplePurposes[$i],$sampleGISAIDaccs[$i]);
			push @download_biosample_template_tsv_content, $download_biosample_template_tsv_string;
			my $download_bioexperiment_template_tsv_string = join("\t", $projNames[$i], $expTitles[$i], $expDesigns[$i],$expLibSelections[$i],$expLibStrategys[$i],$expLibLayouts[$i],$expLibSources[$i],$expPlatforms[$i],$expLibModels[$i],$expSeqPurposes[$i]);
			push @download_bioexperiment_template_tsv_content, $download_bioexperiment_template_tsv_string;
		}
		
		if ($action eq "batch-biosample-template-download" ){
			&write_tsv($download_biosample_template_txt,join("\t",@download_biosample_template_tsv_header), join("\n",@download_biosample_template_tsv_content));
			my $cmd = "zip -j $zip_file $download_biosample_template_txt 2>/dev/null";
			`$cmd`;
			$msg->{PATH} = $zip_file_rel;
			addMessage("BATCH-TEMPLATE-DOWNLOAD","failure","failed to donwload metadata biosample template tsv file") unless  (-r "$download_biosample_template_txt" );
		}
		if ($action eq "batch-bioexperiment-template-download" ){
			&write_tsv($download_bioexperiment_template_txt,join("\t",@download_bioexperiment_template_tsv_header), join("\n",@download_bioexperiment_template_tsv_content));
			my $cmd = "zip -j $zip_file2 $download_bioexperiment_template_txt 2>/dev/null";
			`$cmd`;
			$msg->{PATH} = $zip_file_rel2;
			addMessage("BATCH-TEMPLATE-DOWNLOAD","failure","failed to donwload metadata SRA experiment template tsv file") unless  (-r "$download_bioexperiment_template_txt" );
		}
		remove_tree($sra_submit_data_dir);
		&returnParamsStatus();
	}
	
	if ($action eq "batch-download"){
		my $zip_file = "$metadata_out_dir/${date_str2}_edge_covid19_sra_metadata.zip";
		my $zip_file_rel = "$relative_outdir/${date_str2}_edge_covid19_sra_metadata.zip";
	
		my $cmd = "zip -j $zip_file $batch_sra_bioproject $batch_sra_biosample $batch_sra_experiment $batch_sra_filelist $batch_sra_other_info 2>/dev/null";
		`$cmd`;
		$msg->{PATH} = $zip_file_rel;

		addMessage("BATCH-DOWNLOAD","failure","failed to export SRA biosample metadata") unless (-e "$batch_sra_biosample" );
		addMessage("BATCH-DOWNLOAD","failure","failed to export SRA experiments metadata") unless (-e "$batch_sra_experiment" );
		addMessage("BATCH-DOWNLOAD","failure","failed to zip exported file") unless (-e $zip_file );
		remove_tree($sra_submit_data_dir);
		&returnParamsStatus();
	}
	if ($msg->{SUBMISSION_STATUS} eq "success" && $action eq "batch-upload2sra"){
		my $submit_log = "$metadata_out_dir/sra_submit.log";
		my $submit_current_log = "$metadata_out_dir/sra_submit_current.log";
		my $submit_script = "$metadata_out_dir/sra_submit.sh";
		
		## prepare submission.xml file
		my $ncbi_cmd .= "$EDGE_HOME/scripts/sra_submit/sra_xml.py --input-dir $sra_submit_data_dir --listfile $batch_sra_filelist --projectfile $batch_sra_bioproject --metadatafile $batch_sra_biosample ";
		$ncbi_cmd .= "--experimentfile $batch_sra_experiment --datatype 'Raw sequence reads' --samplescope eMultiisolate   | tee -a $submit_log ";
		## file transfer commands;
		my $ncbi_sra_dir = ($debug)? "submit/Test/": "submit/Production/" ;
		$ncbi_cmd .= "\n cd $metadata_out_dir \n";
		$ncbi_cmd .= "\n echo 'Transfer files to NCBI' \n";
		$ncbi_cmd .= "$EDGE_HOME/scripts/sra_submit/sra_ascp.py --ncbi-sra-dir \'$ncbi_sra_dir\' --input-dir \'". basename($sra_submit_data_dir) . "\/\' ";
		$ncbi_cmd .= " --ncbi-username ". $sys->{'sra_submission_account'} ;
		$ncbi_cmd .= " --ncbi-private-key ". $sys->{'sra_acsp_keyfile'} .  " | tee -a $submit_log \n";
		
		SetSubmitScript($submit_script, $submit_log, \@projCompleteReport_cacheFiles, \@SRADoneFiles, $ncbi_cmd);
		my $script_fh;
		my $pid = open ($script_fh, "-|")
			or exec ("/bin/bash $submit_script > $submit_current_log 2>\&1 &");
		if (not defined $pid ){
			$msg->{SUBMISSION_STATUS}="failure";
		}else{
			$msg->{PID} = ++$pid;
			$msg->{PATH} = "$relative_outdir/sra_submit_current.log";
			addMessage("BATCH-SUBMIT",'success',"See the log window for status.");
		}
		#if ( ! $upload_success_flag ){
		#		addMessage("BATCH-SUBMISSION","failure","Please see $relative_outdir/gisaid_submit.log ncbi_submit.log");
		#}

	}
}else {
	addMessage("$action","failure","$action not recognized.");
}

if ($action !~ /create.*form/){
	&returnParamsStatus();
}

######################################################

sub get_input_fastq {
	my $projdir = shift;
	my $transfer_dir = shift;
	make_dir($transfer_dir);
	my @PEFILES;
	my @SEFILES;
	open(my $sumfh, "$projdir/process.log") or die $!;
	while(<$sumfh>) {
		chomp;
		#parse input files
		if( /runPipeline/ ) {
			undef @PEFILES, @SEFILES; 
        }
		if( /runPipeline .*-p (.*) -\w/ || /runPipeline .*-p (.*) >/ || /runPipeline .*-p (.*)$/) {
			push @PEFILES, split /\s+/,$1;
		}
		if(/runPipeline .*-u (.*) -\w/ || /runPipeline .*-u (.*) >/ || /runPipeline .*-u (.*)$/){
			push @SEFILES, split /\s+/,$1;
		}
		if(/SRA_id=(.*)/){
            # this should not be the case since it already in SRA db.  Just for test.
        	my @SRAFILES = glob("$projdir/SRA_Download/*fastq.gz");
        	my %pair;
        	if (scalar(@SRAFILES)==1){
        		@SEFILES = @SRAFILES;
        	}else{
				foreach (@SRAFILES){
					if (/(\S+)_?.?[12]/){
						push @PEFILES, "$_";
					}else{
						push @SEFILES, "$_";
					}
				}
			}
        }
    }
	close $sumfh;
	
	foreach my $file_i (0..$#PEFILES){
		my $basename = basename($PEFILES[$file_i]);
		symlink($PEFILES[$file_i], "$transfer_dir/$basename");
		$PEFILES[$file_i] = "$transfer_dir/$basename";
	}
	foreach my $file_i (0..$#SEFILES){
		my $basename = basename($SEFILES[$file_i]);
		symlink($SEFILES[$file_i], "$transfer_dir/$basename");
		$SEFILES[$file_i] = "$transfer_dir/$basename";
	}
	
	return (\@PEFILES,\@SEFILES);
}

sub write_tsv{
        my $outfile = shift;
        my $header = shift;
        my $content = shift;
        open (my $ofh, ">", $outfile) or die "Cannot write to $outfile\n";
        print $ofh join("\n",$header,$content);
        close $ofh;
}
sub touchFile{
    my $file=shift;
    open (my $fh,">",$file) or die "$!";
    close $fh;
}

sub wrtieAdditional{
	my $outfile = shift;

	open my $ofh,  ">$outfile";
	print $ofh <<"OUTMSG";
sra-submitter=$opt{'metadata-sra-submitter'}
sra-release-date=$opt{'metadata-sra-release-date'}
OUTMSG
	close $ofh;
}

sub writeBioProject{
	my $outfile = shift;
	

	open my $ofh,  ">$outfile";
	print $ofh <<"OUTMSG";
Last	BPH
First	LANL
Center	LANL Biosecurity and Public Health
Type	Center
Website	https://edge-covid19.edgebioinformatics.org/
Organism	Severe acute respiratory syndrome coronavirus 2
Email	$opt{'metadata-sra-submitter'} 
Release	$opt{'metadata-sra-release-date'} 
Resource	{"$opt{'metadata-sra-bioproject-link-desc'}":"$opt{'metadata-sra-bioproject-link-url'}"}
OUTMSG

	if ($opt{'metadata-sra-bioproject-title'}){
		print $ofh "ProjectName\t$opt{'metadata-sra-bioproject-title'}\n"; 
	}
	if ($opt{'metadata-sra-bioproject-title'}){
		print $ofh "ProjectTitle\t$opt{'metadata-sra-bioproject-title'}\n";
	}
	if ($opt{'metadata-sra-bioproject-desc'}){
		print $ofh "Description\t$opt{'metadata-sra-bioproject-desc'}\n";
	}
	close $ofh;
}

sub writeBioSamples{
	my $outdir = shift;
	my $outfile = shift;
	my $index = shift;
	my $append = shift;
	my @sample_string;
	
	if (!$append){ unlink $outfile; }
	open (my $ofh,  ">>", "$outfile") or die "Cannot write to $outfile";
	my @header = ("sample_name", "sample_title", "organism", "isolate", "collected_by", "collection_date", 
	              "geo_loc_name", "isolation_source", "lat_lon", "host", "host_disease", 
	              "host_health_state","host_age","host_sex",
	              "passage_history", "description", "purpose_of_sampling",
	              "purpose_of_sequencing","GISAID_accession");
	if ($opt{'metadata-sra-bioproject-sw'}) { push @header, "bioproject_accession"; }
	
	$sample_string[0] = $opt{'metadata-sra-biosample-name'} ;
	$sample_string[1] = "";
	$sample_string[2] = "Severe acute respiratory syndrome coronavirus 2";
	$sample_string[3] = $opt{'metadata-sra-biosample-isolate'};
	$sample_string[4] = ($opt{'metadata-sra-biosample-collect-by'})? $opt{'metadata-sra-biosample-collect-by'} : "not collected";
	$sample_string[5] = ($opt{'metadata-sra-biosample-collection-date'})? $opt{'metadata-sra-biosample-collection-date'} : "not collected";
	$sample_string[6] = ($opt{'metadata-sra-biosample-location'})? $opt{'metadata-sra-biosample-location'} : "not collected";
	$sample_string[7] = ($opt{'metadata-sra-biosample-isolate-source'})? $opt{'metadata-sra-biosample-isolate-source'} : "not collected";
	$sample_string[8] = ($opt{'metadata-sra-biosample-latlon'})? $opt{'metadata-sra-biosample-latlon'} : "not collected";
	$sample_string[9] = ($opt{'metadata-sra-biosample-host'})? $opt{'metadata-sra-biosample-host'} : "not collected";
	$sample_string[10] = "Severe acute respiratory syndrome";
	$sample_string[11] = ($opt{'metadata-sra-biosample-status'})? $opt{'metadata-sra-biosample-status'} : "not collected";
	$sample_string[12] = ($opt{'metadata-sra-biosample-age'})? $opt{'metadata-sra-biosample-age'} : "not collected";
	$sample_string[13] = ($opt{'metadata-sra-biosample-gender'})? $opt{'metadata-sra-biosample-gender'} : "not collected";
	$sample_string[14] = ($opt{'metadata-sra-biosample-passage'})? $opt{'metadata-sra-biosample-passage'} : "not collected";
	$sample_string[15] = "not collected";
	$sample_string[16] = ($opt{'metadata-sra-biosample-purpose'})? $opt{'metadata-sra-biosample-purpose'} : "not collected";
	$sample_string[17] = ($opt{'metadata-sra-meta-sequencing-purpose'})? $opt{'metadata-sra-meta-sequencing-purpose'} : "not collected";
	$sample_string[18] = ($opt{'metadata-sra-biosample-gisaid'})? $opt{'metadata-sra-biosample-gisaid'} : "not collected";
	if ($opt{'metadata-sra-bioproject-sw'}) {
		$sample_string[19] = ($opt{'metadata-sra-bioproject-id'})? $opt{'metadata-sra-bioproject-id'} : "not collected";
	}
	
	if ($index==0){
		print $ofh "#Pathogen.cl.1.0\n";
		print $ofh join("\t",@header),"\n";
	}
	print $ofh join("\t",@sample_string), "\n";
	close $ofh;
}

sub writeExperiment{
	my $outdir = shift;
	my $outfile = shift;
	my $outfilelist = shift;
	my $input_pe_fastq = shift;
	my $input_se_fastq = shift;
	my $index = shift;
	my $append = shift;
	my @sample_string;
	my @file_string;
	if (!$append){ unlink $outfile; unlink $outfilelist;}
	open my $ofh,  ">>", "$outfile";
	open my $oflh, ">>", "$outfilelist";
	my @header = ("sample_name","library_ID","title","library_strategy","library_source",
				  "library_selection","library_layout","platform","instrument_model",
				  "design_description","filetype","filename","filename2","filename3",
				  "filename4");
	
	$sample_string[0] = $opt{'metadata-sra-biosample-name'} ;
	$sample_string[1] = $opt{'metadata-sra-biosample-name'}."_lib";
	$sample_string[2] = ($opt{'metadata-sra-meta-title'})? $opt{'metadata-sra-meta-title'} : "not collected";
	$sample_string[3] = $opt{'metadata-sra-meta-libstrategy'};
	$sample_string[4] = $opt{'metadata-sra-meta-libsource'};
	$sample_string[5] = $opt{'metadata-sra-meta-libselection'};
	$sample_string[6] = $opt{'metadata-sra-meta-liblayout'};
	$sample_string[7] = $opt{'metadata-sra-meta-platform'};
	$sample_string[8] = $opt{'metadata-sra-meta-libmodel'};
	$sample_string[9] = ($opt{'metadata-sra-meta-design-description'})? $opt{'metadata-sra-meta-design-description'} : "not collected";
	$sample_string[10] = "fastq";
	foreach my $file(@$input_pe_fastq){
		push @sample_string, basename($file);
		push @file_string, basename($file);
	}
	foreach my $file(@$input_se_fastq){
		push @sample_string, basename($file);
		push @file_string, basename($file);
	}
	while ( (15 - scalar(@sample_string)) > 0){
		push @sample_string,"";
	}
	if ($index==0){
		print $ofh join("\t",@header),"\n";
	}
	print $ofh join("\t",@sample_string), "\n";
	print $oflh join(",",$sample_string[0],@file_string),"\n";
	close $ofh;
	close $oflh;
}

sub SetSubmitScript {
	my $submit_script = shift;
	my $submit_log = shift;
	my $projCompleteReport_cache = shift;
	my $sra_done = shift;
	my $ncbi_cmd = shift;
	open (my $fh,">","$submit_script");
	my $cmd = "cd $edgeui_wwwroot\nexport PYTHONUNBUFFERED=on\n";
	$cmd .= "if [ -f $submit_log ]; then \n";
  
	$cmd .= "if ( ! grep -q 'NCBI submit Completed' $submit_log ) then\n";
	$cmd .= "  $ncbi_cmd  \n";
	$cmd .= "else\n";
	$cmd .= "  echo '# NCBI submit Completed already,'\n\n";
	$cmd .= "fi\n";
	$cmd .= "else\n";
	$cmd .= "  $ncbi_cmd  \n";
	$cmd .= "fi\n";
	$cmd .= "if ( grep -q 'NCBI submit Completed' $submit_log ) then\n";
	if ( ref($sra_done) eq 'ARRAY'){
		$cmd .= "  touch ". join(" ",@$sra_done) . "\n";
	}else{
		$cmd .= "  touch $sra_done\n";
	}
	$cmd .= "fi\n";
	if ( ref($projCompleteReport_cache) eq 'ARRAY'){
		$cmd .= "rm -f ". join(" ",@$projCompleteReport_cache) . "\n";
	}else{
		$cmd .= "rm -f  $projCompleteReport_cache\n";
	}
	$cmd .= "rm -f  $submit_script\n" if (! $debug );
	print $fh "source $EDGE_HOME/thirdParty/Anaconda3/bin/activate base\n";
	print $fh $cmd,"\n";
	close $fh;
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

sub make_dir{
        my $dir=shift;
        make_path($dir,{chmod => 0755,});
        return 0;
}

sub checkParams {
	my $projname = shift;
	$projname .= ": " if ($projname);
	my $index = shift;
	# bioproject
	if ($opt{'metadata-sra-bioproject-sw'}){
		&addMessage("PARAMS", "metadata-sra-bioproject-id", "$projname Bioproject ID required." , $index) unless ( $opt{'metadata-sra-bioproject-id'});
	}else{
		&addMessage("PARAMS", "metadata-sra-bioproject-title", "$projname Bioproject Title required." , $index) unless ( $opt{'metadata-sra-bioproject-title'});
		&addMessage("PARAMS", "metadata-sra-bioproject-desc", "$projname Bioproject Description required." , $index) unless ( $opt{'metadata-sra-bioproject-desc'});
		#&addMessage("PARAMS", "metadata-sra-bioproject-link-desc", "$projname Bioproject External Link Desc required." , $index) unless ( $opt{'metadata-sra-bioproject-link-desc'});
		#&addMessage("PARAMS", "metadata-sra-bioproject-link-url", "$projname Bioproject External Link URL required." , $index) unless ( $opt{'metadata-sra-bioproject-link-url'});
	}
	# biosamples
	&addMessage("PARAMS", "metadata-sra-biosample-name", "$projname Biosample name is required." , $index) unless ( $opt{'metadata-sra-biosample-name'}); 
	&addMessage("PARAMS", "metadata-sra-biosample-isolate", "$projname Biosample isolate name is required." , $index) unless ( $opt{'metadata-sra-biosample-isolate'});
	&addMessage("PARAMS", "metadata-sra-biosample-isolate-source", "$projname Biosample isolate source is required." , $index) unless ( $opt{'metadata-sra-biosample-isolate-source'}); 
	&addMessage("PARAMS", "metadata-sra-biosample-location", "$projname Location is required.", $index) unless ( $opt{'metadata-sra-biosample-location'});
	&addMessage("PARAMS", "metadata-sra-biosample-passage", "$projname Passage is required.", $index) unless ( $opt{'metadata-sra-biosample-passage'});
	&addMessage("PARAMS", "metadata-sra-biosample-collection-date", "$projname Collection date is required." , $index) unless ( $opt{'metadata-sra-biosample-collection-date'}); 
	&addMessage("PARAMS", "metadata-sra-biosample-collect-by", "$projname Collect by is required." , $index) unless ( $opt{'metadata-sra-biosample-collect-by'});
	&addMessage("PARAMS", "metadata-sra-biosample-latlon", "$projname Lat Lon by is required." , $index) unless ( $opt{'metadata-sra-biosample-latlon'}); 
	&addMessage("PARAMS", "metadata-sra-biosample-host", "$projname Host is required.", $index) unless ( $opt{'metadata-sra-biosample-host'});
	&addMessage("PARAMS", "metadata-sra-biosample-gender", "$projname Gender is required.", $index) unless ( $opt{'metadata-sra-biosample-gender'});
	&addMessage("PARAMS", "metadata-sra-biosample-age", "$projname Age is required.", $index) unless ( $opt{'metadata-sra-biosample-age'});
	&addMessage("PARAMS", "metadata-sra-biosample-status", "$projname health state is required.", $index) unless ( $opt{'metadata-sra-biosample-status'});
	&addMessage("PARAMS", "metadata-sra-biosample-purpose", "$projname Purpose of sampling is required.", $index) unless ( $opt{'metadata-sra-biosample-purpose'});
	#&addMessage("PARAMS", "metadata-sra-biosample-gisaid", "$projname GISAID accession is required.", $index) unless ( $opt{'metadata-sra-biosample-gisaid'} );

	# SRA experiments
	&addMessage("PARAMS", "metadata-sra-meta-title", "Experiment Title is required." , $index) unless ( $opt{'metadata-sra-meta-title'});
	&addMessage("PARAMS", "metadata-sra-meta-design-description", "Experiment Design is required." , $index) unless ( $opt{'metadata-sra-meta-design-description'});
	&addMessage("PARAMS", "metadata-sra-meta-libselection", "Library Selection is required." , $index) unless ( $opt{'metadata-sra-meta-libselection'});
	&addMessage("PARAMS", "metadata-sra-meta-libstrategy", "Library Strategy is required." , $index) unless ( $opt{'metadata-sra-meta-libstrategy'});
	&addMessage("PARAMS", "metadata-sra-meta-liblayout", "Library Layout is required." , $index) unless ( $opt{'metadata-sra-meta-liblayout'});
	&addMessage("PARAMS", "metadata-sra-meta-libsource", "Library Source is required." , $index) unless ( $opt{'metadata-sra-meta-libsource'});
	&addMessage("PARAMS", "metadata-sra-meta-platform", "Platform is required." , $index) unless ( $opt{'metadata-sra-meta-platform'});
	&addMessage("PARAMS", "metadata-sra-meta-libmodel", "Library Model is required." , $index) unless ( $opt{'metadata-sra-meta-libmodel'});
	&addMessage("PARAMS", "metadata-sra-meta-sequencing-purpose", "Purpose of sequencing is required." , $index) unless ( $opt{'metadata-sra-meta-sequencing-purpose'});
	
	if ($action eq 'upload2sra' or $action eq 'batch-upload2sra'){
		my $today_str = strftime "%Y-%m-%d", localtime;
		&addMessage("PARAMS", "metadata-sra-submitter", "Submitter is required.") unless ( $opt{'metadata-sra-submitter'});
		&addMessage("PARAMS", "metadata-sra-release-date", "Release date is required.") unless ( $opt{'metadata-sra-release-date'});
		if ($opt{'metadata-sra-release-date'} lt $today_str){
			&addMessage("PARAMS", "metadata-sra-release-date", "Release date cannot be the past date.");
		}
	}
}
