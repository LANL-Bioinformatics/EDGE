#!/usr/bin/env perl

use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use JSON;
use File::Basename;
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
my $info;
my $msg;



my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $um_url      = $sys->{edge_user_management_url};



# session check
if( $sys->{user_management} ){
	my $valid = verifySession($sid);
	if(!$valid){
		my $html = "Invalid session. Please login again!";
		print "Content-Type: text/html\n\n",
		$html;
		exit;
	}
}

my $domain      = $ENV{'HTTP_HOST'} || 'edge-bsve.lanl.gov';
my ($webhostname) = $domain =~ /^(\S+?)\./;
$sys->{edgeui_output} = "$sys->{edgeui_output}"."/$webhostname" if ( -d "$sys->{edgeui_output}/$webhostname");

my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};
my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;
my ($proj_list)=&scanProjToList($edgeui_output);

if( $proj_list->{$pname} ){
	$prealname=$proj_list->{$pname}->{projname};
	if($action eq "create-form") {
		my $projDir = $relpath . "/". $proj_list->{$pname}->{dir};
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
	} elsif($action eq "upload2gisaid" or $action eq "download") {
		$msg->{SUBMISSION_STATUS}="success";

		# validating parameters
		&checkParams();
		if($msg->{SUBMISSION_STATUS} eq "failure") {
			returnParamsStatus();
		}

		my $projDir = $edgeui_output . "/". $proj_list->{$pname}->{dir};
		my $projDir_rel = $relpath . "/". $proj_list->{$pname}->{dir};
		if(-e $projDir) {
			#sample metadata
			my $metadata_out = "$projDir/metadata_gisaid.txt";
			my $selected_consensus = $opt{'metadata-sample-consensus'};
			$selected_consensus =~ s/::/ /g;
			my ($consensus_ref,$consensus_cov,$consensus_depth) = split(/\s+/,$selected_consensus);
			$consensus_ref =~ s/_?\d?$//g;
			$selected_consensus = "Consensus to ". $selected_consensus;
			open OUT,  ">$metadata_out";
			print OUT "virus_name=".$opt{'metadata-virus-name'}."\n"; 
			print OUT "virus_passage=".$opt{'metadata-virus-passage'}."\n"; 
			print OUT "collection_date=".$opt{'metadata-sample-collection-date'}."\n";
			print OUT "location=".$opt{'metadata-sample-location'}."\n";
			print OUT "host=".$opt{'metadata-sample-host'}."\n";
			print OUT "gender=".$opt{'metadata-sample-gender'}."\n";
			print OUT "age=".$opt{'metadata-sample-age'}."\n";
			print OUT "status=".$opt{'metadata-sample-status'}."\n";
			print OUT "sequencing_technology=".$opt{'metadata-sample-sequencing-tech'}."\n";
			print OUT "coverage=$selected_consensus\n";
			close OUT;
			#gisaid profile
			$metadata_out = $sys->{edgeui_input}."/$userDir/gisaid_submission_profile.txt";
			open OUT,  ">$metadata_out";
			print OUT "originating_lab=".$opt{'metadata-orig-lab'}."\n"; 
			print OUT "originating_address=".$opt{'metadata-orig-address'}."\n"; 
			print OUT "submitting_lab=".$opt{'metadata-sub-lab'}."\n";
			print OUT "submitting_address=".$opt{'metadata-sub-address'}."\n";
			print OUT "authors=".$opt{'metadata-authors'}."\n";
			print OUT "submitter=".$opt{'metadata-submitter'}."\n";
			print OUT "gisaid_id=".$opt{'metadata-gisaid-id'}."\n";
			close OUT;

			#create an input file for submission pipeline
			$metadata_out = "$projDir/gisaid_submission.txt";
			open OUT,  ">$metadata_out";
			print OUT "virus_name=".$opt{'metadata-virus-name'}."\n"; 
			print OUT "virus_passage=".$opt{'metadata-virus-passage'}."\n"; 
			print OUT "collection_date=".$opt{'metadata-sample-collection-date'}."\n";
			print OUT "location=".$opt{'metadata-sample-location'}."\n";
			print OUT "host=".$opt{'metadata-sample-host'}."\n";
			print OUT "gender=".$opt{'metadata-sample-gender'}."\n";
			print OUT "age=".$opt{'metadata-sample-age'}."\n";
			print OUT "status=".$opt{'metadata-sample-status'}."\n";
			print OUT "sequencing_technology=".$opt{'metadata-sample-sequencing-tech'}."\n";
			print OUT "assembly_method=$opt{'metadata-sample-assembly-method'}\n";
			print OUT "originating_lab=".$opt{'metadata-orig-lab'}."\n"; 
			print OUT "originating_address=".$opt{'metadata-orig-address'}."\n"; 
			print OUT "submitting_lab=".$opt{'metadata-sub-lab'}."\n";
			print OUT "submitting_address=".$opt{'metadata-sub-address'}."\n";
			print OUT "authors=".$opt{'metadata-authors'}."\n";
			print OUT "submitter=".$opt{'metadata-submitter'}."\n";
			print OUT "gisaid_id=".$opt{'metadata-gisaid-id'}."\n";
			#print OUT "gisaid_pw=".$opt{'metadata-gisaid-pw'}."\n";
			print OUT "coverage=$selected_consensus\n";
			close OUT;
			my $cmd;
			if ($action eq "download"){
				# download
				my $zip_file = "$projDir/GISAID/${prealname}_gisaid_data.zip";
				my $zip_file_rel = "$projDir_rel/GISAID/${prealname}_gisaid_data.zip";
				$cmd = "zip -j $zip_file $metadata_out $projDir/ReadsBasedAnalysis/readsMappingToRef/${consensus_ref}*_consensus.fasta 2>/dev/null";
				`$cmd`;
				$msg->{PATH} = $zip_file_rel;
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
					`touch $projDir/gisaid_submission.done`;
				}
			}
			if ( ! -e "$projDir/gisaid_submission.done"){
					$msg->{SUBMISSION_STATUS}="failure";
			}
			unlink "$projDir/GISAID/submit.sh";
			
		} else {
			$msg->{SUBMISSION_STATUS}="failure";
		}

		#delete HTML_Report/.complete_report_web
		unlink "$projDir/HTML_Report/.complete_report_web";

		# return json
		my $json = encode_json($msg);
		print $cgi->header('application/json'), $json;
	} else {
		my $msg;
		$msg->{SUBMISSION_STATUS}="failure";
	
		# return json
		my $json = encode_json($msg);
		print $cgi->header('application/json'), $json;
	}
}

exit;

######################################################
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
		$msg->{ERROR}=$1;
                return;
        }
        my $tmp_r = from_json($result_json);
        $user_info{$_} = $tmp_r->{$_} foreach (keys %$tmp_r);
        return \%user_info;
}

sub scanProjToList {
	my $out_dir = shift;
        my $list;
	my $projname;
        opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
        while( defined (my $file = readdir BIN) ) {
                next if $file eq '.' or $file eq '..';
		my $projid;
		my $projCode;
                if ( -d "$out_dir/$file" && -r "$out_dir/$file/config.txt"  ) {
			open ( CONFIG, "$out_dir/$file/config.txt") or die "Cannot open $out_dir/$file/config.txt\n";
			while(<CONFIG>){
				last if (/^\[Down/);
				if (/^projid=(\S+)/){
					$projid=$1;
				}
				if (/^projcode=(\S+)/){
					$projCode=$1;
				}
				if (/^projname=(\S+)/){
					$projname=$1;
				}
			}
			close CONFIG;
                        $projid ||= $file;
			$list->{$projid}->{dir} = $file;
			$list->{$file}->{dir}  = $file;
			$list->{$projCode}->{dir} = $file;
			$list->{$projCode}->{projname} = $projname;
			$list->{$projid}->{projname} = $projname;
			$list->{$file}->{projname}  = $projname;
                }
        }
        closedir(BIN);
        return ($list);
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


sub returnStatus {
    my $json;
    $json = to_json( $info ) if $info;
    $json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
    print $cgi->header('application/json'), $json;
    exit;
}

sub stringSanitization{
	my $opt=shift;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		if ($key =~ /metadata|edgesite|locality|administrative|country|lat|lng/){
			$opt->{$key} =~ s/[`";']/ /g;
			next;
		}
		if($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/){
			$info->{INFO} = "Invalid characters detected \'$str\'.";
			&returnStatus();
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
	my ($job, $status, $note) = @_;
	if( $job eq "PARAMS" ){
		#e.g.: $msg->{"PARAMS"}->{"edge-ref-file"}="File not found";
		$msg->{$job}->{$status}=$note;
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
	&addMessage("PARAMS", "metadata-virus-name", "Virus name is required.") unless ( $opt{'metadata-virus-name'}); 
	&addMessage("PARAMS", "metadata-virus-passage", "Passage details/history is required.") unless ( $opt{'metadata-virus-passage'}); 
	&addMessage("PARAMS", "metadata-sample-collection-date", "Collection date is required.") unless ( $opt{'metadata-sample-collection-date'}); 
	&addMessage("PARAMS", "metadata-sample-location", "Location is required.") unless ( $opt{'metadata-sample-location'}); 
	&addMessage("PARAMS", "metadata-sample-host", "Host is required.") unless ( $opt{'metadata-sample-host'});
	&addMessage("PARAMS", "metadata-sample-gender", "Gender is required.") unless ( $opt{'metadata-sample-gender'});
	&addMessage("PARAMS", "metadata-sample-age", "Age is required.") unless ( $opt{'metadata-sample-age'});
	&addMessage("PARAMS", "metadata-sample-status", "Status is required.") unless ( $opt{'metadata-sample-status'});
	&addMessage("PARAMS", "metadata-sample-sequencing-tech", "Sequencing Technology is required.") unless ( $opt{'metadata-sample-sequencing-tech'});
	&addMessage("PARAMS", "metadata-sample-consensus", "Consensus FASTA is required.") unless ( $opt{'metadata-sample-consensus'} );


	&addMessage("PARAMS", "metadata-orig-lab", "Originating lab is required.") unless ( $opt{'metadata-orig-lab'});
	&addMessage("PARAMS", "metadata-orig-address", "Originating address is required.") unless ( $opt{'metadata-orig-address'});
	&addMessage("PARAMS", "metadata-org-lab", "Submitting lab is required.") unless ( $opt{'metadata-sub-lab'});
	&addMessage("PARAMS", "metadata-sub-address", "Submitting address is required.") unless ( $opt{'metadata-sub-address'});
	&addMessage("PARAMS", "metadata-authors", "Authors is required.") unless ( $opt{'metadata-authors'});
	&addMessage("PARAMS", "metadata-submitter", "Submitter is required.") unless ( $opt{'metadata-submitter'});
	&addMessage("PARAMS", "metadata-gisaid-id", "GISAID id is required.") unless ( $opt{'metadata-gisaid-id'});
	if ($action eq 'upload2gisaid'){
		&addMessage("PARAMS", "metadata-gisaid-pw", "GISAID password is required.") unless ( $opt{'metadata-gisaid-pw'});
	}
}
