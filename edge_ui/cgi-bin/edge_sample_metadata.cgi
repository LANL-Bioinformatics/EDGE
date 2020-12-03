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
my $info;

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

&stringSanitization(\%opt);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $um_url      = $sys->{edge_user_management_url};


if ($action eq "check"){
	if ($sys->{edge_sample_metadata}){
		$info->{metadata} = "on";
	}else{
		$info->{metadata}  = "off ";
	}

	&returnStatus();
}

if($action eq "edgesite-save") {
	my $msg;
	$msg->{SUBMISSION_STATUS}="success";# return json
	# session check
	if( $sys->{user_management} ){
		my $valid = verifySession($sid);
		if(!$valid){
			$msg->{SUBMISSION_STATUS}="failure";
			$msg->{ERROR} .="Invalid session. Please login again!\n";

			my $json = encode_json($msg);
			print $cgi->header('application/json'), $json;
			exit;
		}else{
			my ($username,$password) = getCredentialsFromSession($sid);
			my $user_info=&getUserInfo($username,$password,$msg);
			$userType=$user_info->{type};
		}
		if($userType ne "admin") {
			
			$msg->{SUBMISSION_STATUS}="failure";
			$msg->{ERROR} .="Permission denied. Only admin can submit this form!\n";

			my $json = encode_json($msg);
			print $cgi->header('application/json'), $json;
			exit;
		}
	}
	if(!$opt{'edgesite-org'} || !$opt{'edgesite-location'}) {
		$msg->{SUBMISSION_STATUS}="failure";
		$msg->{ERROR} .="Organization Full Name and Location must be filled!\n";

		my $json = encode_json($msg);
		print $cgi->header('application/json'), $json;
		exit;
	}
	if(-w $sysconfig) {
		my ($apiKey, $apiToken);

		if( $opt{'edgesite-share-metadata'} eq "yes") {
			$sys->{'edgesite-organization'} = $opt{'edgesite-org'};
			$sys->{'edgesite-acronym'} = $opt{'edgesite-acronym'};
			$sys->{'edgesite-location'} = $opt{'edgesite-location'};
			$sys->{'edgesite-city'} = $opt{'locality'};
			$sys->{'edgesite-state'} = $opt{'administrative_area_level_1'};
			$sys->{'edgesite-country'} = $opt{'country'};
			$sys->{'edgesite-lat'} = $opt{'lat'};
			$sys->{'edgesite-lng'} = $opt{'lng'};

			#get api key/token from BSVE API server
			($apiKey, $apiToken) = pushEdgeSite($sys);
		}

		if($apiKey eq "error") {
			$msg->{SUBMISSION_STATUS}="failure";
			$msg->{ERROR} .="Failed to submit to BSV API server:$apiToken\n";
		} else {
			open OUT,  ">>$sysconfig";
			print OUT "##edgesite settings####\n";
			print OUT "edgesite-organization=".$opt{'edgesite-org'}."\n";
			print OUT "edgesite-acronym=".$opt{'edgesite-acronym'}."\n";
			print OUT "edgesite-location=".$opt{'edgesite-location'}."\n";
			print OUT "edgesite-city=".$opt{'locality'}."\n";
			print OUT "edgesite-state=".$opt{'administrative_area_level_1'}."\n";
			print OUT "edgesite-country=".$opt{'country'}."\n";
			print OUT "edgesite-lat=".$opt{'lat'}."\n";
			print OUT "edgesite-lng=".$opt{'lng'}."\n";
		
			if( $opt{'edgesite-enable-metadata'} eq "yes") {
				print OUT "edge_sample_metadata=1\n";
				if( $opt{'edgesite-share-metadata'} eq "yes") {
					print OUT "edge_sample_metadata_share2bsve=1\n";
					if( $opt{'edgesite-autoshare-metadata'} eq "yes") {
						print OUT "edge_sample_metadata_autosubmit2bsve=1\n";
					} else {
						print OUT "edge_sample_metadata_autosubmit2bsve=0\n";
					}
				} else {
					print OUT "edge_sample_metadata_share2bsve=0\n";
				}
			} else {
				print OUT "edge_sample_metadata=0\n";
			}

			print OUT "bsve_api_key=$apiKey\n" if $apiKey;
			print OUT "bsve_api_token=$apiToken\n" if $apiToken;
		
			print OUT "##END edgesite settings####\n\n";
			close OUT;
		}
	} else {
		$msg->{SUBMISSION_STATUS}="failure";
		$msg->{ERROR} .="Have no write permission on $sysconfig.\n";
		
	}

	if($msg->{SUBMISSION_STATUS} ne "failure") {
		my $dir = dirname($sysconfig);
		my $doneFile = "$dir/edgesite.installation.done";
		`touch $doneFile`;
		if(!-e $doneFile) {
			$msg->{SUBMISSION_STATUS}="failure";
			$msg->{ERROR}.="Have no write permission to create a new file $doneFile.\n";
		}
	}

	my $json = encode_json($msg);
	print $cgi->header('application/json'), $json;
	exit;
}


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
my $proj_list=&scanProjToList($edgeui_output);

if( $proj_list->{$pname} ){
	if($action eq "edit") {
		my $projDir = $relpath . "/". $proj_list->{$pname};
		chdir $edgeui_wwwroot;
		my $outHtml = "$projDir/Metadata/update.html";
		my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/metadata/outputMetadata_w_temp.pl $projDir $outHtml $prealname";
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
	} elsif($action eq "save") {
		my $msg;
		$msg->{SUBMISSION_STATUS}="success";

		my $projDir = $edgeui_output . "/". $proj_list->{$pname};

		if(-e $projDir) {
			#travels
			my $travel_out = "$projDir/metadata_travels.txt";
			my @travels = split /[\x0]/, $opt{"metadata-travel-location"};
			my @travel_df = split /[\x0]/, $opt{"metadata-travel-date-f"};
			my @travel_dt = split /[\x0]/, $opt{"metadata-travel-date-t"};
			my @cities = split /[\x0]/, $opt{"locality"};
			my @states = split /[\x0]/, $opt{'administrative_area_level_1'};
			my @countries = split /[\x0]/, $opt{'country'};
			my @lats = split /[\x0]/, $opt{'lat'};
			my @lngs = split /[\x0]/, $opt{'lng'};
			my $geoLoc_cnt = 0;
			if (@travels > 0) {
				open TOUT, ">$travel_out";
				foreach my $travel (@travels) {
					print TOUT "travel-date-from=".$travel_df[$geoLoc_cnt]."\n";
					print TOUT "travel-date-to=".$travel_dt[$geoLoc_cnt]."\n";
					print TOUT "travel-location=$travel\n";
					print TOUT "city=".$cities[$geoLoc_cnt]."\n";
					print TOUT "state=".$states[$geoLoc_cnt]."\n";
					print TOUT "country=".$countries[$geoLoc_cnt]."\n";
					print TOUT "lat=".$lats[$geoLoc_cnt]."\n";
					print TOUT "lng=".$lngs[$geoLoc_cnt]."\n";
					$geoLoc_cnt++;
				}
				close(TOUT);
			} else {
				#delete old travels if exist
				`rm -f $travel_out`;
			}

			#symptoms
			my $symptom_out = "$projDir/metadata_symptoms.txt";
			my @symptom_cats = split /[\x0]/, $opt{"metadata-symptom-cat"};
			my @symptom_catIDs = split /[\x0]/, $opt{"metadata-symptom-catID"};
			my ($cat,$catID);
			my $symptomLines;
			for(my $cat_cnt=0; $cat_cnt<@symptom_cats;$cat_cnt++) {
				$cat = $symptom_cats[$cat_cnt];
				$catID = $symptom_catIDs[$cat_cnt];
				my @symptoms = split /[\x0]/, $opt{"metadata-symptom-$catID"};
				foreach my $symptom(@symptoms) {
					$symptomLines .= "$cat\t$symptom\n";
				}
			}
			if($symptomLines) {
				open SOUT, ">$symptom_out";
				print SOUT "$symptomLines";
				close(SOUT);
			} else {
				#delete old symptoms if exist
				`rm -f $symptom_out`;
			}

			#sample metadata
			my $metadata_out = "$projDir/metadata_sample.txt";
			open OUT,  ">$metadata_out";
			my $id = `perl edge_db.cgi study-title-add "$opt{'metadata-study-title'}"`;
			print OUT "study_id=$id\n";
			print OUT "study_title=".$opt{'metadata-study-title'}."\n" if ( $opt{'metadata-study-title'} );
			`perl edge_db.cgi study-type-add "$opt{'metadata-study-type'}"`;
			print OUT "study_type=".$opt{'metadata-study-type'}."\n" if ( $opt{'metadata-study-type'} ); ; 
			print OUT "sample_name=".$opt{'metadata-sample-name'}."\n" if ( $opt{'metadata-sample-name'} ); 
			print OUT "sample_type=".$opt{'metadata-sample-type'}."\n" if ( $opt{'metadata-sample-type'} ); 

			if( $opt{'metadata-sample-type'} eq "human" || $opt{'metadata-sample-type'} eq "animal") {
				if($opt{'metadata-sample-type'} eq "human") {
					if($opt{'metadata-host-gender-cb'}) {
						print OUT "gender=".$opt{'metadata-host-gender'}."\n";
					}
					if($opt{'metadata-host-age-cb'}) {
						print OUT "age=".$opt{'metadata-host-age'}."\n";
					}
					print OUT "host=human\n";
				} else {
					`perl edge_db.cgi animal-host-add "$opt{'metadata-host'}"`;
					print OUT "host=".$opt{'metadata-host'}."\n";
					`rm -f $travel_out`;
					`rm -f $symptom_out`;
				}
				print OUT "host_condition=".$opt{'metadata-host-condition'}."\n";
				`perl edge_db.cgi isolation-source-add "$opt{'metadata-isolation-source'}"`;
				print OUT "isolation_source=".$opt{'metadata-isolation-source'}."\n";
			} else {
				`perl edge_db.cgi isolation-source-add "$opt{'metadata-isolation-source'}" "environmental"`;
				print OUT "isolation_source=".$opt{'metadata-isolation-source'}."\n";
				`rm -f $travel_out`;
				`rm -f $symptom_out`;
			}
			
			print OUT "collection_date=".$opt{'metadata-sample-collection-date'}."\n" if ( $opt{'metadata-sample-collection-date'} );
			print OUT "location=".$opt{'metadata-sample-location'}."\n" if ( $opt{'metadata-sample-location'} );
			print OUT "city=".$cities[$geoLoc_cnt]."\n" if($cities[$geoLoc_cnt]);
			print OUT "state=".$states[$geoLoc_cnt]."\n" if($states[$geoLoc_cnt]);
			print OUT "country=".$countries[$geoLoc_cnt]."\n" if($countries[$geoLoc_cnt]);
			print OUT "lat=".$lats[$geoLoc_cnt]."\n" if($lats[$geoLoc_cnt]);
			print OUT "lng=".$lngs[$geoLoc_cnt]."\n" if($lngs[$geoLoc_cnt]);
			print OUT "experiment_title=".$opt{'metadata-exp-title'}."\n";
			`perl edge_db.cgi seq-center-add "$opt{'metadata-seq-center'}"`;
			print OUT "sequencing_center=".$opt{'metadata-seq-center'}."\n";
			`perl edge_db.cgi sequencer-add "$opt{'metadata-sequencer'}"`;
			print OUT "sequencer=".$opt{'metadata-sequencer'}."\n";
			print OUT "sequencing_date=".$opt{'metadata-seq-date'}."\n" if ( $opt{'metadata-seq-date'} );
			close OUT;
			
			#user defined metadata
			if($opt{'metadata-other'}) {
				my $other_out = "$projDir/metadata_other.txt";
				open OUT,  ">$other_out";
				print OUT $opt{'metadata-other'};
				close OUT;
			}
		} else {
			$msg->{SUBMISSION_STATUS}="failure";
		}

		#delete HTML_Report/.complete_report_web
		`rm $projDir/HTML_Report/.complete_report_web`;

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
