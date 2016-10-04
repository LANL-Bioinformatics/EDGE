#!/usr/bin/env perl

use strict;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
use JSON;
require "edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
my $sid         = $opt{'sid'};
my $prealname = $opt{projname};
my $action = $opt{'action'};

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

my $info;
if ($action eq "check"){
	if ($sys->{edge_sample_metadata}){
		$info->{metadata} = "on";
	}else{
		$info->{metadata}  = "off ";
	}

	&returnStatus();
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
		my $metadata_out = "$projDir/sample_metadata.txt";
		
		#get bsveid
		my $bsveId;
		if( -e $metadata_out ){
			$bsveId = `grep -a "bsve_id=" $metadata_out | awk -F'=' '{print \$2}'`;
			chomp $bsveId;
		}
		#save changes 
		if(open OUT,  ">$metadata_out") {
			print OUT "type=".$opt{'edge-sample-type-edit'}."\n" if ( $opt{'edge-sample-type-edit'} ); 
			if( $opt{'edge-sample-type-edit'} eq "human") {
				if($opt{'pg-cb-gender-edit'}) {
					print OUT "gender=".$opt{'edge-pg-gender-edit'}."\n";
				}
				if($opt{'pg-cb-age-edit'}) {
					print OUT "age=".$opt{'edge-pg-age-edit'}."\n";
				}
			}

			if( $opt{'edge-sample-type-edit'} eq "human" || $opt{'edge-sample-type-edit'} eq "animal") {
				print OUT "host=".$opt{'edge-pg-host-edit'}."\n";
				print OUT "host_condition=".$opt{'edge-pg-host-condition-edit'}."\n";
				print OUT "source=".$opt{'edge-sample-source-host-edit'}."\n";
			} else {
				print OUT "source=".$opt{'edge-sample-source-nonhost-edit'}."\n";
			}
			print OUT "source_detail=".$opt{'edge-pg-sample-source-detail-edit'}."\n";
			print OUT "collection_date=".$opt{'edge-pg-collection-date-edit'}."\n" if ( $opt{'edge-pg-collection-date-edit'} );
			print OUT "city=".$opt{'locality'}."\n" if ( $opt{'locality'} );
			print OUT "state=".$opt{'administrative_area_level_1'}."\n" if ( $opt{'administrative_area_level_1'} );
			print OUT "country=".$opt{'country'}."\n" if ( $opt{'country'} );
			print OUT "lat=".$opt{'lat'}."\n" if ( $opt{'lat'} );
			print OUT "lng=".$opt{'lng'}."\n" if ( $opt{'lng'} );
			print OUT "seq_platform=".$opt{'edge-pg-seq-platform-edit'}."\n" if ( $opt{'edge-pg-seq-platform-edit'} );

			if($opt{'edge-pg-seq-platform-edit'} eq "Illumina") {
				print OUT "sequencer=".$opt{'edge-pg-sequencer-ill-edit'}."\n";
			} elsif($opt{'edge-pg-seq-platform-edit'} eq "IonTorrent") {
				print OUT "sequencer=".$opt{'edge-pg-sequencer-ion-edit'}."\n";
			} elsif($opt{'edge-pg-seq-platform-edit'} eq "Nanopore-edit") {
				print OUT "sequencer=".$opt{'edge-pg-sequencer-nan-edit'}."\n";
			} elsif($opt{'edge-pg-seq-platform-edit'} eq "PacBio") {
				print OUT "sequencer=".$opt{'edge-pg-sequencer-pac-edit'}."\n";
			}

			print OUT "seq_date=".$opt{'edge-pg-seq-date-edit'}."\n" if ( $opt{'edge-pg-seq-date-edit'} );
		
			print OUT "bsve_id=".$bsveId."\n" if ($bsveId);

			close OUT;
		} else {
			$msg->{SUBMISSION_STATUS}="failure";

		}

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


sub returnStatus {
    my $json;
    $json = to_json( $info ) if $info;
    $json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
    print $cgi->header('application/json'), $json;
    exit;
}
