#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($RealBin);
use File::Basename;
use lib "$RealBin/../../lib";
use HTML::Template;
use POSIX qw{strftime};

my $out_dir       = $ARGV[0];
my $html_outfile  = $ARGV[1];
my $projname = $ARGV[2];
my $userDir = $ARGV[3];
my @out_dir_parts = split('/', $out_dir);
my $projid = $out_dir_parts[-1];

## Instantiate the variables
my $vars;

eval {
	&pull_sampleMetadata();
	&pull_submissionData();
	&pull_consensusInfo();
};

output_html();

sub output_html {
	$vars->{OUTPUTDIR}  = $out_dir;
	$vars->{PROJNAME} ||= $projname;
	$vars->{PROJID} ||= $projid;
	
	my $template = HTML::Template->new(filename => "$RealBin/edge_gisaid_upload.tmpl",
		                               strict => 0,
								       die_on_bad_params => 0);
	$template->param(%$vars);

	system("mkdir -p $out_dir/"."GISAID");

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
}
sub pull_consensusInfo{
	# coverage info is at ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt
	open (my $cov_fh, "<", "$out_dir/ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt");
	my $cov_string;
	while(<$cov_fh>){
		chomp;  
		next if (/^Ref/);
		my @array = split(/\t/,$_);
		if( scalar @array > 8){
			my $id=$array[0];
			my $linear_cov=sprintf("%.2f%%",$array[4]);
			my $depth_cov=sprintf("%dX",$array[5]);
			my $value="$id"."::"."$linear_cov"."::"."$depth_cov";
			$vars->{CON_LIST} .= "<option value=$value>$id ($linear_cov, $depth_cov)</option>";
		}
	}       
	close $cov_fh;
	open (my $log_fh,"<", "$out_dir/ReadsBasedAnalysis/readsMappingToRef/mapping.log");
	my ($tool_version, $align_tool);
	while(<$log_fh>){
		if (/Version:\s+(\S+)/){$tool_version=$1;}
		if (/CMD:\s(\S+)/){$align_tool=$1;}
	}
	close $log_fh;
	$vars->{ASM_METHOD}="$align_tool $tool_version";
}
sub pull_sampleMetadata {
	my $metadata = "$out_dir/metadata_gisaid.txt";
	if(-e $metadata) {
        	open CONF, $metadata or die "Can't open $metadata $!";
        	while(<CONF>){
      			chomp;
               	 	next if(/^#/);
           		if ( /(.*)=(.*)/ ){
             			$vars->{VIR_NAME} =$2 if ($1 eq "virus_name");
             			$vars->{VIR_PASSAGE} =$2 if ($1 eq "virus_passage");
             			$vars->{SM_CDATE} =$2 if ($1 eq "collection_date");
             			$vars->{SM_LOC} =$2 if ($1 eq "location");
             			$vars->{SM_HOST} =$2 if ($1 eq "host");
             			$vars->{SM_GENDER} =$2 if ($1 eq "gender");
             			$vars->{SM_AGE} =$2 if ($1 eq "age");
              		}
      		  }
        	close CONF;
	}
}

sub pull_submissionData {
	my $metadata = "$userDir/gisaid_submission_profile.txt";
	if(-e $metadata) {
		open CONF, $metadata or die "Can't open $metadata $!";
        	while(<CONF>){
      			chomp;
               	 	next if(/^#/);
           		if ( /(.*)=(.*)/ ){
             			$vars->{ORIG_LAB} =$2 if ($1 eq "originating_lab");
             			$vars->{ORIG_ADDRESS} =$2 if ($1 eq "originating_address");
             			$vars->{SUB_LAB} =$2 if ($1 eq "submitting_lab");
             			$vars->{SUB_ADDRESS} =$2 if ($1 eq "submitting_address");
             			$vars->{AUTHORS} =$2 if ($1 eq "authors");
             			$vars->{SUBMITTER} =$2 if ($1 eq "submitter");
             			$vars->{ID} =$2 if ($1 eq "gisaid_id");
              		}
      		  }
        	close CONF;
	} 
}
