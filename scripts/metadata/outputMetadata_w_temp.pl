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
my @out_dir_parts = split('/', $out_dir);
my $projid = $out_dir_parts[-1];

## Instantiate the variables
my $vars;

eval {
	&pull_sampleMetadata();
	&pull_otherMetadata();
};

output_html();

sub output_html {
	$vars->{OUTPUTDIR}  = $out_dir;
	$vars->{PROJNAME} ||= $projname;
	$vars->{PROJID} ||= $projid;
	
	my $template = HTML::Template->new(filename => "$RealBin/edge_html_metadata.tmpl",
		                               strict => 0,
								       die_on_bad_params => 0);
	$template->param(%$vars);

	system("mkdir -p $out_dir/"."Metadata");

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
}

sub pull_sampleMetadata {
	my $metadata = "$out_dir/metadata_sample.txt";
	if(-e $metadata) {
        	open CONF, $metadata or die "Can't open $metadata $!";
        	while(<CONF>){
      			chomp;
               	 	next if(/^#/);
           		if ( /(.*)=(.*)/ ){
             			$vars->{SMD_STUDY_TITLE} =$2 if ($1 eq "study_title");
             			$vars->{SMD_STUDY_TYPE} =$2 if ($1 eq "study_type");
             			$vars->{SMD_NAME} =$2 if ($1 eq "sample_name");
             			$vars->{SMD_TYPE} =$2 if ($1 eq "sample_type");
             			$vars->{OUT_SMD_GENDER} =1 if ($1 eq "gender");
             			$vars->{SMD_GENDER} =$2 if ($1 eq "gender");
             			$vars->{OUT_SMD_AGE} =1 if ($1 eq "age");
             			$vars->{SMD_AGE} =$2 if ($1 eq "age");
             			$vars->{OUT_SMD_HOST} =1 if ($1 eq "host");
             			$vars->{SMD_HOST} =$2 if ($1 eq "host");
             			$vars->{OUT_SMD_HOST_CONDITION} =1 if ($1 eq "host_condition");
             			$vars->{SMD_HOST_CONDITION} =$2 if ($1 eq "host_condition");
             			$vars->{SMD_SOURCE} =$2 if ($1 eq "isolation_source");
             			$vars->{SMD_COLLECTION_DATE} =$2 if ($1 eq "collection_date");
             			$vars->{SMD_LOCATION} =$2 if ($1 eq "location");
             			$vars->{SMD_CITY} =$2 if ($1 eq "city");
             			$vars->{SMD_STATE} =$2 if ($1 eq "state");
             			$vars->{SMD_COUNTRY} =$2 if ($1 eq "country");
             			$vars->{SMD_LAT} =$2 if ($1 eq "lat");
             			$vars->{SMD_LNG} =$2 if ($1 eq "lng");
             			$vars->{SMD_EXP_TITLE} =$2 if ($1 eq "experiment_title");
             			$vars->{SMD_SEQ_CENTER} =$2 if ($1 eq "sequencing_center");
             			$vars->{SMD_SEQUENCER} =$2 if ($1 eq "sequencer");
             			$vars->{SMD_SEQ_DATE} =$2 if ($1 eq "sequencing_date");
              		}
      		  }
        	close CONF;

		if($vars->{SMD_TYPE} && $vars->{SMD_TYPE} eq "human") {
			$vars->{SMD_HOST} = "";
		}
		#get list options
		if($vars->{SMD_TYPE}) {
			my $type = $vars->{SMD_TYPE};
			setSampleType($type);
		}
		if($vars->{SMD_GENDER} ) {
			my $gender = $vars->{SMD_GENDER};
			setGender($gender);
		} else {
			setGender("male");
		}
		if($vars->{SMD_HOST_CONDITION}) {
			my $type = $vars->{SMD_HOST_CONDITION};
			setHostCondition($type);
		} else {
			setHostCondition("unknown");
		}
	} else {
		setSampleType("human");
		setGender("male");
		setHostCondition("unknown");
	}
}

sub setSampleType {
	my $type = shift;
	my @types=('animal', 'environmental', 'human');
			
	for (my $i=0; $i<@types; $i++) {
		my $item;
		$item->{SMD_TYPE} = $types[$i];
		$item->{SMD_TYPE_LABEL} = ucfirst($types[$i]);
		$item->{SMD_TYPE_ID} = "metadata-sample-type$i";
		if($type eq $types[$i]) {
			$item->{SMD_TYPE_CHECKED} = "checked";
		} 
		push @{$vars->{LOOP_SMD_TYPES}}, $item;
	}
}

sub setGender {
	my $gender = shift;
	my @genders = ('male', 'female');
			
	for (my $i=0; $i<@genders; $i++) {
		my $item;
		$item->{SMD_GENDER} = $genders[$i];
		$item->{SMD_GENDER_LABEL} = ucfirst($genders[$i]);
		$item->{SMD_GENDER_ID} = "metadata-host-gender$i";
		if($gender eq $genders[$i]) {
			$item->{SMD_GENDER_CHECKED} = "checked";
		} 
		push @{$vars->{LOOP_SMD_GENDERS}}, $item;
	}
}

sub setHostCondition {
	my $type = shift;
	my @types=('healthy', 'diseased', 'unknown');
			
	for (my $i=0; $i<@types; $i++) {
		my $item;
		$item->{SMD_HOST_CONDITION} = $types[$i];
		$item->{SMD_HOST_CONDITION_LABEL} = ucfirst($types[$i]);
		$item->{SMD_HOST_CONDITION_ID} = "metadata-host-condition$i";
		if($type eq $types[$i]) {
			$item->{SMD_HOST_CONDITION_CHECKED} = "checked";
		} 
		push @{$vars->{LOOP_SMD_HOST_CONDITIONS}}, $item;
	}
}

sub pull_otherMetadata {
	my $metadata = "$out_dir/metadata_other.txt";
	if(-e $metadata) {
		open my $fh, '<', $metadata or die "error opening $metadata: $!";
		$vars->{SMD_OTHER}= do { local $/; <$fh> };	
		close $fh;
	} 
}

