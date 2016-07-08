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
	my $metadata = "$out_dir/sample_metadata.txt";
	$vars->{SMD_TYPE_H} = 1; #default setting
	$vars->{SMD_GENDER_M} = 1; #default setting

	if(-e $metadata) {
        	open CONF, $metadata or die "Can't open $metadata $!";
        	while(<CONF>){
      			chomp;
               	 	next if(/^#/);
           		if ( /(.*)=(.*)/ ){
             			$vars->{SMD_TYPE} =$2 if ($1 eq "type");
             			$vars->{SMD_GENDER} =$2 if ($1 eq "gender");
             			$vars->{SMD_AGE} =$2 if ($1 eq "age");
             			$vars->{SMD_HOST} =$2 if ($1 eq "host");
             			$vars->{SMD_HOST_CONDITION} =$2 if ($1 eq "host_condition");
             			$vars->{SMD_SOURCE} =$2 if ($1 eq "source");
             			$vars->{SMD_SOURCE_DETAIL} =$2 if ($1 eq "source_detail");
             			$vars->{SMD_COLLECTION_DATE} =$2 if ($1 eq "collection_date");
             			$vars->{SMD_CITY} =$2 if ($1 eq "city");
             			$vars->{SMD_STATE} =$2 if ($1 eq "state");
             			$vars->{SMD_COUNTRY} =$2 if ($1 eq "country");
             			$vars->{SMD_LAT} =$2 if ($1 eq "lat");
             			$vars->{SMD_LNG} =$2 if ($1 eq "lng");
             			$vars->{SMD_SEQ_PLATFORM} =$2 if ($1 eq "seq_platform");
             			$vars->{SMD_SEQUENCER} =$2 if ($1 eq "sequencer");
             			$vars->{SMD_SEQ_DATE} =$2 if ($1 eq "seq_date");
             			$vars->{SMD_ID} =$2 if ($1 eq "id");
              		}
      		  }
        	close CONF;

		#get drop down options
		if($vars->{SMD_TYPE}) {
			my $sampleType = $vars->{SMD_TYPE};
			if($sampleType eq "human") {
				$vars->{SMD_TYPE_H} = 1;
			} elsif($sampleType eq "animal") {
				$vars->{SMD_TYPE_A} = 1;
				$vars->{SMD_TYPE_H} = 0;
			} elsif($sampleType eq "environmental") {
				$vars->{SMD_TYPE_H} = 0;
				$vars->{SMD_TYPE_E} = 1;
			} 
		}
		if($vars->{SMD_GENDER} ) {
			my $gender = $vars->{SMD_GENDER};
			if($gender eq "male") {
				$vars->{SMD_GENDER_M} = 1;
			} else {
				$vars->{SMD_GENDER_M} = 0;
				$vars->{SMD_GENDER_F} = 1;
			}
		}
		if($vars->{SMD_HOST_CONDITION}) {
			my $val = $vars->{SMD_HOST_CONDITION};
			if($val eq "healthy") {
				$vars->{SMD_HOST_CONDITION_H} = 1;
			} elsif($val eq "diseased") {
				$vars->{SMD_HOST_CONDITION_D} = 1;
			} elsif($val eq "unknown") {
				$vars->{SMD_HOST_CONDITION_U} = 1;
			}
		}
		if($vars->{SMD_SOURCE}) {
			my $val = $vars->{SMD_SOURCE};
			if($val eq "blood") {
				$vars->{SMD_SOURCE_blood} = 1;
			} elsif($val eq "nasal") {
				$vars->{SMD_SOURCE_nasal} = 1;
			} elsif($val eq "saliva") {
				$vars->{SMD_SOURCE_saliva} = 1;
			} elsif($val eq "skin") {
				$vars->{SMD_SOURCE_skin} = 1;
			} elsif($val eq "sputum") {
				$vars->{SMD_SOURCE_sputum} = 1;
			} elsif($val eq "stool") {
				$vars->{SMD_SOURCE_stool} = 1;
			} elsif($val eq "throat") {
				$vars->{SMD_SOURCE_throat} = 1;
			} elsif($val eq "vaginal") {
				$vars->{SMD_SOURCE_vaginal} = 1;
			} elsif($val eq "wound") {
				$vars->{SMD_SOURCE_wound} = 1;
			} elsif($val eq "other") {
				$vars->{SMD_SOURCE_other} = 1;
			} elsif($val eq "unknown") {
				$vars->{SMD_SOURCE_unknown} = 1;
			} elsif($val eq "air") {
				$vars->{SMD_SOURCE_air} = 1;
			} elsif($val eq "built-environment") {
				$vars->{SMD_SOURCE_be} = 1;
			} elsif($val eq "microbial mat/biofilm") {
				$vars->{SMD_SOURCE_mb} = 1;
			} elsif($val eq "plant") {
				$vars->{SMD_SOURCE_plant} = 1;
			} elsif($val eq "sediment") {
				$vars->{SMD_SOURCE_sediment} = 1;
			} elsif($val eq "soil") {
				$vars->{SMD_SOURCE_soil} = 1;
			} elsif($val eq "water") {
				$vars->{SMD_SOURCE_water} = 1;
			} elsif($val eq "wastewater/sludge") {
				$vars->{SMD_SOURCE_ws} = 1;
			} 
		}
		if($vars->{SMD_SEQ_PLATFORM}) {
			my $val = $vars->{SMD_SEQ_PLATFORM};
			if($val eq "Illumina") {
				$vars->{SMD_SEQ_PLATFORM_ILL} = 1;
			} elsif($val eq "IonTorrent") {
				$vars->{SMD_SEQ_PLATFORM_ION} = 1;
			}  elsif($val eq "Nanopore") {
				$vars->{SMD_SEQ_PLATFORM_NAN} = 1;
			}  elsif($val eq "PacBio") {
				$vars->{SMD_SEQ_PLATFORM_PAC} = 1;
			} 
		}
		if($vars->{SMD_SEQUENCER}) {
			my $val = $vars->{SMD_SEQUENCER};
			if($val eq "HiSeq") {
				$vars->{SMD_SEQUENCER_ILL_Hi} = 1;
			} elsif($val eq "HiSeq X") {
				$vars->{SMD_SEQUENCER_ILL_HiX} = 1;
			} elsif($val eq "MiniSeq") {
				$vars->{SMD_SEQUENCER_ILL_Min} = 1;
			} elsif($val eq "MiSeq") {
				$vars->{SMD_SEQUENCER_ILL_Mi} = 1;
			} elsif($val eq "NextSeq") {
				$vars->{SMD_SEQUENCER_ILL_Next} = 1;
			} elsif($val eq "Ion S5") {
				$vars->{SMD_SEQUENCER_ION_S5} = 1;
			} elsif($val eq "Ion PGM") {
				$vars->{SMD_SEQUENCER_ION_PGM} = 1;
			} elsif($val eq "Ion Proton") {
				$vars->{SMD_SEQUENCER_ION_Proton} = 1;
			} elsif($val eq "RS II") {
				$vars->{SMD_SEQUENCER_PAC_RS} = 1;
			} elsif($val eq "Sequel") {
				$vars->{SMD_SEQUENCER_PAC_Sequel} = 1;
			} elsif($val eq "MinTon") {
				$vars->{SMD_SEQUENCER_NAN_MinIon} = 1;
			}
		}
	} 
}

