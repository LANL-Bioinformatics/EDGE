#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($RealBin);
use File::Basename;
use lib "$RealBin/../../lib";
use HTML::Template;
use POSIX qw{strftime};
use Getopt::Long;
use File::Path;

my $html_outfile;
my $project_dir_names;
my $userDir;
my $usage = qq{
Usage: $0
        Required
                -out             html output file
                -projects        project_dir_names, separated by comma
                -udir            user dir
};

GetOptions(
                "out=s"        =>  \$html_outfile,
                "projects=s"       =>  \$project_dir_names,
                "udir=s"       => \$userDir,
                "help|?"           => sub{print "$usage\n";exit;} 
        );

if (!$project_dir_names && !$html_outfile){ print "$usage\n";exit;}

mkpath(dirname($html_outfile));

my $info;
foreach my $proj_dir (split /,/,$project_dir_names){
	## Instantiate the variables
	my $vars={};
	eval {
		&pull_sampleMetadata($proj_dir,$vars);
		&pull_consensusInfo($proj_dir,$vars);
		&check_submission_status($proj_dir,$vars);
	};
	push @{$info->{LOOP_METADATA}}, $vars;
}

eval{
	&pull_submissionData($info);
};

output_html($info);

sub output_html {
	my $template = HTML::Template->new(filename => "$RealBin/edge_gisaid_projectTable.tmpl",
		                               strict => 0,
								       die_on_bad_params => 0);
	$template->param(%$info);

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
}
sub check_submission_status{
	my $out_dir=shift;
	my $vars = shift;
	my $gisaid_done = "$out_dir/gisaid_submission.done";
	if ( -e $gisaid_done){
		my $gisaid_submit_date = strftime "%F",localtime((stat("$gisaid_done"))[9]);
                $vars->{GISAID_SUBMIT_TIME} = $gisaid_submit_date;
	}
}

sub pull_consensusInfo{
	my $out_dir=shift;
	my $vars = shift;
	my $consensus_length_recommand=25000;
	my $consensus_dpcov_recommand=10;
	my $consensus_Nper_recommand=0.05;
	# coverage info is at ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt
	my $con_comp = "$out_dir/ReadsBasedAnalysis/readsMappingToRef/*consensus.fasta.comp";
	open (my $cov_fh, "<", "$out_dir/ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt");
	my $cov_string;
	while(<$cov_fh>){
		chomp;  
		next if (/^Ref/);
		my @array = split(/\t/,$_);
		if( scalar @array > 8){
			my $id=$array[0];
			my $consense_comp_str =`grep $id $con_comp`;
			my @consense_info = split /\s+/,$consense_comp_str;
			my $linear_cov=sprintf("%.2f%%",$array[4]);
			my $depth_cov=sprintf("%dX",$array[5]);
			my $value="$id"."::"."$linear_cov"."::"."$depth_cov";
			my $selected = ( $vars->{SM_COV} && $vars->{SM_COV} =~ /$id/ )? 'selected':'';
			if ($consense_info[1] >= $consensus_length_recommand && $array[5] >= $consensus_dpcov_recommand && ($consense_info[8]/$consense_info[1]) <= $consensus_Nper_recommand){
				$vars->{CON_LIST} .= "<option value='$value' $selected>$id ($linear_cov, $depth_cov), Ready to Submit</option>";
			}else{
				$vars->{CON_LIST} .= "<option value='$value' $selected>$id ($linear_cov, $depth_cov)</option>";
			}
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
	$vars->{ASM_METHOD}="EDGE-covid19: $align_tool $tool_version.";

	my $confFile = "$out_dir/config.txt";
	my $configuration = &pull_EDGEConfig($confFile);

	my $con_min_mapQ = $configuration->{r2g_consensus_min_mapQ};
	my $con_min_cov = $configuration->{r2g_consensus_min_cov};
	my $con_alt_prop = $configuration->{r2g_consensus_alt_prop}*100 . "%";
	my $con_alt_indel = $configuration->{r2g_consensus_altIndel_prop}* 100 . "%";
	$vars->{PROJNAME} = $configuration->{projname};
	$vars->{PROJCODE} = $configuration->{projcode};
	$vars->{ASM_METHOD} .= " Consensus min coverage: ${con_min_cov}X. min map quality: $con_min_mapQ. Alternate Base > $con_alt_prop. Indel > $con_alt_indel.";
	$vars->{PROCHECKBOX} = "<input type='checkbox' class='edge-metadata-projpage-ckb' name='edge-metadata-projpage-ckb' value='$configuration->{projcode}'>";
	#$vars->{SUMITCRITERIA} = sprintf("Length > %d bp, Depth Coverage > %d X, N < %.1f %%",$consensus_length_recommand,$consensus_dpcov_recommand,$consensus_Nper_recommand*100 );
}
sub pull_sampleMetadata {
	my $out_dir=shift;
	my $vars = shift;
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
             			$vars->{SM_GENDER} = $2 if ($1 eq "gender");
             			$vars->{SM_AGE} =$2 if ($1 eq "age");
             			$vars->{SM_STATUS} =$2 if ($1 eq "status");
             			$vars->{SM_SEQUENCING_TECH} =$2 if ($1 eq "sequencing_technology");
				$vars->{SM_COV} = $2 if ($1 eq "coverage");
              		}
      		}
		if ($vars->{SM_GENDER} =~ /Female/i){
			$vars->{SM_GENDER_FEMALE} = "selected";
		}elsif($vars->{SM_GENDER} =~ /Male/i){
			$vars->{SM_GENDER_MALE} = "selected";
		}elsif($vars->{SM_GENDER} =~ /Unknown/i){
			$vars->{SM_GENDER_UNKNOWN} = "selected";
		}elsif($vars->{SM_GENDER} =~ /Other/i){
			$vars->{SM_GENDER_OTHER} = "selected";
		}
		if ($vars->{SM_STATUS} =~ /Hospitalized/i){
			$vars->{SM_STATUS_Hospitalized} = "selected";
		}elsif($vars->{SM_STATUS} =~ /Released/i){
			$vars->{SM_STATUS_Released} = "selected";
		}elsif($vars->{SM_STATUS} =~ /Live/i){
			$vars->{SM_STATUS_Live} = "selected";
		}elsif($vars->{SM_STATUS} =~ /Deceased/i){
			$vars->{SM_STATUS_Deceased} = "selected";
		}elsif($vars->{SM_STATUS} =~ /Unknown/i){
			$vars->{SM_STATUS_Unknown} = "selected";
		}
        	close CONF;
	}
}

sub pull_submissionData {
	my $vars = shift;
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
sub pull_EDGEConfig
{
    my $file=shift;
    my %hash;
    open (my $fh , $file) or die "No config file $!\n";
    my $head=<$fh>;
    if ($head !~ /project/i){ die "Incorrect config file\n"};
    while (<$fh>)
    {
        chomp;
        next if (/^#/);
        if (/=/)
        {
            my ($key,$value)=split /=/,$_;
            if ( defined $value)
            {
               $value =~ s/\"//g;
               if ($key eq "Host")
               {
                 foreach my $each_host(split(/,/,$value))
                 {
                   push @{$hash{$key}} , $each_host;
                 }
               }
               elsif($key eq "reference")
               {
                 foreach my $each_ref(split(/,/,$value))
                 {
                   push @{$hash{$key}} , $each_ref;
                 }
               }
               else
               {
                   $hash{$key}=$value;
               }
            }
            else
            {
               $hash{$key}="";
            }
        }
    }
    close $fh;
    return \%hash;
}
