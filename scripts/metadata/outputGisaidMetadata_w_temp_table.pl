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
my $index=0;
foreach my $proj_dir (split /,/,$project_dir_names){
	## Instantiate the variables
	my $vars={};
	$vars->{ROWINDEX} = $index;
	$index++;
	my $confFile = "$proj_dir/config.txt";
	my $configuration = &pull_EDGEConfig($confFile);
	eval {
		&pull_biosamples($proj_dir,$vars,$configuration);
		&pull_experiments($proj_dir,$vars);
		&pull_sra_additional($proj_dir,$info,$configuration);
		&pull_sampleMetadata($proj_dir,$vars);
		&pull_consensusInfo($proj_dir,$vars,$configuration);
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
	my $gisaid_done = "$out_dir/UPLOAD/gisaid_ncbi_submission.done";
	if ( -e $gisaid_done){
		my $gisaid_submit_date = strftime "%F",localtime((stat("$gisaid_done"))[9]);
                $vars->{GISAID_SUBMIT_TIME} = $gisaid_submit_date;
	}
}

sub pull_biosamples{
    my $out_dir=shift;
    my $vars=shift;
    my $configuration=shift;
	my $metadata = "$out_dir/UPLOAD/sra_samples.txt";
	if(-e $metadata) {
		open my $fh, $metadata or die "Can't open $metadata $!";
		#my @headers = ("sample_name", "sample_title", "organism", "isolate", "collected_by", "collection_date", 
		#              "geo_loc_name", "isolation_source", "lat_lon", "host", "host_disease", 
		#              "host_health_state","host_age","host_sex",
		#              "passage_history", "description", "purpose_of_sampling",
		#              "purpose_of_sequencing","GISAID_accession","bioproject_accession");
		while(<$fh>){
			chomp;
			next if(/^#/);
			my (@headers, @itmes);
			my @items = split /\t/,$_ ;
			if (/^sample_name/){
				@headers = @items;
			}else{
				$vars->{BIOSAMPLE_NAME} =$items[0];
				$vars->{VIR_NAME} =$items[3];
				$vars->{BIOSAMPLE_CBY} =$items[4];
				$vars->{SM_CDATE} =$items[5];
				$vars->{BIOSAMPLE_LOC} =$items[6];
				$vars->{BIOSAMPLE_ISOLATESOURCE} =$items[7];
				$vars->{BIOSAMPLE_LATLON} =$items[8];
				$vars->{SM_HOST} =$items[9];
				$vars->{"SM_STATUS_"."$items[11]"} ="selected";
				$vars->{SM_AGE} =$items[12];
				my $gender = uc($items[13]);
				$vars->{"SM_GENDER_"."$gender"} ="selected";
				$vars->{VIR_PASSAGE} =$items[14];
				(my $selected_BIOSAMPLE_PS = $items[16]) =~ s/\s/_/g;
				$vars->{"BIOSAMPLE_PS_"."$selected_BIOSAMPLE_PS"} ="selected";
				(my $selected_SRA_PS = $items[17]) =~ s/\s/_/g;
				$vars->{"SRA_PS_"."$selected_SRA_PS"} ="selected";
				$vars->{BIOSAMPLE_GISAIDACC} =$items[18];
				$vars->{SM_BIOPROJECT_ID} = $items[19] if scalar(@items) == 20;
			}
		}
		close $fh;
	}
	$vars->{VIR_NAME} =~ s/SARS-CoV-2\/Homo sapiens/hCoV-19/i;
	$vars->{BIOSAMPLE_NAME} ||= $configuration->{projname};		     	
}

sub pull_experiments{
	my $out_dir=shift;
	my $vars=shift;
	my $metadata = "$out_dir/UPLOAD/sra_experiments.txt";
	my ($selected_SRA_LIBSTRATEGY, $selected_SRA_LIBSOURCE, $selected_SRA_LIBSELECT, $selected_SRA_SEQMODEL)=("","","","");
	if(-e $metadata) {
		open my $fh, $metadata or die "Can't open $metadata $!";
	#my @header = ("sample_name","library_ID","title","library_strategy","library_source",
	#			  "library_selection","library_layout","platform","instrument_model",
	#			  "design_description","filetype","filename","filename2","filename3",
	#			  "filename4");
		
		while(<$fh>){
			chomp;
			next if(/^#/);
			my (@headers, @itmes);
			my @items = split /\t/,$_; 
			if (/^sample_name/){
				@headers = @items;
			}else{
			#  {methodology} of {organism}: {sample info}"
				$vars->{"SRA_LIBTITLE"} =$items[2];
				($selected_SRA_LIBSTRATEGY = $items[3]) =~ s/\s/_/g;
				$vars->{"SRA_LIBSTRATEGY_"."$selected_SRA_LIBSTRATEGY"} = "selected";
				($selected_SRA_LIBSOURCE = $items[4]) =~ s/\s/_/g;
				$vars->{"SRA_LIBSOURCE_"."$selected_SRA_LIBSOURCE"} = "selected";
				($selected_SRA_LIBSELECT = $items[5]) =~ s/\s/_/g;
				$vars->{"SRA_LIBSELECT_"."$selected_SRA_LIBSELECT"} = "selected";
				$vars->{"SRA_LIBLAYOUT_"."$items[6]"} = "selected";
				$vars->{"SRA_PLATFORM_"."$items[7]"} = "selected";
				($selected_SRA_SEQMODEL = $items[8]) =~ s/\s/_/g;
				$vars->{"SRA_SEQMODEL_"."$selected_SRA_SEQMODEL"} = "selected";
				$vars->{"SRA_LIBDESIGN"} =$items[9];
				$vars->{SM_SEQUENCING_TECH} = "$items[7] $selected_SRA_SEQMODEL"; 
			}
		}
		close $fh;
	}	
}
sub pull_sra_additional {
	my $out_dir=shift;
	my $vars=shift;
	my $configuration=shift;
	my $metadata = "$out_dir/UPLOAD/sra_additional_info.txt";
	if(-e $metadata) {
		open my $fh, $metadata or die "Can't open $metadata $!";
		while(<$fh>){
			chomp;
			next if(/^#/);
			if ( /(.*)=(.*)/ ){
				$vars->{SRA_SUBMITTER} =$2 if ($1 eq "metadata-sra-submitter");
				$vars->{SM_RELEASE_DATE} =$2 if ($1 eq "metadata-sra-release-date");
			}
		}
		close $fh
	}
	$vars->{SRA_SUBMITTER} ||= $configuration->{projowner};
	my $today_str = strftime "%Y-%m-%d", localtime;
	$vars->{SM_RELEASE_DATE} ||= $today_str;
}

sub pull_consensusInfo{
	my $out_dir=shift;
	my $vars = shift;
	my $configuration=shift;
	my $consensus_length_recommand=25000;
	my $consensus_dpcov_recommand=10;
	my $consensus_Nper_recommand=0.05;
	# coverage info is at ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt
	my @con_comp_files = glob("$out_dir/ReadsBasedAnalysis/readsMappingToRef/*consensus.fasta.comp");
	my $con_comp_info=&consensus_composition_info(\@con_comp_files);
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
			my $con_fasta=$con_comp_info->{$id}->{file};
			my $con_fasta_prefix = $con_fasta =~ s/.fasta//r;
                        my $con_amb_fasta = $con_fasta_prefix."_w_ambiguous.fasta";
                        my $con_len = $con_comp_info->{$id}->{len} ;
                        my $con_N = $con_comp_info->{$id}->{N} + $con_comp_info->{$id}->{n};
                        my ($selected, $amb_selected)=("","");
                        if ( $vars->{SM_COV} ){
                                if ($vars->{SM_COV} =~ /$con_fasta_prefix/  && $vars->{SM_COV} =~ /w_ambiguous/ ){
                                        $amb_selected = "selected";
                                }elsif  ($vars->{SM_COV} =~ /$con_fasta_prefix/){
                                        $selected = "selected";
                                }
                        }
			if ($con_len >= $consensus_length_recommand && $array[5] >= $consensus_dpcov_recommand && ($con_N/$con_len) <= $consensus_Nper_recommand){
				$vars->{CON_LIST} .= "<option value='${con_fasta}::${linear_cov}::${depth_cov}' $selected>$id ($linear_cov, $depth_cov), Ready to Submit</option>";
				$vars->{CON_LIST} .= "<option value='${con_amb_fasta}::${linear_cov}::${depth_cov}' $amb_selected>$id ($linear_cov, $depth_cov, *With Ambiguous*), Ready to Submit</option>" if ( -r "$out_dir/ReadsBasedAnalysis/readsMappingToRef/$con_amb_fasta");
			}else{
				$vars->{CON_LIST} .= "<option value='${con_fasta}::${linear_cov}::${depth_cov}' $selected>$id ($linear_cov, $depth_cov)</option>";
				$vars->{CON_LIST} .= "<option value='${con_amb_fasta}::${linear_cov}::${depth_cov}' $amb_selected>$id ($linear_cov, $depth_cov, *With Ambiguous*)</option>" if ( -r "$out_dir/ReadsBasedAnalysis/readsMappingToRef/$con_amb_fasta");
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

sub consensus_composition_info{
        my $comp_files = shift;
        my %comp;
        foreach my $file(@$comp_files){
                my ($file_name, $file_path, $file_suffix)=fileparse("$file", qr/\.[^.]*/);
                open (my $fh, "<", $file);
                my $id;
                while(<$fh>){
                        chomp;
                        next if ! /\w/;
                        if (/##(\S+)/){
                                $id = $1;
                                $id =~ s/\w+_consensus_(\S+)/$1/;
                                $comp{$id}->{file} = $file_name;
                                next;
                        }
                        my ($nuc,$num) = split /\t/,$_;
                        if ($nuc =~/Total length/){
                                $comp{$id}->{len} = $num;
                        }else {
                                $comp{$id}->{$nuc} = $num;
                        }
                }
                close $fh;
        }
        return \%comp;
}

sub pull_sampleMetadata {
	my $out_dir=shift;
	my $vars = shift;
	my $metadata = "$out_dir/metadata_gisaid_ncbi.txt";
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
             			$vars->{SM_BIOPROJECT_ID} =$2 if ($1 eq "bioproject");
             			$vars->{SM_RELEASE_DATE} =$2 if ($1 eq "release");
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
	my $metadata = "$userDir/gisaid_ncbi_submission_profile.txt";
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
             			$vars->{NCBIID} =$2 if ($1 eq "ncbi_id");
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
