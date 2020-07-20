#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../../lib";
use Spreadsheet::ParseExcel;
use Spreadsheet::ParseExcel::SaveParser;
use Getopt::Long;
use File::Basename;
use File::Path;

my $um; #user management
my $out;
my $project_dir_names;
my $userDir;
my $usage = qq{
Usage: $0
	Required
		-out        	      output file
		-projects        project_dir_names, separated by comma
		-udir		 user dir
};

GetOptions(
		"um=s"        =>  \$um,
		"out=s"        =>  \$out,
		"projects=s"       =>  \$project_dir_names,
		"udir=s"       => \$userDir,
		"help|?"           => sub{print "$usage\n";exit;} 
	);

if (!$project_dir_names && !$out){ print "$usage\n";exit;}

mkpath(dirname($out));

## read template
my $parser   = Spreadsheet::ParseExcel::SaveParser->new();
my $template = $parser->Parse("$FindBin::Bin/20200501_EpiCoV_BulkUpload_Template.xls");

if ( !defined $template ) {
    die $parser->error(), ".\n";
}
my $in_worksheet1 = $template->worksheet('Submissions');
my $in_worksheet0 = $template->worksheet('Instructions');
my $row=2;
my $col=0;


#write metadata to sheets
foreach my $proj_dir (split /,/,$project_dir_names){
	my $vars={};

	my $confFile = "$proj_dir/config.txt";
	my $metadataFile = "$proj_dir/metadata_sample.txt";
	my $conf = &getParams($confFile);
	my $metadata = &getParams($metadataFile);
	my $travelsFile = "$proj_dir/metadata_travels.txt";
	my $symptomsFile = "$proj_dir/metadata_symptoms.txt";
	my $otherFile = "$proj_dir/metadata_other.txt";

	my $proj_name = $conf->{'projname'};
	my $owner = $conf->{'projowner'};
	eval {  
		&pull_sampleMetadata($proj_dir,$vars);
		&pull_submissionData($proj_dir,$vars);
		&pull_consensusInfo($proj_dir,$vars,$conf);
	};
	$in_worksheet1->AddCell( $row, $col, $vars->{ID});  ## Submitter (login account id) *
	$in_worksheet1->AddCell( $row, $col+1, "all_sequences.fasta");  ## FASTA filename *
	$in_worksheet1->AddCell( $row, $col+2, $vars->{VIR_NAME});  ## Virus Name *
	$in_worksheet1->AddCell( $row, $col+3, "betacoronavirus");  ## Type *
	$in_worksheet1->AddCell( $row, $col+4, $vars->{VIR_PASSAGE}); ## Passage details/history *
	$in_worksheet1->AddCell( $row, $col+5, $vars->{SM_CDATE}); ## Collection date *
	$in_worksheet1->AddCell( $row, $col+6, $vars->{SM_LOC}); ## Location *
	$in_worksheet1->AddCell( $row, $col+7, ""); ## Additional location information
	$in_worksheet1->AddCell( $row, $col+8, $vars->{SM_HOST});  ## Host *
	$in_worksheet1->AddCell( $row, $col+9, "");  ## Additional host information
	$in_worksheet1->AddCell( $row, $col+10, $vars->{SM_GENDER});  ## Gender *
	$in_worksheet1->AddCell( $row, $col+11, $vars->{SM_AGE});  ## Patient Age *
	$in_worksheet1->AddCell( $row, $col+12, $vars->{SM_STATUS});  ## Patient status *
	$in_worksheet1->AddCell( $row, $col+13, "");  ## Specimen source
	$in_worksheet1->AddCell( $row, $col+14, "");  ## Outbreak
	$in_worksheet1->AddCell( $row, $col+15, "");  ## Last vaccinated
	$in_worksheet1->AddCell( $row, $col+16, "");  ## Treatment
	$in_worksheet1->AddCell( $row, $col+17, $vars->{SM_SEQUENCING_TECH});  ## Sequencing technology *
	$in_worksheet1->AddCell( $row, $col+18, $vars->{ASM_METHOD});  ## Assembly method 
	$in_worksheet1->AddCell( $row, $col+19, $vars->{SM_COV});  ## Coverage
	$in_worksheet1->AddCell( $row, $col+20, $vars->{ORIG_LAB});  ## Originating lab *
	$in_worksheet1->AddCell( $row, $col+21, $vars->{ORIG_ADDRESS});  ## Address *
	$in_worksheet1->AddCell( $row, $col+22, "");  ## Sample ID given by the sample provider
	$in_worksheet1->AddCell( $row, $col+23, $vars->{SUB_LAB});  ## Submitting lab *
	$in_worksheet1->AddCell( $row, $col+24, $vars->{SUB_ADDRESS});  ## Address *
	$in_worksheet1->AddCell( $row, $col+25, "");  ## Sample ID given by the submitting laboratory
	$in_worksheet1->AddCell( $row, $col+26, "");  ## Authors
	$in_worksheet1->AddCell( $row, $col+27, "");  ## Comment
	$in_worksheet1->AddCell( $row, $col+28, "");  ## Comment Icon
	$row++;
}

$template->SaveAs($out);

sub pull_consensusInfo{
	my $out_dir=shift;
	my $vars=shift;
	my $conf=shift;
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
                        my $selected = ( $vars->{SM_COV} && $vars->{SM_COV} =~ /$id/ )? 'selected':'';
                        $vars->{CON_LIST} .= "<option value='$value' $selected>$id ($linear_cov, $depth_cov)</option>";
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

        my $con_min_mapQ = $conf->{r2g_consensus_min_mapQ};
        my $con_min_cov = $conf->{r2g_consensus_min_cov};
        my $con_alt_prop = $conf->{r2g_consensus_alt_prop}*100 . "%";
        my $con_alt_indel = $conf->{r2g_consensus_altIndel_prop}* 100 . "%";

        $vars->{ASM_METHOD} .= " Consensus min coverage: ${con_min_cov}X. min map quality: $con_min_mapQ. Alternate Base > $con_alt_prop. Indel > $con_alt_indel.";
	my $platform = $conf->{FASTQ_SOURCE} || "";
        $vars->{FASTQ_SOURCE}=($platform eq "nanopore")?"Nanopore":"Iumina";
}

sub pull_sampleMetadata {
	my $out_dir=shift;
	my $vars=shift;
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
                close CONF;
        }
}
sub pull_submissionData {
	my $out_dir=shift;
	my $vars=shift;
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
sub getParams {
        my $config = shift;
        my $sys;

	if(-e $config) {
		open CONF, $config or die "Can't open $config: $!";
		while(<CONF>){
	      		chomp;
		        next if(/^#/);
		   	if ( /(.*)=(.*)/ ){
		     		$sys->{$1}=$2;
		      	}
		}
		close CONF;
	}
        return $sys;
}

1;
