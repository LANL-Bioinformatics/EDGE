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
my $configuration = &pull_EDGEConfig();

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
	$vars->{ASM_METHOD}="$align_tool $tool_version.";

	my $con_min_mapQ = $configuration->{r2g_consensus_min_mapQ};
	my $con_min_cov = $configuration->{r2g_consensus_min_cov};
	my $con_alt_prop = $configuration->{r2g_consensus_alt_prop}*100 . "%";
	my $con_alt_indel = $configuration->{r2g_consensus_altIndel_prop}* 100 . "%";
	
	$vars->{ASM_METHOD} .= " Consensus min coverage: ${con_min_cov}X. min map quality: $con_min_mapQ. Alternate Base > $con_alt_prop. Indel > $con_alt_indel.";
	open (my $config_fh, "<" , "$out_dir/config.txt");
	while(<$config_fh>){
		if (/fastq_source=(\S+)/i){ my $platform=$1; $vars->{FASTQ_SOURCE}=($platform eq "nanopore")?"Nanopore":"Illumina";}
	}
	close $config_fh
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
             			$vars->{SM_SEQUENCING_TECH} =$2 if ($1 eq "sequencing_technology");
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
sub pull_EDGEConfig
{
    my $file="$out_dir/config.txt";
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
