#!/usr/bin/perl
use strict;
#use warnings;
use FindBin qw($RealBin);
use File::Basename;
use lib "$RealBin/../../lib";
use HTML::Template;
use POSIX qw{strftime};
use JSON;

my $out_dir       = $ARGV[0];
my $html_outfile  = $ARGV[1];
my $settings_file          = $ARGV[2];

my $sysconfig    = "$RealBin/../../edge_ui/sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $settings = &getSettings($settings_file);

my $www_root	= $sys->{edgeui_wwwroot};
(my $report_rel_dir = $out_dir) =~ s/$www_root//;

## Instantiate the variables
my $vars;
$vars->{REPORTNAME} = $settings->{'report-name'};
$vars->{REPORTDESC} = $settings->{'report-desc'};
my $project_out_dir = $sys->{'edgeui_output'};
(my $project_rel_dir = $project_out_dir) =~ s/$www_root//;
my @projects = split(",",$settings->{'projects-selected'});


###project/run info
my %reports_map=();
$vars->{PROJURL} = $project_rel_dir;
$vars->{RUNINFO} = 0;
if($settings->{'run-name'} eq "on") {
	$vars->{RUNINFO} = 1;
	$vars->{RUNNAME} =1;
}
if($settings->{'run-desc'} eq "on") {
	$vars->{RUNINFO} = 1;
	$vars->{RUNDESC} =1;
}
if($settings->{'run-files'} eq "on") {
	$vars->{RUNINFO} = 1;
	$vars->{RUNFILES} =1;
}
if($settings->{'sample-metadata'} eq "on") {
	$vars->{RUNINFO} = 1;
	$vars->{SAMPLEMETADATA} =1;
	my @metadata_fields = split(",",$settings->{'sample-metadata-selected'});
	foreach my $field (@metadata_fields) {
		if($field eq "sample-name") {
			$vars->{SAMPLEMETADATA_NAME} =1;
		} elsif($field eq "sample-type") {
			$vars->{SAMPLEMETADATA_TYPE} =1;
		} elsif($field eq "host") {
			$vars->{SAMPLEMETADATA_HOST} =1;
		} elsif($field eq "isolation-source") {
			$vars->{SAMPLEMETADATA_SOURCE} =1;
		} elsif($field eq "collection-date") {
			$vars->{SAMPLEMETADATA_COLLECTIONDATE} =1;
		} elsif($field eq "location") {
			$vars->{SAMPLEMETADATA_LOCATION} =1;
		} elsif($field eq "sequencing-center") {
			$vars->{SAMPLEMETADATA_SEQCENTER} =1;
		} elsif($field eq "sequencer") {
			$vars->{SAMPLEMETADATA_SEQUENCER} =1;
		} elsif($field eq "sequencing-date") {
			$vars->{SAMPLEMETADATA_SEQDATE} =1;
		}
	}
}

###preprocess
$vars->{PREPROCESS} = 0;
if($settings->{'preprocess-stats'} eq "on") {
	my @pre_stats_fields = split(",",$settings->{'preprocess-stats-selected'});
	foreach my $field (@pre_stats_fields) {
		if($field eq "raw-reads") {
			$vars->{PREPROCESS} = 1;
			$vars->{PRE_RAWREADS} =1;
		} elsif($field eq "quality-trimming") {
			$vars->{PREPROCESS} = 1;
			$vars->{PRE_QUALTRIM} =1;
			$vars->{PRE_QUALTRIM_S} =1;
		} elsif($field eq "host-removal-filter") {
			$vars->{PREPROCESS} = 1;
			$vars->{PRE_HOSTRMF} =1;
			$vars->{PRE_HOSTRMF_S} =1;
		} 
	}
}
if($settings->{'preprocess-figures'} eq "on") {
	my @pre_figs_fields = split(",",$settings->{'preprocess-figures-selected'});
	foreach my $field (@pre_figs_fields) {
		if($field eq "quality-trimming") {
			$vars->{PREPROCESS} = 1;
			$vars->{PRE_QUALTRIM} =1;
			$vars->{PRE_QUALTRIM_F} =1;
		} elsif($field eq "host-removal-filter") {
			$vars->{PREPROCESS} = 1;
			$vars->{PRE_HOSTRMF} =1;
			$vars->{PRE_HOSTRMF_F} =1;
		} 
	}
}

####
$vars->{ASSEMBLYANNOTATION} = 0;
if($settings->{'assembly-annotation-stats'} eq "on") {
	my @stats_fields = split(",",$settings->{'assembly-annotation-stats-selected'});
	foreach my $field (@stats_fields) {
		if($field eq "assembly") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_ASSEMBLY} =1;
			$vars->{AA_ASSEMBLY_S} =1;
		} elsif($field eq "read-mapping") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_MAPPING} =1;
			$vars->{AA_MAPPING_S} =1;
		} elsif($field eq "annotation") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_ANNOTATION} =1;
			$vars->{AA_ANNOTATION_S} =1;
		} 
	}
}
if($settings->{'assembly-annotation-figures'} eq "on") {
	my @figs_fields = split(",",$settings->{'assembly-annotation-figures-selected'});
	foreach my $field (@figs_fields) {
		if($field eq "assembly") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_ASSEMBLY} =1;
			$vars->{AA_ASSEMBLY_F} =1;
		} elsif($field eq "read-mapping") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_MAPPING} =1;
			$vars->{AA_MAPPING_F} =1;
		} elsif($field eq "annotation") {
			$vars->{ASSEMBLYANNOTATION} = 1;
			$vars->{AA_ANNOTATION} =1;
			$vars->{AA_ANNOTATION_F} =1;
		} 
	}
}
###reference based analysis
$vars->{REF} = 0;
if($settings->{'ref-stats'} eq "on") {
	my @pre_stats_fields = split(",",$settings->{'ref-stats-selected'});
	foreach my $field (@pre_stats_fields) {
		if($field eq "mapped-reads") {
			$vars->{REF} = 1;
			$vars->{REF_MAPPED_READS} =1;
			$vars->{REF_MAPPED_READS_S} =1;
		} elsif($field eq "mapped-contigs") {
			$vars->{REF} = 1;
			$vars->{REF_MAPPED_CONTIGS} =1;
			$vars->{REF_MAPPED_CONTIGS_S} =1;
		} 
	}
}
if($settings->{'ref-figures'} eq "on") {
	my @pre_figs_fields = split(",",$settings->{'ref-figures-selected'});
	foreach my $field (@pre_figs_fields) {
		if($field eq "mapped-reads") {
			$vars->{REF} = 1;
			$vars->{REF_MAPPED_READS} =1;
			$vars->{REF_MAPPED_READS_F} =1;
		} elsif($field eq "mapped-contigs") {
			$vars->{REF} = 1;
			$vars->{REF_MAPPED_CONTIGS} =1;
			$vars->{REF_MAPPED_CONTIGS_F} =1;
		} 
	}
}
###Taxonomy Classification
$vars->{TAX} = 0;
if($settings->{'tax-stats'} eq "on") {
	my @stats_fields = split(",",$settings->{'tax-stats-selected'});
	foreach my $field (@stats_fields) {
		if($field eq "read") {
			$vars->{TAX} = 1;
			$vars->{TAX_READ} =1;
			$vars->{TAX_READ_S} =1;
		} elsif($field eq "assembly") {
			$vars->{TAX} = 1;
			$vars->{TAX_ASSEMBLY} =1;
			$vars->{TAX_ASSEMBLY_S} =1;
		} 
	}
}
if($settings->{'tax-figures'} eq "on") {
	my @figs_fields = split(",",$settings->{'tax-figures-selected'});
	foreach my $field (@figs_fields) {
		if($field eq "read") {
			$vars->{TAX} = 1;
			$vars->{TAX_READ} =1;
			$vars->{TAX_READ_F} =1;
		} elsif($field eq "assembly") {
			$vars->{TAX} = 1;
			$vars->{TAX_ASSEMBLY} =1;
			$vars->{TAX_ASSEMBLY_F} =1;
		} 
	}
}
my %tax_tools=();
my $pangia_score;
if($settings->{'tax-tools'} eq "on") {
	my @figs_fields = split(",",$settings->{'tax-tools-selected'});
	$pangia_score = $settings->{'tax-tool-pangia-score'};
	foreach my $field (@figs_fields) {
		if($field eq "gottcha-b-species") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_GBS} =1;
			$tax_tools{"gottcha-speDB-b"} = 1;
		} elsif($field eq "gottcha-v-species") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_GVS} =1;
			$tax_tools{"gottcha-speDB-v"} = 1;
		} elsif($field eq "gottcha2-b-species") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_GBS2} =1;
			$tax_tools{"gottcha2-speDB-b"} = 1;
		}  elsif($field eq "gottcha2-v-species") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_GVS2} =1;
			$tax_tools{"gottcha2-speDB-v"} = 1;
		}  elsif($field eq "pangia") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_PANGIA} =1;
			$tax_tools{"pangia"} = 1;
		}  elsif($field eq "metaphlan") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_METAPHLAN} =1;
			$tax_tools{"metaphlan"} = 1;
		}  elsif($field eq "metaphlan2") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_METAPHLAN} =1;
			$tax_tools{"metaphlan2"} = 1;
		}  elsif($field eq "bwa") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL_BWA} =1;
			$tax_tools{"bwa"} = 1;
		}  elsif($field eq "kraken-mini") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_KRAKEN} =1;
			$tax_tools{"kraken_mini"} = 1;
		}  elsif($field eq "kraken") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_KRAKEN} =1;
			$tax_tools{"kraken"} = 1;
		}  elsif($field eq "diamond") {
			$vars->{TAX} = 1;
			$vars->{TAX_TOOL} = 1;
			$vars->{TAX_TOOL_DIAMOND} =1;
			$tax_tools{"diamond"} = 1;
		} 
	}
}

my $out_file_name;
eval {
	&pull_reports();
	if($vars->{RUNINFO}) {
		$out_file_name = "runs_info.csv";
		&create_run_info_csv($out_file_name);
		$vars->{RUNINFOCSV} = "$report_rel_dir/$out_file_name";
	}
	if($vars->{PREPROCESS}) {
		if($vars->{PRE_RAWREADS}) {
			$out_file_name = "preprocess_raw_reads_stats.csv";
			&create_raw_reads_csv($out_file_name);
			$vars->{RAWREADSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{PRE_QUALTRIM_S}) {
			$out_file_name = "preprocess_trimmed_reads_stats.csv";
			&create_trimmed_reads_csv($out_file_name);
			$vars->{TRIMREADSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{PRE_HOSTRMF_S}) {
			$out_file_name = "preprocess_host_removal_stats.csv";
			&create_hostrmf_csv($out_file_name);
			$vars->{HOSTRMFREADSCSV} = "$report_rel_dir/$out_file_name";
		}
	}
	if($vars->{ASSEMBLYANNOTATION}) {
		if($vars->{AA_ASSEMBLY_S}) {
			$out_file_name = "assembly.csv";
			&create_assembly_csv($out_file_name);
			$vars->{ASSEMBLYCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{AA_MAPPING_S}) {
			$out_file_name = "read_mapping.csv";
			&create_read_mapping_csv($out_file_name);
			$vars->{READMAPPINGCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{AA_ANNOTATION_S}) {
			$out_file_name = "annotation.csv";
			&create_annotation_csv($out_file_name);
			$vars->{ANNOTATIONCSV} = "$report_rel_dir/$out_file_name";
		}
	}
	if($vars->{REF}) {
		if($vars->{REF_MAPPED_READS_S}) {
			$out_file_name = "ref_reads_stats.csv";
			&create_ref_reads_stats_csv($out_file_name);
			$vars->{REFREADSSTATSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "ref_reads_refs.csv";
			&create_ref_reads_refs_csv($out_file_name);
			$vars->{REFREADSREFSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{REF_MAPPED_CONTIGS_S}) {
			$out_file_name = "ref_contigs_stats.csv";
			&create_ref_contigs_stats_csv($out_file_name);
			$vars->{REFCONTIGSSTATSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "ref_contigs_refs.csv";
			&create_ref_contigs_refs_csv($out_file_name);
			$vars->{REFCONTIGSREFSCSV} = "$report_rel_dir/$out_file_name";
		}
	}
	if($vars->{TAX}) {
		if($vars->{TAX_READ_S}) {
			$out_file_name = "tax_read_based_genus_ranks.csv";
			&create_tax_read_genus_ranks_csv($out_file_name);
			$vars->{TAXREADGENUSRANKSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "tax_read_based_species_ranks.csv";
			&create_tax_read_species_ranks_csv($out_file_name);
			$vars->{TAXREADSPECIESRANKSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "tax_read_based_strain_ranks.csv";
			&create_tax_read_strain_ranks_csv($out_file_name);
			$vars->{TAXREADSTRAINRANKSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_ASSEMBLY_S}) {
			$out_file_name = "tax_assembly_based_contigs_stats.csv";
			&create_tax_assembly_contigs_stats_csv($out_file_name);
			$vars->{TAXASSEMBLYCONTIGSSTATSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "tax_assembly_based_hit_count_ranks.csv";
			&create_tax_assembly_count_ranks_csv($out_file_name);
			$vars->{TAXASSEMBLYHITRANKSCSV} = "$report_rel_dir/$out_file_name";
			$out_file_name = "tax_assembly_based_assigned_length_ranks.csv";
			&create_tax_assembly_len_ranks_csv($out_file_name);
			$vars->{TAXASSEMBLYLENRANKSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_GBS}) {
			$out_file_name = "tax_tool_gottcha_b_species.csv";
			&create_tax_tool_gottcha_csv($out_file_name, "tax-tool-gbs");
			$vars->{TAXTOOLGBSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_GVS}) {
			$out_file_name = "tax_tool_gottcha_v_species.csv";
			&create_tax_tool_gottcha_csv($out_file_name, "tax-tool-gvs");
			$vars->{TAXTOOLGVSCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_GBS2}) {
			$out_file_name = "tax_tool_gottcha2_b_species.csv";
			&create_tax_tool_gottcha2_csv($out_file_name, "tax-tool-gbs2");
			$vars->{TAXTOOLGBS2CSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_GVS2}) {
			$out_file_name = "tax_tool_gottcha2_v_species.csv";
			&create_tax_tool_gottcha2_csv($out_file_name, "tax-tool-gvs2");
			$vars->{TAXTOOLGVS2CSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_PANGIA}) {
			$out_file_name = "tax_tool_pangia.csv";
			&create_tax_tool_pangia_csv($out_file_name, "tax-tool-pangia");
			$vars->{TAXTOOLPANGIACSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_METAPHLAN}) {
			$out_file_name = "tax_tool_metaphlan.csv";
			&create_tax_tool_other_csv($out_file_name, "tax-tool-metaphlan");
			$vars->{TAXTOOLMETAPHLANCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_BWA}) {
			$out_file_name = "tax_tool_bwan.csv";
			&create_tax_tool_other_csv($out_file_name, "tax-tool-bwa");
			$vars->{TAXTOOLBWACSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_KRAKEN}) {
			$out_file_name = "tax_tool_kraken.csv";
			&create_tax_tool_other_csv($out_file_name, "tax-tool-kraken");
			$vars->{TAXTOOLKRAKENCSV} = "$report_rel_dir/$out_file_name";
		}
		if($vars->{TAX_TOOL_DIAMOND}) {
			$out_file_name = "tax_tool_diamond.csv";
			&create_tax_tool_other_csv($out_file_name, "tax-tool-diamond");
			$vars->{TAXTOOLDIAMONDCSV} = "$report_rel_dir/$out_file_name";
		}
	}			
};

#print STDOUT "start output html \n";
output_html();
#print STDOUT "end output html\n";
#
my $now_string = strftime "%Y %b %e %H:%M:%S", localtime;
`echo 'Created Time=$now_string' >> $settings_file`;
#

sub output_html {
	#reformat number with thousand separator
	foreach my $var ( keys %$vars ){
		next if ($var eq "REPORTID");
		if( ref($vars->{$var}) eq 'ARRAY' ){
			for( my $i=0; $i < scalar(@{$vars->{$var}}); $i++ ){
				if( ref($vars->{$var}[$i]) eq 'HASH' ){
					foreach my $idx ( keys %{$vars->{$var}[$i]} ){
						$vars->{$var}[$i]->{$idx} = &_reformat_val($vars->{$var}[$i]->{$idx});
					}
				}
				else{
					$vars->{$var}[$i] = &_reformat_val($vars->{$var}[$i]);
				}
			}
		}
		else{
			$vars->{$var} = &_reformat_val($vars->{$var});
		}
	}

	my $template = HTML::Template->new(filename => "$RealBin/edge_projects_report_html.tmpl",strict => 0,die_on_bad_params => 0);
	$template->param(%$vars);

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
		
	system("cp -r $www_root/css $out_dir"); 
	system("cp -r $www_root/images $out_dir");  
	system("cp -r $www_root/javascript $out_dir");  

}

sub getSettings {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		chomp;
		if ( /^([^=]+)=([^=]+)/ ){
			$sys->{$1}=$2;
		}
	}
	close CONF;
	return $sys;
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

sub _reformat_val {
	my $content = shift;
	return unless defined $content;

	my $outdir_rgx = $out_dir;
	$outdir_rgx =~ s/\//\\\//g;

	if( $content =~ /^\d{4,}$/ ){
		$content =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	}
	elsif( $content =~ /^(\d{4,}) \(([\d\.]+ \%)\)$/ ){
		my ($num,$pct) = ($1,$2);
		$num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
		$content = "$num ($pct)";
	}
	elsif( $content =~ /^(\d+\.\d{2})\d+$/ ){
		$content = $1;
	}
	elsif( $content =~ /$outdir_rgx/ ){
		#$content =~ s/$outdir_rgx/\.\./ unless $mode eq "web";
	}
	return $content;
}

#parse selected projects' reports
sub pull_reports {
	foreach my $project (@projects) {
		my $proj;
		my $out_dir = "$project_out_dir/$project";
		$proj->{PROJ} = $project;
		my $report = "$out_dir/HTML_Report/report_web.html";
		my $pangiaReport = "$out_dir/ReadsBasedAnalysis/Taxonomy/report/1_allReads/pangia/allReads-pangia.list.txt";
		if(!-e $report) {
			##call outputMunger_w_temp.pl to generate report
			`$RealBin/../munger/outputMunger_w_temp.pl $out_dir $report`;
		}
		
		### get runhost
		open CONF, "$out_dir/config.txt" or die $!;
		while(<CONF>) {
			chomp;
     			if (/^projrunhost=(.*)/) {
				$proj->{RUNHOST} = $1;
				last;
			}
		}
		close CONF;
		###
		#print "open $report\n";
		open IN, "$report" or die $!;
		while(<IN>) {
			chomp;
			s/^\s+|\s+$//g;
			if(/<h2 id="edge-output-projname" data-pid="(.*)">(.*)<\/h2>/) {
				$proj->{PROJNAME} = $2;
				$reports_map{$project}{'run-name'} = $2;
				next;
			} 
			if(/Description: (.*)<br\/>/) {
				$proj->{PROJDESC} = $1;
				$reports_map{$project}{'run-desc'} = $1;
				next;
			} 
			if(/<th data-priority='1'>Report\/Info<\/th>/) {
				#get input files
				my $get = 0;
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					if(/<td>Input Reads<\/td>/) {
						$get = 1;
						next;
					}
					if($get && /<td>(.*)<\/td>/) {
						my $files = $1;
						$files =~ s/^[,\s+]//g;
						$files =~ s/,\s*$//g;
						$files =~ s/,\s*$//g;
						$proj->{PROJFILES} = $files;
						$reports_map{$project}{'run-files'} = $files;
						last;
					}
				}
				next;
			} 
			### sample metadata
			if(/<h4>Sample Metadata<\/h4>/) {
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					if(/<tr><td>Sample Name<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SAMPLENAME} = $1;
						$reports_map{$project}{'sample-name'} = $1;
						next;
					}
					if(/<tr><td>Sample Type<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SAMPLETYPE} = $1;
						$reports_map{$project}{'sample-type'} = $1;
						next;
					}
					if(/<tr><td>Host<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SAMPLEHOST} = $1;
						$reports_map{$project}{'sample-host'} = $1;
						next;
					}
					if(/<tr><td>Isolation Source<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{ISOLATIONSOURCE} = $1;
						$reports_map{$project}{'isolation-source'} = $1;
						next;
					}
					if(/<tr><td>Collection Date<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{COLLECTIONDATE} = $1;
						$reports_map{$project}{'collection-date'} = $1;
						next;
					}
					if(/<tr><td>Location<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SAMPLELOCATION} = $1;
						$reports_map{$project}{'sample-location'} = $1;
						next;
					}
					if(/<tr><td>Sequencing Center<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SEQCENTER} = $1;
						$reports_map{$project}{'sequencing-center'} = $1;
						next;
					}
					if(/<tr><td>Sequencer<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SEQUENCER} = $1;
						$reports_map{$project}{'sequencer'} = $1;
						next;
					}
					if(/<tr><td>Sequencing Date<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{SEQDATE} = $1;
						$reports_map{$project}{'sequencing-date'} = $1;
						last;
					}
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}
				}
				next;
			}## end sample metadata

			## preprocess data
			if(/<h4>Pre-processing<\/h4>/) {
				my $type;
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					if(/<th data-priority='1'>Raw Reads<\/th>/) {
						$type = "raw-reads";
						next;
					}
					if(/<th data-priority='1'>Trimmed Reads<\/th>/) {
						$type = "trimmed-reads";
						next;
					}
					if(/<tr><td>Reads<\/td><td>(.*)<\/td><\/tr>/) {
						if($type eq "raw-reads") {
							$proj->{PRERAWREADS} = $1;
							$reports_map{$project}{'preprocess-raw-reads'} = $1;
						} elsif($type eq "trimmed-reads") {
							if($1 =~ /(.*) \((.*)\)/) {
								$proj->{PRETRIMREADS} = $1;
								$reports_map{$project}{'preprocess-trimmed-reads'} = $1;
								$proj->{PRETRIMREADSPCT} = $2;
								$reports_map{$project}{'preprocess-trimmed-reads-percent'} = $2;
							}
						} 
						next;
					}
					if(/<tr><td>Total Bases<\/td><td>(.*)<\/td><\/tr>/) {
						if($type eq "raw-reads") {
							$proj->{PRERAWBASES} = $1;
							$reports_map{$project}{'preprocess-raw-bases'} = $1;
						} elsif($type eq "trimmed-reads") {
							if($1 =~ /(.*) \((.*)\)/) {
								$proj->{PRETRIMBASES} = $1;
								$reports_map{$project}{'preprocess-trimmed-bases'} = $1;
								$proj->{PRETRIMBASESPCT} = $2;
								$reports_map{$project}{'preprocess-trimmed-bases-percent'} = $2;
							}
						} 
						next;
					}
					if(/<tr><td>Mean Read Length<\/td><td>(.*)<\/td><\/tr>/) {
						if($type eq "raw-reads") {
							$proj->{PRERAWREADLEN} = $1;
							$reports_map{$project}{'preprocess-raw-read-length'} = $1;
						} elsif($type eq "trimmed-reads") {
							$proj->{PRETRIMREADLEN} = $1;
							$reports_map{$project}{'preprocess-trimmed-read-length'} = $1;
						} 
						next;
					}
					if(/<tr><td>Paired Reads<\/td><td>(.*)<\/td><\/tr>/) {
						if($1 =~ /(.*) \((.*)\)/) {
							$proj->{PRETRIMPAIREDREADS} = $1;
							$reports_map{$project}{'preprocess-trimmed-paired-reads'} = $1;
							$proj->{PRETRIMPAIREDREADSPCT} = $2;
							$reports_map{$project}{'preprocess-trimmed-paired-reads-percent'} = $2;
						}
						next;
					}
					if(/<tr><td>Paired Total Bases<\/td><td>(.*)<\/td><\/tr>/) {
						if($1 =~ /(.*) \((.*)\)/) {
							$proj->{PRETRIMPAIREDBASES} = $1;
							$reports_map{$project}{'preprocess-trimmed-paired-bases'} = $1;
							$proj->{PRETRIMPAIREDBASESPCT} = $2;
							$reports_map{$project}{'preprocess-trimmed-paired-bases-percent'} = $2;
						}
						next;
					}
					if(/<tr><td>Unpaired Reads<\/td><td>(.*)<\/td><\/tr>/) {
						if($1 =~ /(.*) \((.*)\)/) {
							$proj->{PRETRIMUNPAIREDREADS} = $1;
							$reports_map{$project}{'preprocess-trimmed-unpaired-reads'} = $1;
							$proj->{PRETRIMUNPAIREDREADSPCT} = $2;
							$reports_map{$project}{'preprocess-trimmed-unpaired-reads-percent'} = $2;
						}
						next;
					}
					if(/<tr><td>Unpaired Total Bases<\/td><td>(.*)<\/td><\/tr>/) {
						if($1 =~ /(.*) \((.*)\)/) {
							$proj->{PRETRIMUNPAIREDBASES} = $1;
							$reports_map{$project}{'preprocess-trimmed-unpaired-bases'} = $1;
							$proj->{PRETRIMUNPAIREDBASESPCT} = $2;
							$reports_map{$project}{'preprocess-trimmed-unpaired-bases-percent'} = $2;
						}
						next;
					}
					if(/<img class="preview_img" data-src="(.*QC_quality_report\.png)" alt="QC_quality_report">/) {
						$proj->{PREQCQUALPNG} = $1;
						$proj->{PREQCQUALPNG_TXT} = "PNG";
						$reports_map{$project}{'preprocess-qc-qual-png'} = $1;
						next;
					}
					if(/<img class="preview_img" data-src="(.*QC_read_length\.png)" alt="QC_quality_report">/) {
						$proj->{PREQCREADLENPNG} = $1;
						$proj->{PREQCREADLENPNG_TXT} = "PNG";
						$reports_map{$project}{'preprocess-qc-read-length-png'} = $1;
						next;
					}
					if(/<img class="preview_img" data-src="(.*QC_nucleotide_content\.png)" alt="QC_quality_report">/) {
						$proj->{PREQCNUCPNG} = $1;
						$proj->{PREQCNUCPNG_TXT} = "PNG";
						$reports_map{$project}{'preprocess-nucleotide-png'} = $1;
						next;
					}
					if(/<img class="preview_img" data-src="(.*QC_quality_boxplot\.png)" alt="QC_quality_report">/) {
						$proj->{PREQCPLOTPNG} = $1;
						$proj->{PREQCPLOTPNG_TXT} = "PNG";
						$reports_map{$project}{'preprocess-qc-boxplot-png'} = $1;
						next;
					}
					if(/<a data-ajax='false' href='(.*QC_qc_report\.pdf)'> QC Report PDF <\/a>/) {
						$proj->{PREQCREPORTPDF} = $1;
						$proj->{PREQCREPORTPDF_TXT} = "PDF";
						$reports_map{$project}{'preprocess-qc-report-pdf'} = $1;
						last;
					}					
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}
				}
				next;
			}
			if(/<li><span class="li-report-content-title">Host Removal and Filter<\/span><div class="li-report-content">/) {
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					if(/<p>Host file: (.*)\.<\/p>/) {
						$proj->{PREHOSTFILE} = $1;
						$reports_map{$project}{'preprocess-host-file'} = $1;
						next;
					}
					if(/<tr><td>Total reads<\/td><td>(.*)<\/td><\/tr>/) {
						$proj->{PREHOSTRMFREADS} = $1;
						$reports_map{$project}{'preprocess-host-removal-reads'} = $1;
						next;
					}
					if(/<tr><td>Total non-host reads<\/td><td>(.*) \((.*)\)<\/td><\/tr>/) {
						$proj->{PREHOSTRMFNONHOSTREADS} = $1;
						$reports_map{$project}{'preprocess-host-removal-nonhost-reads'} = $1;
						$proj->{PREHOSTRMFNONHOSTREADSPCT} = $2;
						$reports_map{$project}{'preprocess-host-removal-nonhost-reads-percent'} = $2;
						next;
					}
					if(/<img class="preview_img" data-src="(.*HostRemovalStats\.png)" alt="Host Removal Stats">/) {
						$proj->{PREHOSTRMFPNG} = $1;
						$proj->{PREHOSTRMFPNG_TXT} = "PNG";
						$reports_map{$project}{'preprocess-host-removal-stats-png'} = $1;
						next;
					}
					if(/<a data-ajax='false' href='(.*HostRemovalStats\.pdf)'> Report PDF <\/a>/) {
						$proj->{PREHOSTRMFREPORTPDF} = $1;
						$proj->{PREHOSTRMFREPORTPDF_TXT} = "PDF";
						$reports_map{$project}{'preprocess-host-removal-report-pdf'} = $1;
						last;
					}
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}

				}
				next;
			}
			## end preprocess data

			## assembly and annotation			
			if(/<h4>Assembly and Annotation<\/h4>/) {
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					##assembly
					if(/<span class="li-report-content-title">De Novo Assembly by (.*)<\/span><div class="li-report-content">/ || /<li><span class="li-report-content-title">Stats of Provided Contigs<\/span><div class="li-report-content">/) {
						$proj->{AA_ASSEMBLY_TOOL} = $1;
						$reports_map{$project}{'aa-assembly-tool'} = $1;
						while(<IN>) {
							if(/<tr><td>Number of contigs<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ASSEMBLY_CONTIGS} = $1;
								$reports_map{$project}{'aa-assembly-contigs'} = $1;
								next;
							}
							if(/<tr><td>N50<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ASSEMBLY_N50} = $1;
								$reports_map{$project}{'aa-assembly-n50'} = $1;
								next;
							}
							if(/<tr><td>Max contig size<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ASSEMBLY_MAXCONTIGSIZE} = $1;
								$reports_map{$project}{'aa-assembly-maxcontigsize'} = $1;
								next;
							}
							if(/<tr><td>Min contig size<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ASSEMBLY_MINCONTIGSIZE} = $1;
								$reports_map{$project}{'aa-assembly-mincontigsize'} = $1;
								next;
							}
							if(/<tr><td>total assembly size<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ASSEMBLY_SIZE} = $1;
								$reports_map{$project}{'aa-assembly-size'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*Assembly_length\.png)" alt="Assembly length">/) {
								$proj->{AA_ASSEMBLY_LENPNG} = $1;
								$proj->{AA_ASSEMBLY_LENPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-assembly-len-png'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*Assembly_GC_content\.png)" alt="Assembly_GC_content">/) {
								$proj->{AA_ASSEMBLY_GCPNG} = $1;
								$proj->{AA_ASSEMBLY_GCPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-assembly-gc-png'} = $1;
								next;
							}
							if(/<a data-ajax='false' href='(.*contigs_stats\.pdf)'> Report PDF <\/a>/) {
								$proj->{AA_ASSEMBLY_REPORTPDF} = $1;
								$proj->{AA_ASSEMBLY_REPORTPDF_TXT} = "PDF";
								$reports_map{$project}{'aa-assembly-report-pdf'} = $1;
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}
						}
						next;
					}
					##read mapping
					if(/<span class="li-report-content-title">Assembly Validation by Read Mapping<\/span><div class="li-report-content">/) {
						while(<IN>) {
							if(/<tr><td>Number of Mapped Reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_MAPPING_MAPPEDREADS} = $1;
								$reports_map{$project}{'aa-mapping-mappedreads'} = $1;
								next;
							}
							if(/<tr><td>% of Total Reads<\/td><td>(.*)<\/td><\/tr>/) {
								if($proj->{AA_MAPPING_MAPPEDPCT}) {
									$proj->{AA_MAPPING_UNMAPPEDPCT} = $1;
									$reports_map{$project}{'aa-mapping-unmappedpct'} = $1;
								} else {
									$proj->{AA_MAPPING_MAPPEDPCT} = $1;
									$reports_map{$project}{'aa-mapping-mappedpct'} = $1;
								}
								next;
							}
							if(/<tr><td>Number of Unmapped Reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_MAPPING_UNMAPPEDREADS} = $1;
								$reports_map{$project}{'aa-mapping-unmappedreads'} = $1;
								next;
							}
							if(/<tr><td>Average Fold Coverage<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_MAPPING_COV} = $1;
								$reports_map{$project}{'aa-mapping-cov'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*Assembly_CovDepth_vs_Len\.png)" alt="Assembly Depth vs Len">/) {
								$proj->{AA_MAPPING_COVDLPNG} = $1;
								$proj->{AA_MAPPING_COVDLPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-mapping-covdl-png'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*Assembly_Cov_vs_Len\.png)" alt="Assembly Coverage vs Len">/) {
								$proj->{AA_MAPPING_COVLPNG} = $1;
								$proj->{AA_MAPPING_COVLPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-mapping-covl-png'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*Assembly_GC_vs_CovDepth\.png)" alt="Assembly GC vs Depth">/) {
								$proj->{AA_MAPPING_COVGDPNG} = $1;
								$proj->{AA_MAPPING_COVGDPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-mapping-covgd-png'} = $1;
								next;
							}
							if(/<a data-ajax='false' href='(.*readsToContigs_plots\.pdf)'> Report PDF <\/a>/) {
								$proj->{AA_MAPPING_REPORTPDF} = $1;
								$proj->{AA_MAPPING_REPORTPDF_TXT} = "PDF";
								$reports_map{$project}{'aa-mapping-report-pdf'} = $1;
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}
						}
						next;
					}
					##annotation
					if(/<span class="li-report-content-title">Annotation<\/span><div class="li-report-content">/) {
						while(<IN>) {
							if(/<tr><td>CDS<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ANNOTATION_CDS} = $1;
								$reports_map{$project}{'aa-annotation-cds'} = $1;
								next;
							}
							if(/<tr><td>rRNA<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ANNOTATION_RRNA} = $1;
								$reports_map{$project}{'aa-annotation-rrna'} = $1;
								next;
							}
							if(/<tr><td>tRNA<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{AA_ANNOTATION_TRNA} = $1;
								$reports_map{$project}{'aa-annotation-trna'} = $1;
								next;
							}
							if(/<img class="preview_img" data-src="(.*annotation_stats_plots\.png)" alt="Annotation Stats">/) {
								$proj->{AA_ANNOTATION_STATSPNG} = $1;
								$proj->{AA_ANNOTATION_STATSPNG_TXT} = "PNG";
								$reports_map{$project}{'aa-annotation-stats-png'} = $1;
								next;
							}
							if(/<a data-ajax='false' href='(.*\.gff)'> GFF3 <\/a>/) {
								$proj->{AA_ANNOTATION_GFF3} = $1;
								$proj->{AA_ANNOTATION_GFF3_TXT} = "GFF3";
								$reports_map{$project}{'aa-annotation-gff3'} = $1;
								next;
							}
							if(/<a data-ajax='false' href='(.*\.gbk)'> GenBank <\/a>/) {
								$proj->{AA_ANNOTATION_GBK} = $1;
								$proj->{AA_ANNOTATION_GBK_TXT} = "GBK";
								$reports_map{$project}{'aa-annotation-gbk'} = $1;
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}
						}
					}
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}
				}
				next;
			}
			## end aa

			## Reference based analysis			
			if(/<h4>Reference-Based Analysis<\/h4>/) {
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					##reads mapped
					my $img;
					if(/<li><span class="li-report-content-title">Reads Mapped to Reference\(s\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/<tr><td>Number of Mapped Reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_MAPPED_READS} = $1;
								$reports_map{$project}{'ref-mapped-reads'} = $1;
								next;
							}
							if(/<tr><td>% of Total Post-QC Reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_POST_QC_READS_PCT} = $1;
								$reports_map{$project}{'ref-post-qc-reads-pct'} = $1;
								next;
							}
							if(/<tr><td>Average Fold<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_READ_AVG_FOLD} = $1;
								$reports_map{$project}{'ref-read-avg-fold'} = $1;
								next;
							}
							if(/<tr><td>Linear Coverage<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_READ_LINEAR_COV} = $1;
								$reports_map{$project}{'ref-read-linear-cov'} = $1;
								next;
							}
							if(/<tr><td>SNPs<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_READ_SNPS} = $1;
								$reports_map{$project}{'ref-read-snps'} = $1;
								next;
							}
							if(/<tr><td>InDels<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_READ_INDELS} = $1;
								$reports_map{$project}{'ref-read-indels'} = $1;
								next;
							}

							###get references
							if(/<th data-priority='1'>Reference<\/th>/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse ref
										my $ref;
										$ref->{PROJNAME} = $proj->{PROJNAME};
										$ref->{PROJ} = $proj->{PROJ};
										$ref->{PROPROJURL} = $vars->{PROJURL};
										if($str =~ /<td title=.*><a href='(.*)'>(.*)<\/a><\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											$ref->{LS_REF_LINK} = $1;
											$ref->{LS_REF} = $2;
											$ref->{LS_REF_NAME} = $3;
											$ref->{LS_REF_LEN} = $4;
											$ref->{LS_REF_CONSENSUS} = $5;
											$ref->{LS_REF_GC} = $6;
											$ref->{LS_REF_MAPPEDREADS} = $7;
											$ref->{LS_REF_MAPPEDREADS_PCT} = $8;
											$ref->{LS_REF_BASECOV} = $9;
											$ref->{LS_REF_AVGFOLD} = $10;
											$ref->{LS_REF_FOLDSTD} = $11;
											$ref->{LS_REF_GAPS} = $12;
											$ref->{LS_REF_GAPBASES} = $13;
											$ref->{LS_REF_SNPS} = $14;
											$ref->{LS_REF_INDELS} = $15;
											$ref->{LS_REF_MAPPEDREADS} =~ s/<a href.*>(.*)<\/a>/${1}/;
										} elsif($str =~ /<td title=.*><a href='(.*)'>(.*)<\/a><\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											$ref->{LS_REF_LINK} = $1;
											$ref->{LS_REF} = $2;
											$ref->{LS_REF_NAME} = $3;
											$ref->{LS_REF_LEN} = $4;
											$ref->{LS_REF_GC} = $5;
											$ref->{LS_REF_MAPPEDREADS} = $6;
											$ref->{LS_REF_MAPPEDREADS_PCT} = $7;
											$ref->{LS_REF_BASECOV} = $8;
											$ref->{LS_REF_AVGFOLD} = $9;
											$ref->{LS_REF_FOLDSTD} = $10;
											$ref->{LS_REF_GAPS} = $11;
											$ref->{LS_REF_GAPBASES} = $12;
											$ref->{LS_REF_SNPS} = $13;
											$ref->{LS_REF_INDELS} = $14;
											$ref->{LS_REF_MAPPEDREADS} =~ s/<a href.*>(.*)<\/a>/${1}/;
										} elsif($str =~ /<td title=.*><a href='(.*)'>(.*)<\/a><\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											$ref->{LS_REF_LINK} = $1;
											$ref->{LS_REF} = $2;
											$ref->{LS_REF_NAME} = $3;
											$ref->{LS_REF_LEN} = $4;
											$ref->{LS_REF_GC} = $5;
											$ref->{LS_REF_MAPPEDREADS} = $6;
											$ref->{LS_REF_BASECOV} = $7;
											$ref->{LS_REF_AVGFOLD} = $8;
											$ref->{LS_REF_FOLDSTD} = $9;
											$ref->{LS_REF_GAPS} = $10;
											$ref->{LS_REF_GAPBASES} = $11;
											$ref->{LS_REF_SNPS} = $12;
											$ref->{LS_REF_INDELS} = $13;
											$ref->{LS_REF_MAPPEDREADS} =~ s/<a href.*>(.*)<\/a>/${1}/;
										}
										push @{$proj->{REF_MAPPEDREADS_REFS_LOOP}}, $ref;
										push @{$reports_map{$project}{'ref-mappedreads-refs'}}, $ref;
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
								next;
							}
							
							### get images
							if(/<img class="preview_img" data-src="(.*\.png)" alt="Coverage">/) {
								undef $img;
								$img->{PROJNAME} = $proj->{PROJNAME};
								$img->{PROJ} = $proj->{PROJ};
								$img->{PROPROJURL} = $vars->{PROJURL};
								$img->{REF_COV_PNG} = $1;
								$img->{REF_COV_PNG_TXT} = "PNG";
								$proj->{REF_MAPPEDREADS_PNGS} = 1;
								if($img->{REF_COV_PNG} =~ /readsToRef_(.*)_base_coverage\.png/) {
									$img->{REF_PNG_REF} = $1;
								}
								next;
							}
							if(/<img class="preview_img" data-src="(.*\.png)" alt="Fold histogram">/) {
								$img->{REF_FOLDHIST_PNG} = $1;
								$img->{REF_FOLDHIST_PNG_TXT} = "PNG";
								push @{$proj->{REF_MAPPEDREADS_PNGS_LOOP}}, $img;
								$proj->{REF_MAPPEDREADS_PNGS} = 1;
								next;
							}
							if(/<a data-ajax='false' href='(.*readsToRef_plots\.pdf)'> All Plots PDF <\/a>/) {
								$proj->{REF_MAPPEDREADS_ALLPLOTS_PDF} = $1;
								$proj->{REF_MAPPEDREADS_ALLPLOTS_PDF_TXT} = "PDF";
								next;
							}
							if(/<a data-ajax='false' href='(.*readsToRef.SNPs_report\.txt)'> SNP Report <\/a>/) {
								$proj->{REF_MAPPEDREADS_SNPREPORT} = $1;
								$proj->{REF_MAPPEDREADS_SNPREPORT_TXT} = "TXT";
								next;
							}

							##unmapped reads
							if(/<tr><td>Number of Unmapped Reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_UNMAPPED_READS} = $1;
								$reports_map{$project}{'ref-unmapped-reads'} = $1;
								next;
							}
							if(/<tr><td>% of total reads<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_UNMAPPED_READS_PCT} = $1;
								$reports_map{$project}{'ref-unmapped-reads-pct'} = $1;
								next;
							}
							if(/<a data-ajax='false' href='(.*UnmappedReads_coverage\.txt)'> Mapping Result <\/a>/) {
								$proj->{REF_UNMAPPEDREADS_COV_FILE} = $1;
								$proj->{REF_UNMAPPEDREADS_COV_FILE_TXT} = "TXT";
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}
							
							if(/<li><span class="li-report-content-title">Contigs Mapped to Reference\(s\)<\/span><div class="li-report-content">/) {
								last;
							}
						}
					}
					##contigs mapped
					if(/<li><span class="li-report-content-title">Contigs Mapped to Reference\(s\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							if(/<tr><td>Number of Mapped Contigs<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_MAPPED_CONTIGS} = $1;
								$reports_map{$project}{'ref-mapped-contigs'} = $1;
								next;
							}
							if(/<tr><td>Proportion<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_PROPORTION_PCT} = $1;
								$reports_map{$project}{'ref-proportion-pct'} = $1;
								next;
							}
							if(/<tr><td>Average Fold<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_CONTIG_AVG_FOLD} = $1;
								$reports_map{$project}{'ref-contig-avg-fold'} = $1;
								next;
							}
							if(/<tr><td>Linear Coverage<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_CONTIG_LINEAR_COV} = $1;
								$reports_map{$project}{'ref-contig-linear-cov'} = $1;
								next;
							}
							if(/<tr><td>Average Identity<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_CONTIG_AVG_IDT_PCT} = $1;
								$reports_map{$project}{'ref-contig-avg-idt-pct'} = $1;
								next;
							}
							if(/<tr><td>SNPs<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_CONTIG_SNPS} = $1;
								$reports_map{$project}{'ref-contig-snps'} = $1;
								next;
							}
							if(/<tr><td>InDels<\/td><td>(.*)<\/td><\/tr>/) {
								$proj->{REF_CONTIG_INDELS} = $1;
								$reports_map{$project}{'ref-contig-indels'} = $1;
								next;
							}

							###get references
							if(/<th data-priority='1'>Reference<\/th>/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse ref
										if($str =~ /<td title=.*><a href='(.*)'>(.*)<\/a><\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $ref;
											$ref->{PROJNAME} = $proj->{PROJNAME};
											$ref->{PROJ} = $proj->{PROJ};
											$ref->{PROPROJURL} = $vars->{PROJURL};
											$ref->{LS_REF_LINK} = $1;
											$ref->{LS_REF} = $2;
											$ref->{LS_REF_NAME} = $3;
											$ref->{LS_REF_LEN} = $4;
											$ref->{LS_REF_GC} = $5;
											$ref->{LS_REF_MAPPEDCONTIGS} = $6;
											$ref->{LS_REF_MAPPEDCONTIGS_PCT} = $7;
											$ref->{LS_REF_BASECOV} = $8;
											$ref->{LS_REF_AVGFOLD} = $9;
											$ref->{LS_REF_GAPS} = $10;
											$ref->{LS_REF_GAPBASES} = $11;
											$ref->{LS_REF_SNPS} = $12;
											$ref->{LS_REF_INDELS} = $13;
											$ref->{LS_REF_MAPPEDCONTIGS} =~ s/<a href.*>(.*)<\/a>/${1}/;
											push @{$proj->{REF_MAPPEDCONTIGS_REFS_LOOP}}, $ref;
											push @{$reports_map{$project}{'ref-mappedcontigs-refs'}}, $ref;
										}
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
							}
							##figures and files
							if(/<a data-ajax='false' href='(.*contigsToRef_plot\.pdf)'> All Plots PDF<\/a>/) {
								$proj->{REF_MAPPEDCONTIGS_ALLPLOTS_PDF} = $1;
								$proj->{REF_MAPPEDCONTIGS_ALLPLOTS_PDF_TXT} = "PDF";
								next;
							}
							if(/<a data-ajax='false' href='(.*UnmappedContigs\.ctg_class\.top\.csv)'> Result <\/a>/) {
								$proj->{REF_MAPPEDCONTIGS_UNMAPPED_CSV} = $1;
								$proj->{REF_MAPPEDCONTIGS_UNMAPPED_CSV_TXT} = "CSV";
								last;
							}
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src=.*> Directory <\/a><\/p>/) {
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}							
						}
						next;
					}
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}
				}
				next;
			}
			## end ref

			## Taxonomy Classification			
			if(/<h4>Taxonomy Classification<\/h4>/) {
				while(<IN>) {
					chomp;
					s/^\s+|\s+$//g;
					##Read-based
					my $img;
					if(/<span class="li-report-content-title">Read-based Taxonomy Classification<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							###get tool results
							if(/<li><span class="li-report-content-title">Comparative Overview of Selected Tools<\/span><div class="li-report-content">/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse ref
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $ref;
											$ref->{PROJNAME} = $proj->{PROJNAME};
											$ref->{PROJ} = $proj->{PROJ};
											$ref->{PROPROJURL} = $vars->{PROJURL};
											$ref->{LS_TOOL} = $1;
											$ref->{LS_READS} = $2;
											$ref->{LS_READS_PCT} = $3;
											$ref->{LS_RANK} = $4;
											$ref->{LS_TOP1} = $5;
											$ref->{LS_TOP2} = $6;
											$ref->{LS_TOP3} = $7;
											$ref->{LS_TOP4} = $8;
											$ref->{LS_TOP5} = $9;
											my $add = 0;
											if($vars->{TAX_TOOL}) {
												if($tax_tools{$ref->{LS_TOOL}}){
													$add = 1;
												}
											} else {
												$add = 1;
											}
											if($add) {
												if($ref->{LS_RANK} eq 'genus') {
													push @{$proj->{TAX_READ_GENUS_LOOP}}, $ref;
													push @{$reports_map{$project}{'tax-read-genus'}}, $ref;
												}
												elsif($ref->{LS_RANK} eq 'species') {
													push @{$proj->{TAX_READ_SPECIES_LOOP}}, $ref;
													push @{$reports_map{$project}{'tax-read-species'}}, $ref;
												}
												elsif($ref->{LS_RANK} eq 'strain') {
													push @{$proj->{TAX_READ_STRAIN_LOOP}}, $ref;
													push @{$reports_map{$project}{'tax-read-strain'}}, $ref;
												}
											}
										}
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
							}
							
							### get images/files
							if(/<span class="iframe_label">Heatmap <a data-ajax='false' href='(.*heatmap_DATASET-allReads\.species\.pdf)'>[full]<\/a><\/span>/) {
								$proj->{TAX_READ_SPECIES_HEATMAP_PDF} = $1;
								$proj->{TAX_READ_SPECIES_HEATMAP_PDF_TXT} = "PDF";
								($proj->{TAX_READ_GENUS_HEATMAP_PDF} = $1) =~ s/species/genus/;
								$proj->{TAX_READ_GENUS_HEATMAP_PDF_TXT} = "PDF";
								($proj->{TAX_READ_STRAIN_HEATMAP_PDF} = $1) =~ s/species/strain/;
								$proj->{TAX_READ_STRAIN_HEATMAP_PDF_TXT} = "PDF";
								next;
							}
							if(/<img class="preview_img" data-src="(.*heatmap_DATASET-allReads\.species\.png)" alt="Coverage">/) {
								$proj->{TAX_READ_SPECIES_HEATMAP_PNG} = $1;
								$proj->{TAX_READ_SPECIES_HEATMAP_PNG_TXT} = "PNG";
								($proj->{TAX_READ_GENUS_HEATMAP_PNG} = $1) =~ s/species/genus/;
								$proj->{TAX_READ_GENUS_HEATMAP_PNG_TXT} = "PNG";
								($proj->{TAX_READ_STRAIN_HEATMAP_PNG} = $1) =~ s/species/strain/;
								$proj->{TAX_READ_STRAIN_HEATMAP_PNG_TXT} = "PNG";
								next;
							}
							if(/<iframe style='height:510px;' data-src='(.*radarchart_DATASET_allReads\.species\.html)'>Your browser does not support iframes<\/iframe>/) {
								$proj->{TAX_READ_SPECIES_RADARCHART_HTML} = $1;
								$proj->{TAX_READ_SPECIES_RADARCHART_HTML_TXT} = "HTML";
								($proj->{TAX_READ_GENUS_RADARCHART_HTML} = $1) =~ s/species/genus/;
								$proj->{TAX_READ_GENUS_RADARCHART_HTML_TXT} = "HTML";
								($proj->{TAX_READ_STRAIN_RADARCHART_HTML} = $1) =~ s/species/strain/;
								$proj->{TAX_READ_STRAIN_RADARCHART_HTML_TXT} = "HTML";
								next;
							}
					
							if(/<a data-ajax='false' href='(.*report_summary\.xlsx)'> Summary table <\/a>/) {
								$proj->{TAX_READ_REPORT_SUM} = $1;
								$proj->{TAX_READ_REPORT_SUM_TXT} = "xlsx";
								next;
							}
					
							if(/<a data-ajax='false' href='(.*report_SEQ1_allReads\.xlsx)'> Abundance <\/a>/) {
								$proj->{TAX_READ_REPORT_ABUNDANCE} = $1;
								$proj->{TAX_READ_REPORT_ABUNDANCE_TXT} = "xlsx";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/<div .*data-role='collapsible' .*>/) {
								last;
							}
						}
						next;
					}

					##assembly tools
					if(/<span class="li-report-content-title">GOTTCHA \(bacterial species database\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_PLASMID_PCT} = $4;
											$tool->{LS_ABUNDANCE} = $5;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_GBS_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-gbs'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
								next;
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-b\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GBS_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_GBS_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-b\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GBS_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_GBS_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-b\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_GBS_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_GBS_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if

					if(/<span class="li-report-content-title">GOTTCHA2 \(BacteriaViruses species database\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_LEN} = $4;
											$tool->{LS_LDOC} = $5;
											$tool->{LS_RDOC} = $6;
											$tool->{LS_ABUNDANCE} = $7;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_GBS2_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-gbs2'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-b\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GBS2_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_GBS2_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-b\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GBS2_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_GBS2_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-b\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_GBS2_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_GBS2_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#GOTTCHA viral
					if(/<span class="li-report-content-title">GOTTCHA \(viral species database\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_PLASMID_PCT} = $4;
											$tool->{LS_ABUNDANCE} = $5;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_GVS_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-gvs'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-v\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GVS_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_GVS_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-v\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GVS_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_GVS_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-gottcha-speDB-v\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_GVS_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_GVS_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#GOTTCHA2 viral
					if(/<span class="li-report-content-title">GOTTCHA2 \(viral species database\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_LEN} = $4;
											$tool->{LS_LDOC} = $5;
											$tool->{LS_RDOC} = $6;
											$tool->{LS_ABUNDANCE} = $7;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_GVS2_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-gvs2'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-v\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GVS2_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_GVS2_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-v\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_GVS2_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_GVS2_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-gottcha2-speDB-v\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_GVS2_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_GVS2_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#PanGIA
					if(/<span class="li-report-content-title">PanGIA<\/span><div class="li-report-content">/) {
						#get pangia result list
						my @parts;
						open(PR, $pangiaReport);				
						while(<PR>) {
							chomp;
							@parts = split(/\t/);
							next if($parts[0] ne "species");
							my $cols = @parts;
							my $tool;

							if($cols > 70) {
								#new pangia result
								next if($parts[13] < $pangia_score);
								$tool->{PROJNAME} = $proj->{PROJNAME};
								$tool->{PROJ} = $proj->{PROJ};
								$tool->{PROPROJURL} = $vars->{PROJURL};
								$tool->{LS_TAX} = $parts[1];
								$tool->{LS_TAX_LINK} = 'http://www.ncbi.nlm.nih.gov/genome/?term="'.$parts[1].'"';
								$tool->{LS_READS} = $parts[5];
								$tool->{LS_COV} = $parts[8];
								$tool->{LS_DOC} = $parts[9];
								$tool->{LS_RS} = $parts[7];
								$tool->{LS_SCORE} = $parts[13];
								my $abu = $parts[14];
								$tool->{LS_ABUNDANCE} = $abu*100;
							} else {							
								#old pangia result
								next if($parts[13] < $pangia_score);
								$tool->{PROJNAME} = $proj->{PROJNAME};
								$tool->{PROJ} = $proj->{PROJ};
								$tool->{PROPROJURL} = $vars->{PROJURL};
								$tool->{LS_TAX} = $parts[1];
								$tool->{LS_TAX_LINK} = 'http://www.ncbi.nlm.nih.gov/genome/?term="'.$parts[1].'"';
								$tool->{LS_READS} = $parts[5];
								$tool->{LS_COV} = $parts[16];
								$tool->{LS_DOC} = $parts[11];
								$tool->{LS_RS} = $parts[7];
								$tool->{LS_SCORE} = $parts[13];
								my $abu = $parts[14];
								$tool->{LS_ABUNDANCE} = $abu*100;
							}
																
							push @{$proj->{TAX_TOOL_PANGIA_LOOP}}, $tool;
							push @{$reports_map{$project}{'tax-tool-pangia'}}, $tool;
						}
						close(PR);
						
						#end
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-pangia\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_PANGIA_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_PANGIA_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-pangia\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_PANGIA_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_PANGIA_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-pangia\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_PANGIA_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_PANGIA_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#Metaphlan2
					if(/<span class="li-report-content-title">Metaphlan\d?<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_ABUNDANCE} = $4;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_METAPHLAN_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-metaphlan'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-metaphlan\d?\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_METAPHLAN_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_METAPHLAN_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-metaphlan\d?\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_METAPHLAN_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_METAPHLAN_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-metaphlan\d?\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_METAPHLAN_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_METAPHLAN_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#BWA
					if(/<span class="li-report-content-title">BWA \(reads mapping\)<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_ABUNDANCE} = $4;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_BWA_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-bwa'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-bwa\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_BWA_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_BWA_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-bwa\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_BWA_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_BWA_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-bwa\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_BWA_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_BWA_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if 
					#Kraken mini
					if(/<span class="li-report-content-title">Kraken.*<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_ABUNDANCE} = $4;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_KRAKEN_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-kraken'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-kraken.*\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_KRAKEN_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_KRAKEN_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-kraken.*\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_KRAKEN_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_KRAKEN_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-kraken.*\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_KRAKEN_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_KRAKEN_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
					#Diamond
					if(/<span class="li-report-content-title">Diamond<\/span><div class="li-report-content">/) {
						while(<IN>) {
							chomp;
							s/^\s+|\s+$//g;
							if(/edge-output-tool-summary/ && /table/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse tool summary
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $tool;
											$tool->{PROJNAME} = $proj->{PROJNAME};
											$tool->{PROJ} = $proj->{PROJ};
											$tool->{PROPROJURL} = $vars->{PROJURL};
											$tool->{LS_TAX} = $2;
											$tool->{LS_READS} = $3;
											$tool->{LS_ABUNDANCE} = $4;
											if($tool->{LS_TAX} =~ /<a href='(.*)' target='_blank'>(.*)<\/a>/) {
												$tool->{LS_TAX} = $2;
												$tool->{LS_TAX_LINK} = $1;
											}
											push @{$proj->{TAX_TOOL_DIAMOND_LOOP}}, $tool;
											push @{$reports_map{$project}{'tax-tool-diamond'}}, $tool;
										} 
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}#end inner while
							}

							##get files/figures					
							if(/<a data-ajax='false' href='(.*allReads-diamond\.tree\.svg)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_DIAMOND_TREE_PLOT} = $1;
								$proj->{TAX_TOOL_DIAMOND_TREE_PLOT_TXT} = "SVG";
								next;
							}					
							if(/<a data-ajax='false' href='(.*allReads-diamond\.krona\.html)'>\[fullwindow\]<\/a><\/span>/) {
								$proj->{TAX_TOOL_DIAMOND_KRONA_PLOT} = $1;
								$proj->{TAX_TOOL_DIAMOND_KRONA_PLOT_TXT} = "HTML";
								next;
							}								
							if(/<a data-ajax='false' href='(.*allReads-diamond\.list\.txt)'> Abundance <\/a>/) {
								$proj->{TAX_TOOL_DIAMOND_ABUNDANCE} = $1;
								$proj->{TAX_TOOL_DIAMOND_ABUNDANCE_TXT} = "TXT";
								last;
							}
					
							if(/<a href='#edge-outputfile-dialog' class='edge-outputfile-tree' data-rel='popup' dir-src='.*'> Directory <\/a><\/p>/) {
								last;
							}
							if(/\[Not available\]/) {
								last;
							}
						}#end while
						next;						
					}#end if
				
					##assembly-based
					if(/<span class="li-report-content-title">Assembly-based Community Profiling<\/span><div class="li-report-content">/) {
						while(<IN>) {
							###get stats
							if(/<div id="edge-output-ccp-summary">/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr><td>Total Contigs<\/td><td>(.*)<\/td><td>(.*)<\/td><\/tr>/) {
										$proj->{TAX_ASSEMBLY_TOTALCONTIGS_STATS} = $1;
										$reports_map{$project}{'tax-assembly-total-contigs'} = $1;
										$proj->{TAX_ASSEMBLY_TOTALCONTIGS_BASES} = $2;
										$reports_map{$project}{'tax-assembly-total-contigs-bases'} = $2;
										next;
									}
									if(/<tr><td>Classified Contigs<\/td><td>(.*)<\/td><td>(.*)<\/td><\/tr>/) {
										$proj->{TAX_ASSEMBLY_CLASSCONTIGS_STATS} = $1;
										$reports_map{$project}{'tax-assembly-classified-contigs'} = $1;
										$proj->{TAX_ASSEMBLY_CLASSCONTIGS_BASES} = $2;
										$reports_map{$project}{'tax-assembly-classified-contigs-bases'} = $2;
										next;
									}
									if(/<tr><td>Unclassified Contigs<\/td><td>(.*)<\/td><td>(.*)<\/td><\/tr>/) {
										$proj->{TAX_ASSEMBLY_UNCLASSCONTIGS_STATS} = $1;
										$reports_map{$project}{'tax-assembly-unclassified-contigs'} = $1;
										$proj->{TAX_ASSEMBLY_UNCLASSCONTIGS_BASES} = $2;
										$reports_map{$project}{'tax-assembly-unclassified-contigs-bases'} = $2;
										next;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
							}
							###get ranks
							if(/<div id="edge-output-ccp-rank-by-count-table">/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse ref
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $ref;
											$ref->{PROJNAME} = $proj->{PROJNAME};
											$ref->{PROJ} = $proj->{PROJ};
											$ref->{PROPROJURL} = $vars->{PROJURL};
											$ref->{LS_RANK} = $1;
											$ref->{LS_TOP1} = $2;
											$ref->{LS_TOP2} = $3;
											$ref->{LS_TOP3} = $4;
											$ref->{LS_TOP4} = $5;
											$ref->{LS_TOP5} = $6;
											$ref->{LS_RANK} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP1} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP2} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP3} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP4} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP5} =~ s/<a href.*>(.*)<\/a>/${1}/;
											
											push @{$proj->{TAX_ASSEMBLY_COUNT_RANK_LOOP}}, $ref;
											push @{$reports_map{$project}{'tax-assembly-count-ranks'}}, $ref;
											
										}
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
								next;
							}
							###get ranks
							if(/<div id="edge-output-ccp-rank-by-length-table">/) {
								while(<IN>) {
									chomp;
									if(/<tbody>/) {
										last;
									}
								}
								my $str;
								while(<IN>) {
									chomp;
									s/^\s+|\s+$//g;
									if(/<tr>/) {
										$str = '';
									} elsif(/<\/tr>/) {
										##parse ref
										if($str =~ /<td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td><td>(.*)<\/td>/) {
											my $ref;
											$ref->{PROJNAME} = $proj->{PROJNAME};
											$ref->{PROJ} = $proj->{PROJ};
											$ref->{PROPROJURL} = $vars->{PROJURL};
											$ref->{LS_RANK} = $1;
											$ref->{LS_TOP1} = $2;
											$ref->{LS_TOP2} = $3;
											$ref->{LS_TOP3} = $4;
											$ref->{LS_TOP4} = $5;
											$ref->{LS_TOP5} = $6;
											$ref->{LS_RANK} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP1} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP2} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP3} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP4} =~ s/<a href.*>(.*)<\/a>/${1}/;
											$ref->{LS_TOP5} =~ s/<a href.*>(.*)<\/a>/${1}/;
											
											push @{$proj->{TAX_ASSEMBLY_LEN_RANK_LOOP}}, $ref;
											push @{$reports_map{$project}{'tax-assembly-len-ranks'}}, $ref;
											
										}
									} else {
										$str .= $_;
									}
									if(/<\/tbody>/) {
										last;
									}
								}
								next;
							}
							##figures and files
							if(/<img class="preview_img" data-src="(.*contigClassification_phylum_barplot\.png)" alt="Contigs Classification at Phylum Level">/) {
								$proj->{TAX_ASSEMBLY_LEN_BARPLOT} = $1;
								$proj->{TAX_ASSEMBLY_LEN_BARPLOT_TXT} = "PNG";
								next;
							}
							if(/<img class="preview_img" data-src="(.*contigClassification_phylum_barplot2\.png)" alt="Contigs Classification at Phylum Level">/) {
								$proj->{TAX_ASSEMBLY_COUNT_BARPLOT} = $1;
								$proj->{TAX_ASSEMBLY_COUNT_BARPLOT_TXT} = "PNG";
								next;
							}
							if(/<img class="preview_img" data-src="(.*contigClassification_phylum_scatterplot\.png)" alt="Contigs Classification at Phylum Level">/) {
								$proj->{TAX_ASSEMBLY_LEN_SCATTERPLOT} = $1;
								$proj->{TAX_ASSEMBLY_LEN_SCATTERPLOT_TXT} = "PNG";
								next;
							}
					
							if(/<a data-ajax='false' href='(.*\.contigsClassification\.pdf)'> Report PDF <\/a>/) {
								$proj->{TAX_ASSEMBLY_REPORT_PDF} = $1;
								$proj->{TAX_ASSEMBLY_REPORT_PDF_TXT} = "PDF";
								next;
							}
							if(/<a href='.*' data-src='(.*LCA\.json)'> Table <\/a>/) {
								$proj->{TAX_ASSEMBLY_RESULT_JSON} = $1;
								$proj->{TAX_ASSEMBLY_RESULT_JSON_TXT} = "JSON";
								last;
							}
							
						}
					}
					if(/<div .*data-role='collapsible' .*>/) {
						last;
					}
				}
				next;
			}
			## end Taxonomy Classification

		}
		close IN;

		$proj->{PROPROJURL} = $vars->{PROJURL};
		$proj->{PRORUNNAME} = $vars->{RUNNAME};
		$proj->{PRORUNDESC} = $vars->{RUNDESC};
		$proj->{PRORUNFILES} = $vars->{RUNFILES};
		$proj->{PROSAMPLEMETADATA_NAME} = $vars->{SAMPLEMETADATA_NAME};
		$proj->{PROSAMPLEMETADATA_TYPE} = $vars->{SAMPLEMETADATA_TYPE};
		$proj->{PROSAMPLEMETADATA_HOST} = $vars->{SAMPLEMETADATA_HOST};
		$proj->{PROSAMPLEMETADATA_SOURCE} = $vars->{SAMPLEMETADATA_SOURCE};
		$proj->{PROSAMPLEMETADATA_COLLECTIONDATE} = $vars->{SAMPLEMETADATA_COLLECTIONDATE};
		$proj->{PROSAMPLEMETADATA_LOCATION} = $vars->{SAMPLEMETADATA_LOCATION};
		$proj->{PROSAMPLEMETADATA_SEQCENTER} = $vars->{SAMPLEMETADATA_SEQCENTER};	
		$proj->{PROSAMPLEMETADATA_SEQUENCER} = $vars->{SAMPLEMETADATA_SEQUENCER};	
		$proj->{PROSAMPLEMETADATA_SEQDATE} = $vars->{SAMPLEMETADATA_SEQDATE};
		
		push @{$vars->{RUNINFO_LS}}, $proj;
	}
}

sub create_run_info_csv {
	my $file_name = shift;
	#get header
	my $content;
	if($settings->{'run-name'} eq "on") {
		$content .= "Run Name,";
	}
	if($settings->{'run-desc'} eq "on") {
		$content .= "Description,";
	}
	if($settings->{'run-files'} eq "on") {
		$content .= "Input Files,";
	}
	if($settings->{'sample-metadata'} eq "on") {
		my @metadata_fields = split(",",$settings->{'sample-metadata-selected'});
		foreach my $field (@metadata_fields) {
			if($field eq "sample-name") {
				$content .= "Sample Name,";
			} elsif($field eq "sample-type") {
				$content .= "Sample Type,";
			} elsif($field eq "host") {
				$content .= "Sample Host,";
			} elsif($field eq "isolation-source") {
				$content .= "Isolation Source,";
			} elsif($field eq "collection-date") {
				$content .= "Collection Date,";
			} elsif($field eq "location") {
				$content .= "Sample Location,";
			} elsif($field eq "sequencing-center") {
				$content .= "Sequencing Center,";
			} elsif($field eq "sequencer") {
				$content .= "Sequencer,";
			} elsif($field eq "sequencing-date") {
				$content .= "Sequencing Date,";
			}
		}
	}
	$content =~ s/,$/\n/;
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		if($settings->{'run-name'} eq "on") {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";
		}
		if($settings->{'run-desc'} eq "on") {
			($celltext = $reports_map{$project}{'run-desc'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";
		}
		if($settings->{'run-files'} eq "on") {
			($celltext = $reports_map{$project}{'run-files'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";
		}
		if($settings->{'sample-metadata'} eq "on") {
			my @metadata_fields = split(",",$settings->{'sample-metadata-selected'});
			foreach my $field (@metadata_fields) {
				if($field eq "sample-name") {
					($celltext = $reports_map{$project}{'sample-name'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "sample-type") {
					($celltext = $reports_map{$project}{'sample-type'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "host") {
					($celltext = $reports_map{$project}{'sample-host'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "isolation-source") {
					($celltext = $reports_map{$project}{'isolation-source'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "collection-date") {
					($celltext = $reports_map{$project}{'collection-date'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "location") {
					($celltext = $reports_map{$project}{'sample-location'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "sequencing-center") {
					($celltext = $reports_map{$project}{'sequencing-center'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "sequencer") {
					($celltext = $reports_map{$project}{'sequencer'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				} elsif($field eq "sequencing-date") {
					($celltext = $reports_map{$project}{'sequencing-date'}) =~ s/"/""/g;
					$content .= "\"$celltext\"".",";
				}
			}
		}
		$content =~ s/,$/\n/;
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_raw_reads_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Reads,Total Bases,Mean Read Length\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-raw-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-raw-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-raw-read-length'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_trimmed_reads_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Reads,Reads %, Total Bases,Total Bases %,Mean Read Length,Paired Reads,Paired Reads %,Paired Total Bases,Paired Total Bases %,Unpaired Reads,Unpaied Reads %,Unpaired Total Bases,Unpaired Total Bases %\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-reads-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-bases-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-read-length'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-paired-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-paired-reads-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-paired-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-paired-bases-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-unpaired-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-unpaired-reads-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-unpaired-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-trimmed-unpaired-bases-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_hostrmf_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Reads,Non-host Reads,Non-host Reads %,Host File(s)\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-host-removal-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-host-removal-nonhost-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'preprocess-host-removal-nonhost-reads-percent'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		$celltext = $reports_map{$project}{'preprocess-host-file'};
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_assembly_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Assembly Tool,Number of Contigs,N50,Max Contig Size,Min Contig Size,Total Assembly Size\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-tool'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-contigs'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-n50'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-maxcontigsize'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-mincontigsize'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-assembly-size'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_read_mapping_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Number of Mapped Reads,% of Total Reads,Number of Unmapped Reads,% of Total Reads,Average Fold Coverage\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-mapping-mappedreads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-mapping-mappedpct'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-mapping-unmappedreads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-mapping-unmappedpct'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-mapping-cov'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_annotation_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,CDS,rRNA,tRNA\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-annotation-cds'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-annotation-rrna'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'aa-annotation-trna'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_ref_reads_stats_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Number of Mapped Reads,% of Total Post-QC Reads,Average Fold,Linear Coverage,SNPs,InDels\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-mapped-reads'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-post-qc-reads-pct'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-read-avg-fold'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-read-linear-cov'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-read-snps'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-read-indels'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_ref_contigs_stats_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Number of Mapped Contigs,Proportion,Average Fold,Linear Coverage,Average Identity,SNPs,InDels\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-mapped-contigs'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-proportion-pct'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-contig-avg-fold'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-contig-linear-cov'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-contig-avg-idt-pct'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-contig-snps'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ref-contig-indels'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_ref_reads_refs_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Reference,Name,Length,GC%,Mapped Reads,Base Coverage,Avg Fold,Fold std.,Gaps,Gap bases,SNPs,INDELs\n";	
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'ref-mappedreads-refs'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_NAME}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_LEN}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GC}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_MAPPEDREADS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_BASECOV}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_AVGFOLD}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_FOLDSTD}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GAPS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GAPBASES}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_SNPS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_INDELS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_ref_contigs_refs_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Reference,Name,Length,GC%,Mapped Contigs,Mapped Contigs %,Base Coverage,Avg Fold,Gaps,Gap bases,SNPs,INDELs\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {		
		foreach my $ref (@{$reports_map{$project}{'ref-mappedcontigs-refs'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_NAME}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_LEN}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GC}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_MAPPEDCONTIGS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_MAPPEDCONTIGS_PCT}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_BASECOV}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_AVGFOLD}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GAPS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_GAPBASES}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_SNPS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_REF_INDELS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_read_genus_ranks_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Tool,Reads,% Reads,Rank,Top1,Top2,Top3,Top4,Top5\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'tax-read-genus'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOOL};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS_PCT}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RANK};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP1};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP2};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP3};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP4};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP5};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_read_species_ranks_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Tool,Reads,% Reads,Rank,Top1,Top2,Top3,Top4,Top5\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'tax-read-species'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOOL};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS_PCT}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RANK};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP1};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP2};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP3};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP4};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP5};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_read_strain_ranks_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Tool,Reads,% Reads,Rank,Top1,Top2,Top3,Top4,Top5\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'tax-read-strain'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOOL};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS_PCT}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RANK};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP1};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP2};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP3};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP4};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP5};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_assembly_contigs_stats_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Total Contigs, Total Bases,Classified Contigs,Classified Bases,Unclassified Contigs,Unclassified Bases\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'tax-assembly-total-contigs'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'tax-assembly-total-contigs-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'ax-assembly-classified-contigs'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'tax-assembly-classified-contigs-bases'}) =~ s/,|bp//g;
		$content .= "\"$celltext\"".",";

		($celltext = $reports_map{$project}{'tax-assembly-unclassified-contigs'}) =~ s/,|\s+|bp//g;

		($celltext = $reports_map{$project}{'tax-assembly-unclassified-contigs-bases'}) =~ s/,|\s+|bp//g;
		$content .= "\"$celltext\""."\n";
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_assembly_count_ranks_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Rank,Top1,Top2,Top3,Top4,Top5\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'tax-assembly-count-ranks'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RANK};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP1};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP2};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP3};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP4};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP5};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_assembly_len_ranks_csv {
	my $file_name = shift;
	#get header
	my $content = "Project/Run Name,Rank,Top1,Top2,Top3,Top4,Top5\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{'tax-assembly-len-ranks'}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RANK};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP1};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP2};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP3};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP4};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TOP5};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_tool_gottcha_csv {
	my $file_name = shift;
	my $tool = shift;
	#get header
	my $content = "Project/Run Name,Taxonomy,Reads,% Plasmid,Abundance\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{$tool}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TAX};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_PLASMID_PCT};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_ABUNDANCE};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_tool_gottcha2_csv {
	my $file_name = shift;
	my $tool = shift;
	#get header
	my $content = "Project/Run Name,Taxonomy,Reads,Linear LEN,Linear DOC, Rollup DOC,Abundance\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{$tool}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TAX};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_LEN}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_LDOC};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RDOC};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_ABUNDANCE};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_tool_pangia_csv {
	my $file_name = shift;
	my $tool = shift;
	#get header
	my $content = "Project/Run Name,Taxonomy,Reads,Linear COV,Linear DOC, RS Norm Reads,Score,Abundance\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{$tool}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TAX};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_COV};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_DOC};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_RS};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_SCORE};
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_ABUNDANCE};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}

sub create_tax_tool_other_csv {
	my $file_name = shift;
	my $tool = shift;
	#get header
	my $content = "Project/Run Name,Taxonomy,Reads,Abundance\n";
	#get runs
	my $celltext;
	foreach my $project (@projects) {
		foreach my $ref (@{$reports_map{$project}{$tool}}) {
			($celltext = $reports_map{$project}{'run-name'}) =~ s/"/""/g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_TAX};
			$content .= "\"$celltext\"".",";

			($celltext = $ref->{LS_READS}) =~ s/,|\s+|bp//g;
			$content .= "\"$celltext\"".",";

			$celltext = $ref->{LS_ABUNDANCE};
			$content .= "\"$celltext\""."\n";
		}
	}
	
	open OUT, ">$out_dir/$file_name" || die "failed to open $out_dir/$file_name: $!";
	print OUT $content;
	close OUT;
}





