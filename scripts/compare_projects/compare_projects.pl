#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use File::Basename;
use Getopt::Long;
use lib "$RealBin/../../lib";
use HTML::Template;
use Cwd;
use Data::Dumper;

my $version=0.1;
my $debug=0;

my $workingDir=getcwd();
my $out_dir;;
my $html_outfile; 
my $project_dir_names;
my $tax_tools;
my $mode;
my $runhost;
my $usage = qq{
Usage: $0
	Required
		-out_dir         output directory
		-projects        project_dir_names, separated by comma

	Optional:
		-html_outfile    compare_project.html
		-html_host	 the project host uri
		-tax_tools       taxanomy_classification tools, separated by comma
		-help            show this help
};

GetOptions(
		"out_dir=s"        =>  \$out_dir,
		"html_outfile=s"   =>  \$html_outfile,
		"projects=s"       =>  \$project_dir_names,
		"tax_tools=s"      =>  \$tax_tools,
		"html_host=s"      =>  \$runhost,
		"version"          => sub{print "Version: $version\n";exit;},
		"help|?"           => sub{print "$usage\n";exit;} 
	);

if (!$project_dir_names && !$out_dir){ print "$usage\n";exit;}

my $EDGE_HOME = $ENV{EDGE_HOME}||Cwd::abs_path("$workingDir/../..");
my $edge_www="$EDGE_HOME/edge_ui"; 
## Instantiate the variables
$tax_tools ||= "gottcha-genDB-b,gottcha-speDB-b,gottcha-strDB-b,gottcha-genDB-v,gottcha-speDB-v,gottcha-strDB-v,gottcha2-speDB-b,gottcha2-genDB-v,gottcha2-speDB-v,gottcha2-speDB-e-pln,gottcha2-speDB-e-ptz,gottcha2-speDB-e-fug,bwa,kraken,kraken2,centrifuge,metaphlan,metaphlan2,pangia,diamond";
$html_outfile ||= "$out_dir/compare_project.html";
my $vars;

my $rel_out_dir=  ($out_dir =~ /^$edge_www\/(.*)$/)? $1:$out_dir;
$vars->{OUTPUTDIR} = $rel_out_dir;
$vars->{RUNHOST} = $runhost;

`mkdir -p $out_dir`;
&check_analysis();
my $count = &runMetaComp();
&output_html();

unlink "Rplots.pdf";

if (!$count){
	my $msg =  "No Taxonomy Classification results of selecting projects";
	`echo "$msg" >> "$out_dir/log.txt"`;
}

########
sub pull_fastqCount {
	my $count_file = shift;
	my $list = shift;
	open( my $fh, "$count_file");
	while(<$fh>){
		 chomp;
		 my @temp = split /\t/, $_;
		 if(@temp){
		 	$list->{INPUTREADSNUM} += $temp[1];
			$list->{INPUTREADSBASE} += $temp[2];
		}
	}
	$list->{INPUTMRL} = sprintf ('%.2f',$list->{INPUTREADSBASE}/$list->{INPUTREADSNUM});
	close $fh;
	return $list;
}

sub check_analysis {
	
	# clean up existing fof files
	unlink "$out_dir/$_.list.fof.txt" foreach (split /,/,$tax_tools);

	foreach my $proj_dir (split /,/,$project_dir_names){
		my $list;
		my $projId;
		my $projCode;
		my $projName;
		my $projDesc;
		my $projTaxa_dir;
		my $reads_type;
		my $projTaxa_dir1 = "$proj_dir/ReadsBasedAnalysis/Taxonomy/report";
		my $projTaxa_dir2 = "$proj_dir/ReadsBasedAnalysis/UnmappedReads/Taxonomy/report";

		open (my $fh, "$proj_dir/config.txt") || die "Cannot open $proj_dir config file\n";
		while(<$fh>){
			last if (/^\[Down/);
			$projId = $1 if (/^projid=(\S+)/);
			$projCode = $1 if (/^projcode=(\S+)/);
			$projName = $1 if (/^projname=(\S+)/);
			$projDesc = $1 if (/^projdesc=(\S+)/);
		}
		close $fh;
		if ( ! -e $projTaxa_dir1 && ! -e $projTaxa_dir2 ){
			$list->{PROJCPTOOLS} = "No Taxonomy Classification Result";
		}elsif(-e $projTaxa_dir1){
			$projTaxa_dir = $projTaxa_dir1;
			$reads_type  = "allReads";
		}else{
			$projTaxa_dir = $$projTaxa_dir2;
			$reads_type = "UnmappedReads";
		}
		
		$list->{PROJNAME}=$projName;
		$list->{PROJID}=$projId;
		$list->{PROJDESC}=$projDesc;
		$list->{PROJCODE}=$projCode;

		my @proj_run_taxa_tools;
		foreach my $tool (split /,/,$tax_tools){
			my $abu_list = "$projTaxa_dir/1_$reads_type/$tool/$reads_type-$tool.list.txt";
			`echo "$projName\t$abu_list" >> $out_dir/$tool.list.fof.txt` if ( -s $abu_list);
			push @proj_run_taxa_tools, $tool if ( -e $abu_list);
		}
		$list->{PROJCPTOOLS} = join(",",@proj_run_taxa_tools);
		
		my $fastq_count_file = "$proj_dir/QcReads/fastqCount.txt";
		$list = &pull_fastqCount($fastq_count_file,$list) if (-e $fastq_count_file);
		push @{$vars->{LOOP_PROJSUMMARY}}, $list; 
	}
}

sub runMetaComp {
	my $kingdom;
	my $rank;
	my $cmd;
	my $log = "$out_dir/log.txt";
	my $tool_list;
	my $count=0;
	my %tool_name = (
		'gottcha2-speDB-b' => 'GOTTCHA2 BacteriaViruses Species',
		'gottcha2-genDB-v' => 'GOTTCHA2 Viral Genus', 
		'gottcha2-speDB-v' => 'GOTTCHA2 Viral Species', 
		'gottcha2-speDB-e-pln' => 'GOTTCHA2 Species Plant',
		'gottcha2-speDB-e-ptz' => 'GOTTCHA2 Species Protozoa',
		'gottcha2-speDB-e-fug' => 'GOTTCHA2 Species Fungi'
	);
	foreach my $tool (split /,/,$tax_tools){
		my $title;
		my $plot_filename;
		if ( -s "$out_dir/$tool.list.fof.txt" ){
			if ($tool =~ /gottcha-(\w+)DB-(\w)/){
				$kingdom = ($2 eq "b")? "Bacteria": "Virus";
				$rank = "species" if ($1 eq "spe");
				$rank = "genus" if ($1 eq "gen");
				$rank = "strain" if ($1 eq "str");
				$title = "GOTTCHA Merged plot ($rank)";
				$plot_filename = "$out_dir/${tool}_merged_heatmap";
				$cmd = "Rscript $RealBin/merge_and_plot.R $out_dir/$tool.list.fof.txt $out_dir/$tool.merged_assignments.txt $rank \'$title\' $plot_filename gottcha";
				$cmd .= " 1>>$log 2>\&1 ";
				system($cmd) if (! -s "$plot_filename.pdf");
				$tool_list->{GOTTCHA}->{$rank}->{CPRANK}=$rank;
				if ($kingdom eq "Bacteria"){
					$tool_list->{GOTTCHA}->{$rank}->{CPKINGDOMB} = $kingdom;
					$tool_list->{GOTTCHA}->{$rank}->{CPTOOLBPDF} = "${tool}_merged_heatmap.pdf";
					$tool_list->{GOTTCHA}->{$rank}->{CPTOOLBSVG} = "${tool}_merged_heatmap.svg";
				}else{
					$tool_list->{GOTTCHA}->{$rank}->{CPKINGDOMV} = $kingdom;
					$tool_list->{GOTTCHA}->{$rank}->{CPTOOLVPDF} = "${tool}_merged_heatmap.pdf";
					$tool_list->{GOTTCHA}->{$rank}->{CPTOOLVSVG} = "${tool}_merged_heatmap.svg";
				}
				$count++;
			}else{
				$rank= "species";
				if ($tool =~ /gottcha2-(\w+)DB/){
					$rank = "species" if ($1 eq "spe");
					$rank = "genus" if ($1 eq "gen");
				}
				my $loop = "LOOP_".uc($tool);
				$title = uc($tool)." Merged plot ($rank)";
				$plot_filename = "$out_dir/${tool}_merged_heatmap";
				my $tool_for_Rscript = ($tool =~/kraken/i)?"kraken":($tool=~/gottcha2/i)?"gottcha2":"$tool";
				$tool_for_Rscript =~ s/metaphlan\d/metaphlan/;
				$tool_for_Rscript = "kraken" if ($tool=~/centrifuge/);
				$cmd = "Rscript $RealBin/merge_and_plot.R $out_dir/$tool.list.fof.txt $out_dir/$tool.merged_assignments.txt $rank \'$title\' $plot_filename $tool_for_Rscript";
				$cmd .= " 1>>$log 2>\&1 ";
				system($cmd) if (! -s "$plot_filename.pdf");
				my $tool_res;
				$tool_res->{CPRANK}=$rank;
				$tool_res->{CPTOOLNAME}=($tool_name{$tool})?$tool_name{$tool}:uc($tool);
				$tool_res->{CPTOOLPDF}="${tool}_merged_heatmap.pdf";
				$tool_res->{CPTOOLSVG}="${tool}_merged_heatmap.svg" if ( -e "$out_dir/${tool}_merged_heatmap.svg");
				$tool_res->{CPTOOLANCHOR}="edge-compare-${tool}";
				$vars->{ANCHARLINK} .= "<a data-ajax='false' href='compare_project.html#edge-compare-${tool}'> ${tool} </a>|";
				push @{$vars->{LOOP_CPTOOL}}, $tool_res;
				$count++;
			}
		}
	}
	foreach my $tool ( keys %$tool_list ){
		next if ($tool !~ /gottcha/i);
		my $res = $tool_list->{$tool};
		foreach my $rank ( sort keys %$res ){
			my $rank_res = $res->{$rank};
			my $loop = "LOOP_".uc($tool);
			push @{$tool_list->{$loop}}, $rank_res;
		}
	}
	unshift @{$vars->{LOOP_CPTOOL}}, $tool_list;
	return $count;
}

sub output_html {

	my $template = HTML::Template->new(filename => "$RealBin/compare_projects_html.tmpl", 
													strict => 0,
	                                                die_on_bad_params => 0);
	$template->param(%$vars);

	system("cp -r $EDGE_HOME/edge_ui/css $out_dir/"); 
	system("cp -r $EDGE_HOME/edge_ui/images $out_dir/");  
	system("cp -r $EDGE_HOME/edge_ui/javascript $out_dir/");

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
}

