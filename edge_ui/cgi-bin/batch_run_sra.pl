#!/usr/bin/env perl

use strict;
use FindBin qw($RealBin);
use POSIX qw(strftime);


if (!@ARGV){ print "$0 [ bsve_sra_tsv_dir | file | ebi ]\n"; exit;}

my $date_string = strftime "%Y-%m-%d", localtime;
my $proxy = $ENV{HTTP_PROXY} || $ENV{http_proxy};
$proxy = "--proxy $proxy " if ($proxy);

my $file;
if ( -d $ARGV[0]){
	my $dir = $ARGV[0];
	$file= "$dir/bsve_$date_string.tsv";
}elsif( -f $ARGV[0]){
	$file=$ARGV[0];
}else{
	`mkdir -p $ENV{HOME}/sra_runs`;
	$file = "$ENV{HOME}/sra_runs/bsve_$date_string.tsv";
	my $cmd = "/usr/bin/curl $proxy -o $file \"http://www.ebi.ac.uk/ena/data/warehouse/search?query=tax_tree(408169)%20AND%20first_public=$date_string%20AND%20library_strategy=WGS%20AND%20library_selection=RANDOM%20AND%20(library_source=GENOMIC%20OR%20library_source=METAGENOMIC)&result=read_run&fields=run_accession,study_title,experiment_title,scientific_name,instrument_model,library_layout,base_count&limit=10000&display=report\" ";
	print "Query ebi using REST URL ...\n";
	#print $cmd,"\n";
	system($cmd);
}

if ( ! -e $file  or -z $file){
	print "The $file not exist or empty\n";
}

open (my $fh,$file) or die "Cannot open $file. $!\n";
while(<$fh>){
	chomp;
	next if ($_ =~ /run_accession/);
	my @array= split /\t/,$_;
	my $sra_id = $array[0];
	my $study = $array[1];
	my $projname = "bsve_$sra_id";
	my $projdesc = "$sra_id $study";
	$projdesc =~ s/(['"])/\\$1/g;
	#$projdesc =~ s/\(/\\(/g;
	#$projdesc =~ s/\)/\\)/g;
	next if ( -d "$RealBin/../EDGE_output/$projname");
	chdir $RealBin;
	my $cmd = "$RealBin/edge_submit.cgi $projname $sra_id \"$projdesc\"";
	#print $cmd,"\n";
	system($cmd);
}
close $fh;
