#!/usr/bin/env perl
#chienchi @ lanl.gov
#20160511

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use JSON;

my $limit=3000;
my $contigsizelimit=700;
my $info;

my $projdir;
my $mode;
my $out;
my $delimit="tab";

GetOptions(
			"mode=s"        => \$mode,
			"project_dir=s" => \$projdir,
			"limit=i"	=> \$limit,
			"delimit=s"     => \$delimit, 
			"size=i"	=> \$contigsizelimit,
			"out=s"         => \$out,
			"help|?"        =>  sub{Usage()} );

sub Usage{
	print <<"END";
	Usage: perl $0 [options] tab_delimited_table 

    Options:
          -project_dir    the project directory
          -mode           'contig' or 'reference' or 'ref_gap' or 'feature_count'
                          the first column of the tab_delimited_table 
                          is contig ID or Reference ID or other ID
          -delimit        tab (default) or comma 
          -limit <NUM>    process # for row for the input table. (default: 3000)
          -out            output filename
END

exit;
}


my $table=$ARGV[0];

my $sep = "\t"; 
$sep = "," if ($delimit =~ /comma/i);

if (!$table){&Usage();}

if ($out){
	unlink $out;
	open (my $fh,">",$out);
	select $fh;
}

my $projname;
if ($projdir){
	$projdir =~ s/\/$//;
	$projname = basename($projdir);
}

open(my $fh, $table) or die "Cannot read $table\n";
my $header = <$fh>;
chomp $header;
my @headers = split /$sep/,$header;
splice @headers, 1, 0,  "NCBI BLAST" if ($mode eq 'contig');
$headers[0]= "CONTIG_ID" if ($mode eq 'contig');
shift @headers if ($mode eq "feature_count");
my ($length_index,$start_index,$end_index);
foreach my $i (0..$#headers){
	my $hash;
	$headers[$i] =~ s/^\"|\"$//g;
	if ($mode eq "feature_count"){
		my @fc_header = split(/[\.\/]/, $headers[$i]);
		if (scalar(@fc_header)>=3){
			#$headers[$i]="$fc_header[-3]"."_"."$fc_header[-2]";
			$headers[$i]="$fc_header[-2]";
		}
	}
	$hash->{title}= $headers[$i];
	$hash->{data}= $headers[$i];
	push @{$info->{columns}},$hash;
	$length_index = $i if ($headers[$i] =~ /length/i);
	$start_index = $i if ($headers[$i] =~ /gap_start/i);
	$end_index = $i if ($headers[$i] =~ /gap_end/i);
}	
my $count=0;
while(<$fh>){
	chomp;
	my $data;
	last if ($limit > 0 && $count >= $limit);

	my @array=split(/$sep/,$_);
	if ($mode eq 'contig'){
		splice @array, 1, 0,  "<a href='#' class='edge-contigBlast' >BLAST</a>";
		my $end = ($length_index)? $array[$length_index] : $contigsizelimit;
		$array[0]="<a href='JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA&loc=$array[0]%3A1..$end' target='_blank'>$array[0]</a>" if ($projname);
		next if ($length_index && $array[$length_index] < $contigsizelimit);
	}elsif ($mode eq 'ref_gap'){
		my $start = ($start_index)? $array[$start_index] : 1;
		my $end = ($end_index)? $array[$end_index]: $contigsizelimit;
		$array[0]="<a href='JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA&loc=$array[0]%3A$start..$end' target='_blank'>$array[0]</a>" if ($projname);
	}
	shift @array if ($mode eq "feature_count");
	foreach my $i (0..$#headers){
		$array[$i] = "" if (!$array[$i]);
		$array[$i] =~ s/^\"|\"$//g;
		$array[$i] = ($array[$i] !~ /\D/ && $array[$i] =~ /\d+\.\d+/)? sprintf("%.4f",$array[$i]):$array[$i];
		$data->{$headers[$i]}=$array[$i];
	}
	push @{$info->{data}},$data;
	$count++;
}

&returnStatus;

sub returnStatus {
	my $json = "{}";
	$json = to_json($info) if $info;
	#$json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
	print $json;
	exit;
}

