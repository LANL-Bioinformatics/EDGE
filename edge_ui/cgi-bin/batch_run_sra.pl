#!/usr/bin/env perl

use strict;
use FindBin qw($RealBin);
use POSIX qw(strftime);
use LWP::Simple; # from CPAN
use JSON qw( decode_json ); # from CPAN

exit if ( $ENV{"REQUEST_METHOD"} );
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
	my $cmd = "/usr/bin/curl $proxy -o $file \"http://www.ebi.ac.uk/ena/data/warehouse/search?query=tax_tree(408169)%20AND%20first_public=$date_string%20AND%20library_strategy=WGS%20AND%20library_selection=RANDOM%20AND%20(library_source=GENOMIC%20OR%20library_source=METAGENOMIC)&result=read_run&fields=run_accession,sample_accession,study_accession,study_title,experiment_title,scientific_name,instrument_model,library_layout,base_count&limit=10000&display=report\" ";
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
        next if ($_ =~ /^\s*$/);
	my @array= split /\t/,$_;
	my $sra_id = $array[0];
	my $study_id = $array[2];
	my $study = $array[3];
	my $experiment = $array[4];
	my $instrument = $array[6];
        next unless($sra_id =~ /[A-Z]{3}\d+/);
	my $projname = "bsve_$sra_id";
	my $projdesc = "$sra_id $study";
	$projdesc =~ s/(['"])/\\$1/g;
	#$projdesc =~ s/\(/\\(/g;
	#$projdesc =~ s/\)/\\)/g;
	next if ( -d "$RealBin/../EDGE_output/$projname");

        #get sample metadata from ebi by sample id
        my $sample_id = $array[1];
        my $cmd = "mkdir -p $ENV{HOME}/sra_runs";
        system($cmd);
        my $tmp = "$ENV{HOME}/sra_runs/bsve_metadata_tmp.txt";
        $cmd = "/usr/bin/curl $proxy -o $tmp \"http://www.ebi.ac.uk/ena/data/warehouse/search?query=accession=$sample_id&result=sample&fields=accession,collection_date,country,description,first_public,isolation_source,location,scientific_name,sample_alias,center_name,environment_material,host,host_status,host_sex&display=report\" ";
	#print "$cmd\n";
	system($cmd);
        open TMP, $tmp;
        my $sequencer = $array[6];
        my ($sampleType, $host, $collectionDate, $city, $state, $country, $lat, $lng,$seqPlatform, $gender, $hostCondition, $source, $sampleName, $center, $seqDate, $location);
        while(<TMP>) {
        	chomp;
		next if($_ =~ /^accession/);
                next if ($_ =~ /^\s*$/);
		my @parts = split /\t/, $_;
	 	$collectionDate = $parts[1];
                if($collectionDate =~ /\//) {
                	my @its = split /\//, $collectionDate;
			$collectionDate = $its[1];	
		}
		$location = $parts[2];
                $country = $location;
                if($country =~ /(.*):\s*(.*?),\s*(.*)\s*/) {
                	$country = $1;
			$city = $2;
			$state = $3;
		}	
                if($country =~ /(.*):\s*(.*)\s*:\s*(.*)\s*/) {
                	$country = $1;
			$city = $2;
		}	
                elsif($country =~ /(.*):\s*(.*)\s*/) {
                	$country = $1;
			$state= $2;
		}
		$sampleName = $parts[3];
		$sampleName = $parts[8] unless $sampleName;
	 	$seqDate = $parts[4];
                if($seqDate =~ /\//) {
                	my @its = split /\//, $seqDate;
			$seqDate = $its[1];	
		}	
                $source = $parts[5];
                $source = $parts[10] unless $source;
                my $latlng = $parts[6];
		$latlng =~ s/^\s+//;
		my @its = split /\s+/, $latlng;
                $lat = $its[0];
		$lng = $its[2];
		if($its[1] eq "S") {
 			$lat = -$lat;
		}
		if($its[3] eq "W") {
			$lng = -$lng;
		}
		$host = $parts[11];
                $sampleType = "environmental";
                my $stype = lc $parts[7];
                if($stype =~ /human/ || lc($host) =~ /human/ || lc($host) =~ /homo/) {
			$sampleType = "human";
		} 
                elsif($stype =~ /mouse|rat|pig|fish|ant|chicken|bee|frog/ || lc($host) =~ /mouse|rat|pig|fish|ant|chicken|bee|frog/) {
			$sampleType = "animal";
		} 
		$center = $parts[9];
		$hostCondition = $parts[12];
		$gender = $parts[13];	
                #get lat,lng from location
                if(!$lat && !$lng && $location) {
                	($lat,$lng) = getLatLong($location);
                }
	}
        close TMP;

        my $cpu = $ARGV[1];
        $cpu = 12 unless $cpu;

	chdir $RealBin;
	$cmd = "$RealBin/edge_submit.cgi \"$projname\" \"$sra_id\" \"$projdesc\"";
	$cmd .= " \"$study_id\" \"$study\" \"$sampleName\" \"$sampleType\" \"$gender\" \"$host\" \"$hostCondition\" \"$source\" \"$collectionDate\" \"$location\" \"$city\" \"$state\" \"$country\" \"$lat\" \"$lng\" \"$experiment\" \"$center\" \"$instrument\" \"$seqDate\" \"$cpu\" \"lanl_only=true\"";
	#print $cmd,"\n";
	system($cmd);
}
close $fh;


#################
sub getLatLong($){
  my ($address) = @_;
  my $format = "json"; #can also to 'xml'
  my $geocodeapi = "https://maps.googleapis.com/maps/api/geocode/";
  my $url = $geocodeapi . $format . "?address=" . $address;
  my $json = get($url);
  my $d_json = decode_json( $json );

  my $lat = $d_json->{results}->[0]->{geometry}->{location}->{lat};
  my $lng = $d_json->{results}->[0]->{geometry}->{location}->{lng};

  return ($lat, $lng);
}

1;
