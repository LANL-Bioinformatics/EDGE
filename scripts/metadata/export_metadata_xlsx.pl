#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../../lib";
use Excel::Writer::XLSX;
use Getopt::Long;
use File::Basename;
use File::Path;

my $um; #user management
my $out;
my $project_dir_names;
my $usage = qq{
Usage: $0
	Required
		-out        	      output file
		-projects        project_dir_names, separated by comma
};

GetOptions(
		"um=s"        =>  \$um,
		"out=s"        =>  \$out,
		"projects=s"       =>  \$project_dir_names,
		"help|?"           => sub{print "$usage\n";exit;} 
	);

if (!$project_dir_names && !$out){ print "$usage\n";exit;}

mkpath(dirname($out));

#create sheets and wrte headers
my ($n1,$n2,$n3,$n4);
my $workbook = Excel::Writer::XLSX->new( $out );
my $headerFormat = $workbook->add_format();
$headerFormat->set_bold();

my $worksheet1 = $workbook->add_worksheet( "sample_metadata" );
my (@header1, @header2, @header3, @header4);
if($um) {
	 @header1 = ("Project/Run Name","Project/Run Owner","Study Title","Study Type", "Sample Name", "Sample Type", "Host", "Host Condition", "Gender", "Age", "Isolation Source", "Collection Date", "Location", "City", "State", "Country", "Lat", "Lng", "Experiment Title", "Sequencing Center", "Sequencer", "Sequencing Date");
	@header2 = ("Project/Run Name","Project/Run Owner", "Date From", "Date To", "Location", "City", "State", "Country", "Lat", "Lng");
	@header3 = ("Project/Run Name","Project/Run Owner", "Category", "Symptom");
	@header4 = ("Project/Run Name","Project/Run Owner", "Field", "Value");
} else {
	 @header1 = ("Project/Run Name","Study Title","Study Type", "Sample Name", "Sample Type", "Host", "Host Condition", "Gender", "Age", "Isolation Source", "Collection Date", "Location", "City", "State", "Country", "Lat", "Lng", "Experiment Title", "Sequencing Center", "Sequencer", "Sequencing Date");
	@header2 = ("Project/Run Name", "Date From", "Date To", "Location", "City", "State", "Country", "Lat", "Lng");
	@header3 = ("Project/Run Name", "Category", "Symptom");
	@header4 = ("Project/Run Name", "Field", "Value");
}
$n1 = 1;
$worksheet1->write( "A$n1", \@header1, $headerFormat );
$n1 ++;

my $worksheet2 = $workbook->add_worksheet( "travels" );
$n2 =1;
$worksheet2->write( "A$n2", \@header2, $headerFormat );
$n2 ++;

my $worksheet3 = $workbook->add_worksheet( "symptoms" );
$n3 =1;
$worksheet3->write( "A$n3", \@header3, $headerFormat );
$n3 ++;

my $worksheet4 = $workbook->add_worksheet( "other" );
$n4 =1;
$worksheet4->write( "A$n4", \@header4, $headerFormat );
$n4 ++;

#write metadata to sheets
foreach my $proj_dir (split /,/,$project_dir_names){
	my $confFile = "$proj_dir/config.txt";
	my $metadataFile = "$proj_dir/metadata_sample.txt";
	my $conf = &getParams($confFile);
	my $metadata = &getParams($metadataFile);
	my $travelsFile = "$proj_dir/metadata_travels.txt";
	my $symptomsFile = "$proj_dir/metadata_symptoms.txt";
	my $otherFile = "$proj_dir/metadata_other.txt";

	my $proj_name = $conf->{'projname'};
	my $owner = $conf->{'projowner'};
	if(-e $metadataFile) {
		writeSampleMetadata($owner, $proj_name, $metadata);
		writeTravelsMetadata($owner, $proj_name, $travelsFile);
		writeSymptomsMetadata($owner, $proj_name, $symptomsFile);
		writeOtherMetadata($owner, $proj_name, $otherFile);
	}
}

$workbook->close() or die "Error closing file: $!";

sub writeSampleMetadata {
	my $owner = shift;
	my $proj = shift;
	my $meta = shift;
	my @row;
	if($um) {
		@row = ($proj, $owner, $meta->{'study_title'},$meta->{'study_type'},$meta->{'sample_name'},$meta->{'sample_type'},$meta->{'host'},$meta->{'host_condition'},$meta->{'gender'},$meta->{'age'},$meta->{'isolation_source'},$meta->{'collection_date'},$meta->{'location'},$meta->{'city'},$meta->{'state'},$meta->{'country'},$meta->{'lat'},$meta->{'lng'},$meta->{'experiment_title'},$meta->{'sequencing_center'},$meta->{'sequencer'},$meta->{'sequencing_date'});
	} else {
		@row = ($proj, $meta->{'study_title'},$meta->{'study_type'},$meta->{'sample_name'},$meta->{'sample_type'},$meta->{'host'},$meta->{'host_condition'},$meta->{'gender'},$meta->{'age'},$meta->{'isolation_source'},$meta->{'collection_date'},$meta->{'location'},$meta->{'city'},$meta->{'state'},$meta->{'country'},$meta->{'lat'},$meta->{'lng'},$meta->{'experiment_title'},$meta->{'sequencing_center'},$meta->{'sequencer'},$meta->{'sequencing_date'});
	}
	
	$worksheet1->write( "A$n1", \@row );
	$n1 ++;
}

sub writeTravelsMetadata {
	my $owner = shift;
	my $proj = shift;
	my $file = shift;
	if(! -e $file) {
		return;
	}
	#parse file
	open TVL, "$file";
	my ($from, $to, $location, $city, $state, $country, $lat, $lng);
	while(<TVL>){
		chomp;
		next if(/^#/);
		if ( /(.*)=(.*)/ ){
			if ($1 eq "travel-date-from") {
				$from = $2;
			} elsif ($1 eq "travel-date-to") {
				$to = $2;
			} elsif ($1 eq "travel-location") {
				$location = $2;
			} elsif ($1 eq "city") {
				$city = $2;
			} elsif ($1 eq "state") {
				$state = $2;
			} elsif ($1 eq "country") {
				$country = $2;
			} elsif ($1 eq "lat") {
				$lat = $2;
			} elsif ($1 eq "lng") {
				$lng = $2;
				my @row;
				if($um) {
					@row =  ($proj, $owner, $from, $to, $location, $city, $state, $country, $lat, $lng);
				} else {
					@row =  ($proj, $from, $to, $location, $city, $state, $country, $lat, $lng);
				}
				$worksheet2->write( "A$n2", \@row );
				$n2 ++;
			}
		}
	}
	close TVL;
}

sub writeSymptomsMetadata {
	my $owner = shift;
	my $proj = shift;
	my $file = shift;
	if(! -e $file) {
		return;
	}
	#parse file
	open SM, "$file";
	while(<SM>) {
		chomp;
		next if(/^#/);
		if ( /(.*)\t(.*)/ ){
			my @row;
			if($um) {
				@row =  ($proj, $owner, $1, $2);
			} else {
				@row =  ($proj, $1, $2);
			}
			$worksheet3->write( "A$n3", \@row );
			$n3 ++;
		}
	}
	close SM;
}

sub writeOtherMetadata {
	my $owner = shift;
	my $proj = shift;
	my $file = shift;
	if(! -e $file) {
		return;
	}
	#parse file
	open SM, "$file";
	while(<SM>) {
		chomp;
		next if(/^#/);
		if ( /(.*)=(.*)/ ){
			my @row;
			if($um) {
				@row =  ($proj, $owner, $1, $2);
			} else {
				@row =  ($proj, $1, $2);
			}
			$worksheet4->write( "A$n4", \@row );
			$n4 ++;
		}
	}
	close SM;
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
