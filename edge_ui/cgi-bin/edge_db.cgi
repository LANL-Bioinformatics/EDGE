#!/usr/bin/env perl

use strict;
use CGI qw(:standard);
use DBI;
use FindBin qw($RealBin);
require "./edge_user_session.cgi";

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $action = $opt{'action'}|| $ARGV[0];
my $name = $opt{'name'} || $ARGV[1];
my $newName = $opt{'new-name'} || $ARGV[2];
my $sampleType = $opt{'sample-type'} || $ARGV[2];
my $proj = $opt{'proj'};

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

&stringSanitization(\%opt);

##db settings
my $host=$sys->{edge_dbhost};
my $db=$sys->{edge_dbname};
my $user=$sys->{edge_dbuser};
my $pw=$sys->{edge_dbpasswd};
# PERL MYSQL CONNECT()
my $dbh   = DBI->connect ("DBI:mysql:database=$db:host=$host", $user, $pw,{PrintError =>0,}) 
                          or die "Can't connect to database: $DBI::errstr\n";
my $stmt;
my $sql;
my $exec;
my $out;

if ($action eq "study-list"){
	$sql = "select name from studies order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "study-type-list"){
	$sql = "select name from study_types order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "seq-center-list"){
	$sql = "select name from seq_centers order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "sequencer-list"){
	$sql = "select name from sequencers order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "isolation-source-list"){
	my $table_prefix = "host";
	if(lc($sampleType) eq "environmental") {
		$table_prefix = "nonhost";
	}
	$sql = "select name from ".$table_prefix."_isolation_source order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "pg-host-list"){
	$sql = "select name from animal_hosts order by name";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	while ( my ($name) = $stmt->fetchrow_array( ) )  {
		 $out .= "$name\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
} elsif ($action eq "symptom-list"){
	$sql = "select c.id, c.name, s.name from symptom_categories c, symptoms s where s.cat_id=c.id";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	my ($id, $cat, $symptom);
	my $old_id;
	my $cnt = 0;
	while ( ($id, $cat, $symptom) = $stmt->fetchrow_array( ) )  {
		 if($id ne $old_id) {
			if($old_id) {
				#not first one
				$out .= "   </fieldset>\n";
				$out .= "</div>\n";
			}
			$old_id = $id;
			$cnt =0;
			$out .= "<div class='ui-field-contain' id='symptom-$id'>\n";
			$out .= "   <fieldset data-role='controlgroup' data-mini='true' data-type='horizontal'>\n";
			$out .= "   <legend>$cat</legend>\n";
			$out .= "   <input type='hidden' name='metadata-symptom-cat' id='metadata-symptom-cat-$id' value='$cat'>\n";
			$out .= "   <input type='hidden' name='metadata-symptom-catID' id='metadata-symptom-catID-$id' value='$id'>\n";
		} 
		$cnt ++;
		$out .= "   <input type='checkbox' name='metadata-symptom-$id' id='metadata-symptom-$id-$cnt' value='$symptom'>\n";
		$out .= "   <label for='metadata-symptom-$id-$cnt'>$symptom</label>\n";
	}
	if($old_id) {
		#not first one
		$out .= "   </fieldset>\n";
		$out .= "</div>\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
}  elsif ($action eq "symptom-edit-list"){
	#get checked symptoms
	my %symptoms=();
	my $symptomFile = $sys->{edgeui_output}."/$proj/metadata_symptoms.txt";

	if(-e $symptomFile) {
		open SYM, $symptomFile or die "Can't open $symptomFile $!";
	 	while(<SYM>){
	      		chomp;
			next if(/^#/);
	     		if ( /(.*)\t(.*)/ ){
				$symptoms{$1}{$2} = 1;
			}
		}
		close SYM;
	}


	$sql = "select c.id, c.name, s.name from symptom_categories c, symptoms s where s.cat_id=c.id";
	#prepare the query
	$stmt = $dbh->prepare( $sql);

	#execute the query
	$stmt->execute( );
	my ($id, $cat, $symptom);
	my $old_id;
	my $cnt = 0;
	while ( ($id, $cat, $symptom) = $stmt->fetchrow_array( ) )  {
		 if($id ne $old_id) {
			if($old_id) {
				#not first one
				$out .= "   </fieldset>\n";
				$out .= "</div>\n";
			}
			$old_id = $id;
			$cnt =0;
			$out .= "<div class='ui-field-contain' id='symptom-$id'>\n";
			$out .= "   <fieldset data-role='controlgroup' data-mini='true' data-type='horizontal'>\n";
			$out .= "   <legend>$cat</legend>\n";
			$out .= "   <input type='hidden' name='metadata-symptom-cat' id='metadata-symptom-cat-$id' value='$cat'>\n";
			$out .= "   <input type='hidden' name='metadata-symptom-catID' id='metadata-symptom-catID-$id' value='$id'>\n";
		} 
		$cnt ++;
		if($symptoms{$cat}{$symptom}) {
			$out .= "   <input type='checkbox' name='metadata-symptom-$id' id='metadata-symptom-$id-$cnt' value='$symptom' checked>\n";
		} else {
			$out .= "   <input type='checkbox' name='metadata-symptom-$id' id='metadata-symptom-$id-$cnt' value='$symptom'>\n";
		}
		$out .= "   <label for='metadata-symptom-$id-$cnt'>$symptom</label>\n";
	}
	if($old_id) {
		#not first one
		$out .= "   </fieldset>\n";
		$out .= "</div>\n";
	}
	warn "Problem in retrieving results", $stmt->errstr( ), "\n" if $stmt->err( );

	$stmt->finish();
}   elsif ($action eq "travel-edit-list"){
	my $travelFile = $sys->{edgeui_output}."/$proj/metadata_travels.txt";

	if(-e $travelFile) {
		open TVL, $travelFile or die "Can't open $travelFile $!";
		my ($from, $to, $location, $city, $state, $country, $lat, $lng);
		my $travels=0;
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
					$travels ++;
					$out .= '<div id="metadata-travel-'.$travels.'">';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label for="remove-travel">Travel #: </label>';
					$out .= '<a href="#" class="ui-btn ui-btn-icon-left ui-icon-delete ui-mini" id="remove-travel">Remove this travel</a>';
					$out .= '</div>';
								
					$out .= '<div class="ui-field-contain">';
					$out .= '<label>From</label>';
					$out .= '<input  data-role="date" type="text" data-mini="true" data-clear-btn="false" name="metadata-travel-date-f" id="metadata-travel-date-f-'.$travels.'" maxlength="30" value="'.$from.'" class="metadata-travel-date">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label>To</label>';
					$out .= '<input  data-role="date" type="text" data-mini="true" data-clear-btn="false" name="metadata-travel-date-t" id="metadata-travel-date-t-'.$travels.'" maxlength="30" value="'.$to.'" class="metadata-travel-date">';
					$out .= '</div>';
					$out .= '<div id="metadata-travel-geo-'.$travels.'">';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label>Location</label>';
					$out .= '<input name="metadata-travel-location" id="geocomplete-travel-'.$travels.'" type="text" placeholder="Type in an address to let system auto fill the location fields below" value="'.$location.'">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label></label>';
					$out .= '<input name="locality" id="metadata-travel-city-'.$travels.'" data-mini="true" data-clear-btn="false" type="text" placeholder="City" value="'.$city.'">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label></label>';
					$out .= '<input name="administrative_area_level_1"  id="metadata-travel-state-'.$travels.'" data-mini="true" data-clear-btn="false" type="text" placeholder="State" value="'.$state.'">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label></label>';
					$out .= '<input name="country" id="metadata-travel-country-'.$travels.'" data-mini="true" data-clear-btn="false" type="text" placeholder="Country" value="'.$country.'">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label></label>';
					$out .= '<input name="lat" id="metadata-travel-lat-'.$travels.'" data-mini="true" data-clear-btn="false" type="text" placeholder="Latitude" value="'.$lat.'">';
					$out .= '</div>';
					$out .= '<div class="ui-field-contain">';
					$out .= '<label></label>';
					$out .= '<input name="lng" id="metadata-travel-lng-'.$travels.'" data-mini="true" data-clear-btn="false" type="text" placeholder="Longitude" value="'.$lng.'">';
					$out .= '</div>';
					$out .= '</div>';
					$out .= '<br><br>';
					$out .= '</div>';
				}
			}
		}
		close TVL;
	}
} elsif($action eq "study-title-add") {
	#get id if exists
	$sql = "select id from studies where name='$name'";
	#prepare the query
	$stmt = $dbh->prepare( $sql);
	#execute the query
	$stmt->execute( );
	
	if($stmt->rows > 0) {
		$out = $stmt->fetchrow_array( );
	} else {
		# insert data into the links table
		$sql = "INSERT INTO studies(name)  VALUES(?)";
	 
		$stmt = $dbh->prepare($sql);
	  	$stmt->execute($name);
		$out = $dbh->{mysql_insertid}
	}
	$stmt->finish();
} elsif($action eq "study-type-add") {
	# insert data into the links table
	$sql = "INSERT INTO study_types(name)  VALUES(?)";
 
	$stmt = $dbh->prepare($sql);
  	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "seq-center-add") {
	# insert data into the links table
	$sql = "INSERT INTO seq_centers(name)  VALUES(?)";
 
	$stmt = $dbh->prepare($sql);
  	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "sequencer-add") {
	# insert data into the links table
	$sql = "INSERT INTO sequencers(name)  VALUES(?)";
 
	$stmt = $dbh->prepare($sql);
  	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "isolation-source-add") {
	my $table_prefix = "host";
	if(lc($sampleType) eq "environmental") {
		$table_prefix = "nonhost";
	}
	# insert data into the links table
	$sql = "INSERT INTO ".$table_prefix."_isolation_source(name)  VALUES(?)";
 
	$stmt = $dbh->prepare($sql);
  	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "animal-host-add") {
	# insert data into the links table
	$sql = "INSERT INTO animal_hosts(name)  VALUES(?)";
 
	$stmt = $dbh->prepare($sql);
  	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "run-add") {
	# insert data into the links table
	$sql = "INSERT INTO runs(name)  VALUES(?)";
	 
	$stmt = $dbh->prepare($sql);
	$stmt->execute($name);
	$stmt->finish();
	$out = $dbh->{mysql_insertid}
} elsif($action eq "run-delete") {
	# insert data into the links table
	$sql = "delete from runs where id=?";
	 
	$stmt = $dbh->prepare($sql);
	$stmt->execute($name);
	$stmt->finish();
} elsif($action eq "run-update") {
	# insert data into the links table
	$sql = "update runs set name='$newName' where id=?";
	 
	$stmt = $dbh->prepare($sql);
	$stmt->execute($name);
	$stmt->finish();
} 

# disconnect from the MySQL database
$dbh->disconnect();

if($action =~ /list$/) {
	print "Content-Type: text/html\n\n", $out;
} else {
	print $out;
}
exit;

#############

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			$flag=1;
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
	die "Incorrect system file\n" if (!$flag);
	return $sys;
}
sub stringSanitization{
	my $opt=shift;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		if($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/){
			print "Content-Type: text/html\n\n", "Invalid characters detected \'$str\'.\n\n";
			exit;
        }
	}
}

1;
