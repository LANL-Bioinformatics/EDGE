#!/user/bin/perl

package SampleMetadata;

sub new {
	my $class = shift;
	my $self = {
		metadata_id 			=> shift,
		sra_run_accession		=> shift,
		study_id	 			=> shift,
		study_title 			=> shift,
		study_type 			=> shift,
		sample_name 			=> shift,
		sample_type 			=> shift,
		host					=> shift,
		host_condition			=> shift,
		gender				=> shift,
		age					=> shift,
		isolation_source		=> shift,
		collection_date		=> shift,
		location 				=> shift,
		city					=> shift,
		state				=> shift,
		country				=> shift,
		lat					=> shift,
		lng					=> shift,
		experiment_title		=> shift,
		sequencing_center		=> shift,
		sequencer			=> shift,
		sequencing_date		=> shift,
	};

   	 bless $self, $class;
    	return $self;
}

sub toJson {
	my ($self, $key, $token, $runHost, $runId, $lanlOnly) = @_;
	my $json = '{"api_key":"'.$key.'","api_token":"'.$token.'"';

	$json .= ',"run_host":"'.$runHost.'"' if $runHost;
	$json .= ',"run_id":"'.$runId.'"' if $runId;
	$json .= ',"lanl_only":"'.$lanlOnly.'"' if $lanlOnly;

	$json .= ',"metadata_id":"'.$self->{metadata_id}.'"' if $self->{metadata_id};
	$json .= ',"sra_run_accession":"'.$self->{sra_run_accession}.'"' if $self->{sra_run_accession};
	$json .= ',"study_id":"'.$self->{study_id}.'"' if $self->{study_id};
	$json .= ',"study_title":"'.$self->{study_title}.'"' if $self->{study_title};
	$json .= ',"study_type":"'.$self->{study_type}.'"' if $self->{study_type};
	$self->{sample_name} =~ s/"/\\"/g;
	$json .= ',"sample_name":"'.$self->{sample_name}.'"' if $self->{sample_name};
	$json .= ',"sample_type":"'.$self->{sample_type}.'"' if $self->{sample_type};
	$json .= ',"host":"'.$self->{host}.'"' if $self->{host};
	$json .= ',"host_condition":"'.$self->{host_condition}.'"' if $self->{host_condition};
	$json .= ',"gender":"'.$self->{gender}.'"' if $self->{gender};
	$json .= ',"age":"'.$self->{age}.'"' if $self->{age};
	$json .= ',"isolation_source":"'.$self->{isolation_source}.'"' if $self->{isolation_source};
	$json .= ',"collection_date":"'.formatDate($self->{collection_date}).'"' if $self->{collection_date};
	$json .= ',"location":"'.$self->{location}.'"' if $self->{location};
	$json .= ',"city":"'.$self->{city}.'"' if $self->{city};
	$json .= ',"state":"'.$self->{state}.'"' if $self->{state};
	$json .= ',"country":"'.$self->{country}.'"' if $self->{country};
	$json .= ',"lat":"'.$self->{lat}.'"' if $self->{lat};
	$json .= ',"lng":"'.$self->{lng}.'"' if $self->{lng};
	$json .= ',"experiment_title":"'.$self->{experiment_title}.'"' if $self->{experiment_title};
	$json .= ',"sequencing_center":"'.$self->{sequencing_center}.'"' if $self->{sequencing_center};
	$json .= ',"sequencer":"'.$self->{sequencer}.'"' if $self->{sequencer};
	$json .= ',"sequencing_date":"'.formatDate($self->{sequencing_date}).'"' if $self->{sequencing_date};
	
	$json .="}";

	return $json;
}

sub formatDate {
	my $date = shift;
	if($date =~ /^[0-9]{4}$/) {
		$date = "$date-01-01";
	}elsif($date =~ /^[0-9]{4}-[0-9]{2}$/) {
		$date = "$date-01";
	}
	if($date !~ /^[0-9]{4}-[0-9]{2}-[0-9]{2}$/) {
		$date = '';
	}
	return $date;
}

1;
