#!/user/bin/perl

package SampleMetadata;

sub new {
	my $class = shift;
	my $self = {
		id 					=> shift,
		sample_type 			=> shift,
		host					=> shift,
		host_condition			=> shift,
		gender				=> shift,
		age					=> shift,
		isolation_source		=> shift,
		source_detailed_info	=> shift,
		collection_date		=> shift,
		city					=> shift,
		state				=> shift,
		country				=> shift,
		lat					=> shift,
		lng					=> shift,
		seq_date				=> shift,
		seq_platform			=> shift,
		sequencer			=> shift,
	};

   	 bless $self, $class;
    	return $self;
}

sub setId {
	my ( $self, $id ) = @_;
	$self->{id} = $id if defined($id);
	return $self->{id};
}

sub toJson {
	my ($self, $key, $token, $runHost, $runId) = @_;
	my $json = '{"key":"'.$key.'","token":"'.$token.'"';

	$json .= ',"run_host":"'.$runHost.'"' if $runHost;
	$json .= ',"run_id":"'.$runId.'"' if $runId;

	$json .= ',"id":"'.$self->{id}.'"' if $self->{id};
	$json .= ',"sample_type":"'.$self->{sample_type}.'"' if $self->{sample_type};
	$json .= ',"host":"'.$self->{host}.'"' if $self->{host};
	$json .= ',"host_condition":"'.$self->{host_condition}.'"' if $self->{host_condition};
	$json .= ',"gender":"'.$self->{gender}.'"' if $self->{gender};
	$json .= ',"age":"'.$self->{age}.'"' if $self->{age};
	$json .= ',"isolation_source":"'.$self->{isolation_source}.'"' if $self->{isolation_source};
	$json .= ',"source_detailed_info":"'.$self->{source_detailed_info}.'"' if $self->{source_detailed_info};
	$json .= ',"collection_date":"'.$self->{collection_date}.'"' if $self->{collection_date};
	$json .= ',"city":"'.$self->{city}.'"' if $self->{city};
	$json .= ',"state":"'.$self->{state}.'"' if $self->{state};
	$json .= ',"country":"'.$self->{country}.'"' if $self->{country};
	$json .= ',"lat":"'.$self->{lat}.'"' if $self->{lat};
	$json .= ',"lng":"'.$self->{lng}.'"' if $self->{lng};
	$json .= ',"sequencing_date":"'.$self->{seq_date}.'"' if $self->{seq_date};
	$json .= ',"sequencing_platform":"'.$self->{seq_platform}.'"' if $self->{seq_platform};
	$json .= ',"sequencer":"'.$self->{sequencer}.'"' if $self->{sequencer};
	
	$json .="}";

	return $json;
}

1;
