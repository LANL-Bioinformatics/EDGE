#!/user/bin/perl

package Travel;

sub new {
	my $class = shift;
	my $self = {
		metadata_id 			=> shift,
		dates	 			=> shift,
		location 				=> shift,
		city					=> shift,
		state				=> shift,
		country				=> shift,
		lat					=> shift,
		lng					=> shift,
	};

   	 bless $self, $class;
    	return $self;
}

sub toJson {
	my ($self, $key, $token, $lanlOnly) = @_;
	my $json = '{"api_key":"'.$key.'","api_token":"'.$token.'"';
	$json .= ',"lanl_only":"'.$lanlOnly.'"' if $lanlOnly;

	$json .= ',"metadata_id":"'.$self->{metadata_id}.'"' if $self->{metadata_id};
	$json .= ',"dates":"'.$self->{dates}.'"' if $self->{dates};
	$json .= ',"location":"'.$self->{location}.'"' if $self->{location};
	$json .= ',"city":"'.$self->{city}.'"' if $self->{city};
	$json .= ',"state":"'.$self->{state}.'"' if $self->{state};
	$json .= ',"country":"'.$self->{country}.'"' if $self->{country};
	$json .= ',"lat":"'.$self->{lat}.'"' if $self->{lat};
	$json .= ',"lng":"'.$self->{lng}.'"' if $self->{lng};
	
	$json .="}";

	return $json;
}

1;
