#!/user/bin/perl

package EdgeSite;

sub new {
	my $class = shift;
	my $self = {
		organization			=> shift,
		acronym				=> shift,
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
	my $self = shift;

	my $json = '{"organization":"'.$self->{organization}.'"';
	$json .= ',"acronym":"'.$self->{acronym}.'"' if $self->{acronym};
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
