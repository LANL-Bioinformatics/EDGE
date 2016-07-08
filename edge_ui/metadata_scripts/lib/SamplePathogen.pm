#!/user/bin/perl

package SamplePathogen;

sub new {
	my $class = shift;
	my $self = {
		metadata_id 			=> shift,
		name	 			=> shift,
		hosts				=> shift,
		diseases				=> shift,
		top					=> shift,
	};

   	 bless $self, $class;
    	return $self;
}

sub toJson {
	my ($self, $key, $token) = @_;
	my $json = '{"key":"'.$key.'","token":"'.$token.'"';

	$json .= ',"metadata_id":"'.$self->{metadata_id}.'"' if $self->{metadata_id};
	$json .= ',"name":"'.$self->{name}.'"' if $self->{name};
	$json .= ',"hosts":"'.$self->{hosts}.'"' if $self->{hosts};
	$json .= ',"diseases":"'.$self->{diseases}.'"' if $self->{diseases};
	$json .= ',"top":"'.$self->{top}.'"' if $self->{top};
	
	$json .="}";

	return $json;
}

1;
