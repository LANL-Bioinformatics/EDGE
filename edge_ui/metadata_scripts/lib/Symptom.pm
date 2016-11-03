#!/user/bin/perl

package Symptom;

sub new {
	my $class = shift;
	my $self = {
		metadata_id 			=> shift,
		category	 			=> shift,
		symptom				=> shift,
	};

   	 bless $self, $class;
    	return $self;
}

sub toJson {
	my ($self, $key, $token, $lanlOnly) = @_;
	my $json = '{"api_key":"'.$key.'","api_token":"'.$token.'"';
	$json .= ',"lanl_only":"'.$lanlOnly.'"' if $lanlOnly;

	$json .= ',"metadata_id":"'.$self->{metadata_id}.'"' if $self->{metadata_id};
	$json .= ',"category":"'.$self->{category}.'"' if $self->{category};
	$json .= ',"symptom":"'.$self->{symptom}.'"' if $self->{symptom};
	
	$json .="}";

	return $json;
}

1;
