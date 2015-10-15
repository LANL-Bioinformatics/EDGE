#!/usr/bin/perl 
no warnings 'all';
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use JSON;

# Package the data in a data structure matching the expected JSON

# Set the request parameters
my $url = 'http://localhost:8080/userManagementWS/user/hello';
my $browser = LWP::UserAgent->new;
my $req = GET $url;
#$req->header('Accept' => 'application/json');
#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning

my $response = $browser->request($req);

print $response->decoded_content."\n";
