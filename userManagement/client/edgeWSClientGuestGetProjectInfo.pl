#!/usr/bin/perl 
no warnings 'all';
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use JSON;

# Package the data in a data structure matching the expected JSON
my %data = (
   email => 'test2@my.com',
   password => 'admin',
   project_id => 100000 
);

# Encode the data structure to JSON
my $json = JSON->new;
my $data = $json->encode(\%data);

# Set the request parameters
my $url = 'http://localhost:8080/userManagementWS/project/getInfo';
my $browser = LWP::UserAgent->new;
my $req = PUT $url;
$req->header('Content-Type' => 'application/json');
$req->header('Accept' => 'application/json');
#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
$req->header( "Content-Length" => length($data) );
$req->content($data);

my $response = $browser->request($req);

print $response->decoded_content."\n";
