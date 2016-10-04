#!/usr/bin/perl 
no warnings 'all';
use strict;
use Getopt::Long;
use LWP::UserAgent;
use HTTP::Request::Common;
use JSON;
use Term::ReadKey;

my $url = 'http://localhost:8080/userManagementWS/user/admin/register';

print "Enter database username: ";
my $http_auth_username=<STDIN>;
chomp $http_auth_username;
ReadMode('noecho');
print "Enter database password: ";
my $http_auth_password=<STDIN>;
chomp $http_auth_password;
ReadMode(0);
print "\n";
my ($admin_password, $admin_email, $admin_firstname, $admin_lastname);
GetOptions(
	'e=s', \$admin_email,
	'fn=s', \$admin_firstname,
	'ln=s', \$admin_lastname,
);
usage() unless $admin_email;
ReadMode('noecho');
print "Enter EDGE Password: ";
$admin_password=<STDIN>;
chomp $admin_password;
ReadMode(0);
print "\n";


# Package the data in a data structure matching the expected JSON
my %data = (
   email => $admin_email,
   password => $admin_password,
   firstname => $admin_firstname,
   lastname => $admin_lastname,
);

# Encode the data structure to JSON
my $json = JSON->new;
my $data = $json->encode(\%data);
#print "$data\n";

# Set the request parameters
my $browser = LWP::UserAgent->new;

my $req = POST $url;
$req->authorization_basic($http_auth_username, $http_auth_password);
$req->header('Content-Type' => 'application/json');
$req->header('Accept' => 'application/json');
#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
$req->header( "Content-Length" => length($data) );
$req->content($data);

my $response = $browser->request($req);

print $response->decoded_content."\n";

#####################
#usage
#####################
sub usage {
print<<'end';
Usage:   createAdminAccount.pl -e <email> -p <password> -fn <first name> -ln <last name>
Required options:
        -e
        -p
end
die;
}

1;

