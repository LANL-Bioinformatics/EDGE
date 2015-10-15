#!/usr/bin/env perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
require "edge_user_session.cgi";
use CGI;
my $query = CGI->new;
my $sid = $query->param( 'sid' );

my $valid = &verifySession($sid);

print "Content-Type: text/html\n\n";
if($valid) {
   print "true";
} else {
   print "false";
}
