#!/usr/bin/env perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
require "./edge_user_session.cgi";
use CGI;
my $query = CGI->new;
my $sid = $query->param( 'sid' );
my $fn = $query->param( 'fn' );
my $ln = $query->param( 'ln' );

print "Content-Type: text/html\n\n";

&stringSanitization($sid);
my $valid = &verifySession($sid);

if($valid or ($fn eq "Guest" and $ln eq "EDGE") ) {
   print "true";
} else {
   print "false";
}

sub stringSanitization{
	my $str=shift;
	if($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/){
		print "Invalid characters detected \'$str\'.\n\n";
		exit;
	}
}
