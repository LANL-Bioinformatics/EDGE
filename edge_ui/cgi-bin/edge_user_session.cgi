#!/usr/bin/env perl
use FindBin qw($RealBin);
use lib "../../lib";
use CGI::Session qw/-ip_match/;
use CGI::Session;
use strict;

my $DEBUG=0;

sub createSession {
	my ($user,$pass,$type) = @_;
	my $session = new CGI::Session->new() or return 0;
	my $sid = $session->id;
	$session->expire('+12h');
	$session->param('username',$user);
	$session->param('password',$pass);
	$session->param('type',$type);

	#print STDERR "IP: $ENV{REMOTE_ADDR}\nUsername: $user\nPassword: $pass\nSid: $sid\n" if $DEBUG;

	return $sid;
}

sub verifySession {
	my $sid = shift;
	return 0 unless $sid;
	my $session = CGI::Session->load($sid) or return 0;
	return 0 if $session->is_expired;
	return 0 if $session->is_empty;
	$session->expire('+12h');
	return 1;
}

sub closeSession {
	my $sid = shift;
	return 0 unless $sid;
	return 0 unless $sid =~ /^[\w\d]+$/;
	my $session = CGI::Session->load($sid) or return 0;
	$session->delete();
	#$session->flush();
	my $e = unlink("/tmp/cgisess_$sid");
	return $e ? 0 : 1 ;
}

sub getCredentialsFromSession {
	my $sid = shift;
	return 0 unless $sid;
	my $session = CGI::Session->load($sid) or return 0;
	return 0 if $session->is_expired;
	my $user = $session->param('username');
	my $pass = $session->param('password');
	my $type = $session->param('type');
	return ($user,$pass,$type) if $user && $pass;
}

1;
