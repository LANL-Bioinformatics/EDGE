#!/usr/bin/env perl
#
# Po-E (Paul) Li
# Los Alamos National Lab.
# 2014-08-14
#

use strict;
#use lib "/Users/paulli/perl5/lib/perl5";
#use JSON;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw($RealBin);
#use Data::Dumper;

my $cgi   = CGI->new;
my %opt   = $cgi->Vars();
my $pname = $opt{proj};
$pname ||= $ARGV[0];
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";

# read system params from config template
my $config_tmpl = "$RealBin/edge_config.tmpl";
my $sys         = &getSysParamFromConfig($config_tmpl);

my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};
my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;

if( -d "$edgeui_output/$pname" ){
	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/munger/outputMunger_w_temp.pl $relpath/$pname $relpath/$pname/HTML_Report/report_web.html";
	`$cmd`;
	print STDERR "$cmd";

	open REP, "$edgeui_output/$pname/HTML_Report/report_web.html" or die "Can't open report_web.html: $!";
	my $pr=0;
	my @htmls;
	while(<REP>){
		last if /<!-- \/content -->/;
		push @htmls, $_ if $pr;
		$pr=1 if /id='edge-content-report'/;
	}

	close REP;

	my $html = join "", @htmls;
	#$html =~ s/>\s+</></mg;
	#$html =~ s/\t/ /mg;
	#$html =~ s/ {2,}//mg;
	#$html =~ s/\n/ /mg;

	print "Content-Type: text/html\n\n",
		  $html;
	
		  #open LOG, ">/tmp/edge_report$$.log";
		  #print LOG $test;
		  #close LOG;
}

exit;

######################################################

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			while(<CONF>){
				chomp;
				last if /^\[/;
				if ( /^([^=]+)=([^=]+)/ ){
					$sys->{$1}=$2;
				}
			}
		}
		last;
	}
	close CONF;
	return $sys;
}


