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

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $edgeui_output  = $sys->{edgeui_output};
my ($relpath)    = $edgeui_output =~ /^$edgeui_wwwroot\/(.*)$/;
my $proj_list=&scanProjToList($edgeui_output);

if( $proj_list->{$pname} ){
	my $projDir = $relpath . "/". $proj_list->{$pname};
	chdir $edgeui_wwwroot;
	my $cmd = "cd $edgeui_wwwroot; $EDGE_HOME/scripts/munger/outputMunger_w_temp.pl $projDir $projDir/HTML_Report/report_web.html";
	`$cmd`;
	#print STDERR "$cmd";

	open REP, "$projDir/HTML_Report/report_web.html" or die "Can't open report_web.html: $!";
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

sub scanProjToList {
	my $out_dir = shift;
        my $list;
        opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
        while( defined (my $file = readdir BIN) ) {
                next if $file eq '.' or $file eq '..';
		my $projid;
		my $projCode;
                if ( -d "$out_dir/$file" && -r "$out_dir/$file/config.txt"  ) {
			open ( CONFIG, "$out_dir/$file/config.txt") or die "Cannot open $out_dir/$file/config.txt\n";
			while(<CONFIG>){
				last if (/^\[Down/);
				if (/^projid=(\S+)/){
					$projid=$1;
				}
				if (/^projcode=(\S+)/){
					$projCode=$1;
				}
			}
			close CONFIG;
                        $projid ||= $file;
			$list->{$projid} = $file;
			$list->{$file} = $file;
			$list->{$projCode} = $file;
                }
        }
        closedir(BIN);
        return $list;
}

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


