#! /usr/bin/env perl

use strict;
use File::Find;
use FindBin qw($RealBin);
use CGI qw(:standard);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $keep_days	= $sys->{edgeui_proj_store_days};
my $out_dir     = $sys->{edgeui_output};
my $edgeui_wwwroot = $sys->{edgeui_wwwroot};
my $input_dir     = $sys->{edgeui_input};
exit if ( $ENV{"REQUEST_METHOD"} );
exit if (!$keep_days or $keep_days =~ /\D+/);
#my @dirs = ($out_dir,$input_dir);
my $ncbi_tmp = "$edgeui_wwwroot/ncbi";
my @dirs = ($out_dir);
push @dirs, $ncbi_tmp if (-d $ncbi_tmp);

my $keep_secs = $keep_days * 24 * 60 * 60;
my $now = time(); 

if (@ARGV){
	find( \&printOldFiles, @dirs);
}else{
	find( \&unwanted, @dirs);
}

sub unwanted {
	$File::Find::name =~ /\.(sam|bam|fastq|fq|gz|tgz|cache)$/i && 
	$File::Find::name !~ /JBrowse|public/i &&  
	($now-(stat $_)[9]) > $keep_secs &&  
	unlink $File::Find::name;
}

sub printOldFiles {
	$File::Find::name =~ /\.(sam|bam|fastq|fq|gz|tgz|cache)$/i && 
	$File::Find::name !~ /JBrowse|public/i &&  
	($now-(stat $_)[9]) > $keep_secs &&  
	print $File::Find::name;
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
