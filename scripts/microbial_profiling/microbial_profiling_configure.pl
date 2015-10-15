#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use HTML::Template;

print STDERR &usage() if !-e $ARGV[0];

my %opt;
my @tools = split /,/, $ARGV[1];
print STDERR &usage() if scalar @tools < 1;

my $bwaScoreCut = $ARGV[2];
$bwaScoreCut ||= 30;

foreach my $tool ( @tools ){
	next unless $tool;
	$opt{"$tool"} = 1;
}

my $template = HTML::Template->new(filename => $ARGV[0], die_on_bad_params => 1 );
$template->param(BWASCORECUT => $bwaScoreCut);
$template->param( %opt );
print $template->output();

sub usage {
print <<USAGE;
$0 [template.tmpl] [tools] > microbial_profiling_configure.settings.ini 
USAGE
exit;
} 

exit 0;

