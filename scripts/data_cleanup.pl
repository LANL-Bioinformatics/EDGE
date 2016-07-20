#! /usr/bin/env perl

use strict;
use File::Find;
use FindBin qw($RealBin);
use Getopt::Long;

if (!@ARGV){print "$0 [-print] <dirs> ...\n ";exit;}

my $printonly;
GetOptions{
	"print"	=>\$printonly
};
# read system params from sys.properties
my @dirs = (@ARGV);

if ($printonly){
	find( \&printOldFiles, @dirs);
}else{
	find( \&unwanted, @dirs);
}
sub unwanted {
	$File::Find::name =~ /\.(sam|bam|fastq|fq|gz|tgz|cache)$/i && 
	$File::Find::name !~ /JBrowse|public/i &&  
	unlink $File::Find::name;
}

sub printOldFiles {
	$File::Find::name =~ /\.(sam|bam|fastq|fq|gz|tgz|cache)$/i && 
	$File::Find::name !~ /JBrowse|public/i &&  
	print $File::Find::name;
}

