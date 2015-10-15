#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDF::API2;

if(!@ARGV){print "$0 <pdf_file>\n"; exit;}

my $pdf = PDF::API2->open( "$ARGV[0]");
my $count = $pdf->pages();

print $count,"\n";

$pdf->end;

