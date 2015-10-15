#!/usr/bin/perl
# contig_classifier.pl
# ver 0.1
# 2013/09/25
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

use strict;
use Getopt::Long;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$ENV{PATH}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'line|l=s',
    'use_upper_name',
    'help|?') || &usage();

if ( $opt{help} || !-e $opt{input} ) { &usage(); }

my ($INPUT, $LINE) = ($opt{input}, $opt{line});

#init options
$LINE ||= 20000;

my $time = time;
my $part=0;
my $flag=0;
my $cnt=0;
my $ctg = "";
my $fh;

open INPUT, $INPUT or die "Can't open input file: $!\n";
open $fh, ">$INPUT.part-$part" or die "Can't write to file: $!\n";

while(<INPUT>){
    my @temp = split /\t/, $_;
    my $name = $temp[0];
    $name =~ s/_\d+$// if defined $opt{use_upper_name};
    if( $name ne $ctg && $cnt>=$LINE ){
        $cnt = 0;
        close $fh;
        $part++;
        open $fh, ">$INPUT.part-$part" or die "Can't write to file: $!\n";
    }
    $ctg = $name;
    print $fh $_;
    $cnt++;
}

close $fh;
close INPUT;

sub timeInterval{
    my $now = shift;
    $now = time - $now;
    return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}

sub usage {
print <<__END__;

$0 [OPTIONS] --input <FILE>

    --input   | -i <STRING>  the sequence of contigs in FASTA format

[OPTIONS]

    --outdir | -o <DIR>      output directory
    --line   | -l <NUM>      split over lines
    --use_upper_name         split over contig/assembly name
                             (E.g.: contig100_1 -> contig100)
    --help/h/?               display this help                   

__END__
exit();
}
