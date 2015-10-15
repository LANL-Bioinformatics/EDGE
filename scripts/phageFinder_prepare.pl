#!/usr/bin/env perl
# Purpose: prepare files for phage finder. 
# This script takes a GFF file from Prokka as input, and produces a
# phage_finder_info.txt (protein table)
# Written by Chien-Chi Lo
# 16 Oct 2014

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $outDir;
my $version=0.1;
GetOptions(
            "o=s"              => \$outDir,
            "version"          => sub{print "Version $version\n";exit;},
            "help|?"           => sub{Usage()} );


if (@ARGV != 2) {&Usage();} 
unless ( -e $ARGV[0] ) { print "GFF File not exist\n"; &Usage();}
unless ( -e $ARGV[1] ) { print "Genome/Contig fasta file not exist\n"; &Usage();}
open(my $fh, $ARGV[0]) or die "Cannot open GFF file\n";

my %len;
my $cds_count=0;
my %id_map;
my $id_map_file="$outDir/id_map.txt";
my $seq_id="Sequence0000001";

## rename fasta file to mapped id
my $new_fasta="$outDir/Assembly.con";
open(my $ofh, ">$new_fasta") or die "Cannot write $new_fasta\n";
open(my $ffh, $ARGV[1]) or die "Cannot open Fasta file\n";
open (my $id_fh, ">$id_map_file") or die "Cannot write $id_map_file\n";
my ($id,$seq);
while(<$ffh>)
{
    chomp;
    if(/^>(\S+)/)
    {
        if ($seq){
           $len{$id}=length($seq);
        }
	$id = $1;
        $id_map{$id}=$seq_id;
        print $id_fh "$seq_id\t$id\n";
        print $ofh ">$seq_id\n";
        $seq_id++;
        $seq="";
    }else{
	$seq .= $_;
    	print $ofh $_,"\n";
    }
}
$len{$id}=length($seq) if ($seq);

close $ffh;
close $id_fh;
close $ofh;

## prepare phage_finder_info file
open (my $ph_fh, ">$outDir/phage_finder_info.txt") or die "Cannot write $outDir/phage_finder_info.txt\n";
while (<$fh>)  # each LOCUS
{
    chomp;
    if (/#sequence-region/)
    {
        my ($tmp, $region_id, $start, $end)=split/\s+/,$_;
        $len{$region_id}=$end-$start+1;
    }
    else
    {
        my ($id,$source,$type,$start,$end,$score,$strand,$phase,$Attributes)=split /\t/,$_;
        if (defined $type and $type eq "CDS")
        {
            my $region_len = $len{$id};
            my %annotations=map { split /=/;} split /;/,$Attributes;
            my $product = $annotations{"product"} || $annotations{"Note"} ||  $annotations{"function"} || "Unknown" ;
            my $locus_tag = $annotations{"locus_tag"} || $annotations{"Name"} || "";
            $product =~ s/\%2C/,/g;
            $product =~ s/\%3B/;/g;
            print $ph_fh join("\t",$id_map{$id},$region_len,$locus_tag,$start,$end,$product),"\n";
            $cds_count++;
        }    
    }    
}
close $ph_fh;
close $fh;


sub Usage 
{
    print <<"END";
    Usage: perl $0 -o outDir GFF_file Fasta_file
    Version $version
    -o      Output directory.
END
    exit;
}
