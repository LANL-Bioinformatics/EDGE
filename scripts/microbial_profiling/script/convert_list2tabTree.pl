#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage;

my $col = $ARGV[0];

$|=1;

# load taxonomy
loadTaxonomy();

my %major_level = (
	'root'         => 0,
	'superkingdom' => 15,
	'phylum'       => 25,
	'class'        => 35,
	'order'        => 45,
	'family'       => 55,
	'genus'        => 65,
	'species'      => 75,
    'strain'       => 85,
    #'replicon'     => 95
);

my %level = (
    'root'             => 0,
    'superkingdom'     => 15,
    'kingdom'          => 16,
    'subkingdom'       => 17,
    'superphylum'      => 24,
    'phylum'           => 25,
    'subphylum'        => 26,
    'infraclass'       => 33,
    'superclass'       => 34,
    'class'            => 35,
    'subclass'         => 36,
    'infraorder'       => 43,
    'superorder'       => 44,
    'order'            => 45,
    'suborder'         => 46,
    'parvorder'        => 47,
    'superfamily'      => 54,
    'family'           => 55,
    'subfamily'        => 56,
    'tribe'            => 58,
    'genus'            => 65,
    'subgenus'         => 66,
    'species group'    => 73,
    'species subgroup' => 74,
    'species'          => 75,
    'subspecies'       => 76,
    'strain'           => 85,
    'no rank'          => 90,
    'unclassified'     => 100
);

my $taxa;
my $count=0;

while(<STDIN>)
{
	chomp;

    if( /^LEVEL/ && !$col ){
        $col=0;
        my @header_titles = split /\t/, $_;
        for my $header_line (@header_titles) {
            last if $header_line =~ m/^(TAXID|TAXAID|TAXOID|TAXA_ID)$/;
            $col++;
        }
        next;
    }

    my @fields = split /\t/, $_;

    my $taxID = $fields[$col];
	my $name = $fields[1];
	$name =~ s/^ +//;
    
    my $rank = getTaxRank($taxID);
    my $depth = getTaxDepth($taxID);
    
    $rank = "unclassified" if $taxID == 0;
    $rank = "root" if $taxID == 1;
    $rank = "strain" if $rank eq "no rank" && $depth>6;

    #transfer minor rank to major rank
    my $close_major_rank = $rank;
    if( !defined $major_level{$rank} ){
        foreach my $temp ( sort {$major_level{$a}<=>$major_level{$b}} keys %major_level ){
            last if $level{$rank} > $major_level{$temp};
            $close_major_rank = $temp;
        }
    }

    my @lineage;
    foreach my $rank ( sort {$major_level{$a}<=>$major_level{$b}} keys %major_level ){
        my $rname = taxid2rank($taxID,$rank);
        $rname = $name if $rank eq "strain";
        $rname = "no rank" unless $rname;
        push @lineage, $rname;
        last if $rank eq $close_major_rank;
    }

    my $lineage_idx = join "\t", @lineage;
    $lineage_idx =~ s/([\(\),:])//g;
    $taxa->{$lineage_idx} ||= 0;
    $taxa->{$lineage_idx} += $fields[2];
}

foreach my $lineage_idx ( sort keys %$taxa )
{
    my $score = $taxa->{$lineage_idx};
    print "$score\t$lineage_idx\n";
}
