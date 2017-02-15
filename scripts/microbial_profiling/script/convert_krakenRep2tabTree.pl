#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage;

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

#  23.85  1564844 1564844 U   0   unclassified
#  76.15  4997221 7839    -   1   root
#  76.03  4989140 105 -   131567    cellular organisms
#  75.47  4952327 9220    -   2       Bacteria
#  28.79  1889195 0   P   1297          Deinococcus-Thermus
#  28.79  1889195 661 C   188787          Deinococci
#  28.78  1888527 106 O   118964            Deinococcales
#  28.78  1888419 0   F   183710              Deinococcaceae
#  28.78  1888419 14154   G   1298                  Deinococcus
#  28.56  1873957 0   S   1299                    Deinococcus radiodurans
#  28.56  1873957 1873957 -   243230                    Deinococcus radiodurans R1
#   0.00  0   0   -   1408434                   Deinococcus radiodurans ATCC 13939
#   0.00  145 0   S   502394                  Deinococcus gobiensis
#   0.00  145 145 -   745776                    Deinococcus gobiensis I-0
#   0.00  64  0   S   55148                   Deinococcus proteolyticus
#   0.00  64  64  -   693977                    Deinococcus proteolyticus MRP
#   0.00  38  0   S   310783                  Deinococcus deserti

while(<STDIN>)
{
	chomp;
	my @fields = split /\t/, $_;
    next if $fields[1] == 0;
    next if $fields[2] == 0;

    my $taxID = $fields[4];
	my $name = $fields[5];
	$name =~ s/^ +//;
    
    my $rank = getTaxRank($taxID)||'no rank';
    my $depth = getTaxDepth($taxID)||'10';
    
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
