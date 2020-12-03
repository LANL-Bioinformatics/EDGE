#! /usr/bin/perl -w
use strict;
use lib (`ktGetLibPath`);
use KronaTools;

$|=1;

# load taxonomy
my $time = time;
my $period = &timeInterval($time);
loadTaxonomy();
$period = &timeInterval($time);

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
	'replicon'     => 95
);

my %level = (
    'root'             => 0,
    'superkingdom'     => 14,
    'kingdom'          => 15,
    'subkingdom'       => 16,
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
    $count++;
    if( $count%5000 == 0 ){
        $period = &timeInterval($time);
        my $number = $count;
        $number =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
        print STDERR "[$period] $number taxonomy processed.\r";
    }
    
	chomp;
	my @fields = split /\t/, $_;
    next if $fields[1] == 0;
    
    my $taxID = $fields[4];
	my $name = $fields[5];
	$name =~ s/^ +//;
    
    my $rank = getTaxRank($taxID) || 'no rank';
    my $depth = getTaxDepth($taxID) || '10';
    
    $rank = "unclassified" if $taxID == 0;
    $rank = "root" if $taxID == 1;
    $rank = "strain" if $rank eq "no rank" && $depth>6;
   
    #skip extra strain levels
    #no rank -- Escherichia coli O111:H-
    #no rank --  Escherichia coli O111:H- str. 11128
    next if $rank eq "strain" && $fields[2] == 0;

    $taxa->{$rank}->{$name}->{READ_COUNT} = $fields[2];
	$taxa->{$rank}->{$name}->{ROLLUP}   = $fields[1];
	$taxa->{$rank}->{$name}->{TAXID}      = $taxID;
}
$period = &timeInterval($time);
my $number = $count;
$number =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
print STDERR "[$period] $number sequences processed.\n";
 
print "LEVEL\tTAXA\tROLLUP\tASSIGNED\tTAXID\n";
foreach my $rank ( sort {$level{$a}<=>$level{$b}} keys %$taxa )
{
	foreach my $name ( sort { $taxa->{$rank}->{$b}->{ROLLUP} <=> $taxa->{$rank}->{$a}->{ROLLUP} } keys %{$taxa->{$rank}}  ) {
		printf "%s\t%s\t%s\t%s\t%s\n",
			$rank,
			$name,
			$taxa->{$rank}->{$name}->{ROLLUP},
			defined $taxa->{$rank}->{$name}->{READ_COUNT} ? $taxa->{$rank}->{$name}->{READ_COUNT} : "",
			$taxa->{$rank}->{$name}->{TAXID}
        ;
	}
}

sub timeInterval{
    my $now = shift;
    $now = time - $now;
    return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}
