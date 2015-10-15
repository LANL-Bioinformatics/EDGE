#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use gi2lineage;

$|=1;

my $min_count;
$min_count = defined $ARGV[0] ? $ARGV[0] : 0;

# load taxonomy
loadTaxonomy();

my %major_level = (
	'root'         => 0,
	'superkingdom' => 10,
	'phylum'       => 20,
	'class'        => 30,
	'order'        => 40,
	'family'       => 50,
	'genus'        => 60,
	'species'      => 70,
	'strain'       => 80,
    'misc'         => 90
);

my $taxID = 0;
my $taxa;

while(<STDIN>)
{
	chomp;
	next if /^#/;
	next if /^$/;

    # taxaID: 28035
    #     Bacteria; Firmicutes; Bacilli; Bacillales; Staphylococcus.
    #     genus: Staphylococcus
    #     species: lugdunensis
    #     read count: 3
    #     reported gene count: 1
    
    if( /^taxaID: (\d+)$/ ){
        $taxID = $1;
    }
    next unless $taxID;

    if( $taxID > 1 && /^\tread count: (\d+)$/ ){
        my $read_count = $1;

        if( $read_count < $min_count ){
            $taxID=0;
            next;
        }

	    if( $taxID )
	    {
	    	my $rank = getTaxRank($taxID);

            #deal with non-ncbi taxid
	        unless($rank){
                $taxa->{"misc"}->{"non-ncbi taxid"}->{READ_COUNT} ||= 0;
	    	    $taxa->{"misc"}->{"non-ncbi taxid"}->{READ_COUNT} += $read_count;
                $taxa->{"misc"}->{"non-ncbi taxid"}->{ROLLUP} ||= 0;
	    	    $taxa->{"misc"}->{"non-ncbi taxid"}->{ROLLUP} += $read_count;
                next;
            }

            $rank = "strain" if $rank eq "no rank";
            my $name = getTaxName($taxID);

	    	my $p_id = $taxID;
	    	while( !defined $major_level{$rank} ){
	    		$p_id = getTaxParent($p_id);
	    		$rank = getTaxRank($p_id);
                $name = getTaxName($p_id);
	    	}
            
            $taxa->{$rank}->{$name}->{READ_COUNT} ||= 0;
	    	$taxa->{$rank}->{$name}->{READ_COUNT} += $read_count;
	    	
            while( $taxID )
	    	{
	    		$rank = getTaxRank($taxID);
	    		$rank = "strain" if $rank eq "no rank";

	    		$name = getTaxName($taxID);
                last if $name eq 'root';
	    
	    		if( defined $major_level{$rank} ){
	    			$taxa->{$rank}->{$name}->{ROLLUP} ||= 0;
	    			$taxa->{$rank}->{$name}->{ROLLUP} += $read_count;
	    	        $taxa->{$rank}->{$name}->{TAXID} = $taxID;
	    		}
	    		$taxID = getTaxParent($taxID);
	    	}
	    }
    
        $taxID=0;
    }
}

print "LEVEL\tTAXA\tROLLUP\tASSIGNED\tTAXID\n";
foreach my $rank ( sort {$major_level{$a}<=>$major_level{$b}} keys %$taxa )
{
	foreach my $name ( sort { $taxa->{$rank}->{$b}->{ROLLUP} <=> $taxa->{$rank}->{$a}->{ROLLUP} } keys %{$taxa->{$rank}}  ) {
		my $result = join "\t", @{$taxa->{$rank}->{$name}->{RESULT}} if defined $taxa->{$rank}->{$name}->{RESULT};
		printf "%s\t%s\t%s\t%s\t%s\n",
			$rank,
			$name,
			$taxa->{$rank}->{$name}->{ROLLUP},
			defined $taxa->{$rank}->{$name}->{READ_COUNT} ? $taxa->{$rank}->{$name}->{READ_COUNT} : "",
            $taxa->{$rank}->{$name}->{TAXID};
	}
}
