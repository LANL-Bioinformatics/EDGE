#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage; 

$|=1;

# load taxonomy
print STDERR "Loading taxonomy...\n";
loadTaxonomy();

my %major_level = (
	'superkingdom' => 10,
	'phylum'       => 20,
	'class'        => 30,
	'order'        => 40,
	'family'       => 50,
	'genus'        => 60,
	'species'      => 70,
	'strain'       => 80,
	'replicon'     => 90
);

my @headers;
my $taxa;

while(<STDIN>)
{
	chomp;
	next if /^#/;
	next if /^$/;

	my @fields = split /\t/, $_;

	next unless $fields[0];

	my $acc = $fields[0];
	my $taxID = getTaxIDFromAcc($acc);

	print STDERR "[WARNING] Can't find Accession#$acc\n" unless $taxID;

	my $rank = "replicon";
	my $name = $acc;

	$taxa->{$rank}->{$name}->{READ_COUNT} = 0 unless defined $taxa->{$rank}->{$name}->{READ_COUNT};
	$taxa->{$rank}->{$name}->{READ_COUNT} += $fields[1];
	$taxa->{$rank}->{$name}->{ROLLUP} = 0 unless defined $taxa->{$rank}->{$name}->{ROLLUP};
	$taxa->{$rank}->{$name}->{ROLLUP} += $fields[1];

	if( $taxID )
	{
		#print STDERR "Processing GI:$gi TAXID:$taxID NAME:'$fields[8]'...\n";

		$rank = getTaxRank($taxID);
		$rank ||= "no rank";
		$rank = "strain" if $rank eq "no rank";
		$name = getTaxName($taxID);

		my $p_id = $taxID;
		while( !defined $major_level{$rank} ){
			$p_id = getTaxParent($p_id);
			$rank = getTaxRank($p_id);
			$name = getTaxName($p_id);
		}
		$taxa->{$rank}->{$name}->{READ_COUNT} = 0 unless defined $taxa->{$rank}->{$name};
		$taxa->{$rank}->{$name}->{READ_COUNT} += $fields[1];

		while( $taxID )
		{
			$rank = getTaxRank($taxID);
			$rank = "strain" if $rank eq "no rank";

			$name = getTaxName($taxID);
			last if $name eq 'root';
	
			if( defined $major_level{$rank} ){
				$taxa->{$rank}->{$name}->{ROLLUP} = 0 unless defined $taxa->{$rank}->{$name}->{ROLLUP};
				$taxa->{$rank}->{$name}->{ROLLUP} += $fields[1];
			}
			$taxID = getTaxParent($taxID);
		}
	}
}

print STDERR "\n";
print "LEVEL\tTAXA\tROLLUP\tASSIGNED\n";
foreach my $rank ( sort {$major_level{$a}<=>$major_level{$b}} keys %$taxa )
{
	foreach my $name ( sort { $taxa->{$rank}->{$b}->{ROLLUP} <=> $taxa->{$rank}->{$a}->{ROLLUP} } keys %{$taxa->{$rank}}  ) {
		printf "%s\t%s\t%s\t%s\n",
			$rank,
			$name,
			$taxa->{$rank}->{$name}->{ROLLUP},
			defined $taxa->{$rank}->{$name}->{READ_COUNT} ? $taxa->{$rank}->{$name}->{READ_COUNT} : "";
	}
}
