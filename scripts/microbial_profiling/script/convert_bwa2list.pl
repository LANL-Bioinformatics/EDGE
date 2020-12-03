#! /usr/bin/perl -w
use strict;
#use Storable;
use lib (`ktGetLibPath`);
use KronaTools;

$|=1;

# load taxonomy
print STDERR "Loading taxonomy...\n";
loadTaxonomy();

# parse BLAST results
my $fileName=$ARGV[0];

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

open TAXA, $fileName or die "Can't open $fileName\n";
while(<TAXA>)
{
	chomp;
	next if /^#/;
	next if /^$/;

	#0 Accession
	#1 Length
	#2 GC%
	#3 Avg_fold
	#4 Fold_std
	#5 Base_Coverage%
	#6 Mapped_reads
	#7 Linear_length
	#8 ID

	my @fields = split /\t/, $_;

	if ( $fields[0] eq "GI" ){
		@headers = @fields[0..7];
		next;
	}

	my $acc = $fields[0];
	my $taxID = getTaxIDFromGI($acc);

	print STDERR "[WARNING] Can't find Accession#$acc: $fields[8]\n" unless $taxID;

	my $rank = "replicon";
	my $name = $fields[8];

	$taxa->{$rank}->{$name}->{READ_COUNT} = 0 unless defined $taxa->{$rank}->{$name}->{READ_COUNT};
	$taxa->{$rank}->{$name}->{READ_COUNT} += $fields[6];
	$taxa->{$rank}->{$name}->{ROLLUP} = 0 unless defined $taxa->{$rank}->{$name}->{ROLLUP};
	$taxa->{$rank}->{$name}->{ROLLUP} += $fields[6];
	@{$taxa->{$rank}->{$name}->{RESULT}} = @fields[0..7];

	if( $taxID )
	{
		#print STDERR "Processing GI:$gi TAXID:$taxID NAME:'$fields[8]'...\n";

		$rank = getTaxRank($taxID);
		$rank = "strain" if $rank eq "no rank";
		$name = getTaxName($taxID);

		my $p_id = $taxID;
		while( !defined $major_level{$rank} ){
			$p_id = getTaxParent($p_id);
			$rank = getTaxRank($p_id);
			$name = getTaxName($p_id);
		}
		$taxa->{$rank}->{$name}->{READ_COUNT} = 0 unless defined $taxa->{$rank}->{$name};
		$taxa->{$rank}->{$name}->{READ_COUNT} += $fields[6];

		while( $taxID )
		{
			$rank = getTaxRank($taxID);
			$rank = "strain" if $rank eq "no rank";

			$name = getTaxName($taxID);
			last if $name eq 'root';
	
			if( defined $major_level{$rank} ){
				$taxa->{$rank}->{$name}->{ROLLUP} = 0 unless defined $taxa->{$rank}->{$name}->{ROLLUP};
				$taxa->{$rank}->{$name}->{ROLLUP} += $fields[6];
			}
			$taxID = getTaxParent($taxID);
		}
	}
}
close TAXA;
print STDERR "\n";

my $header = join "\t", @headers;
print "LEVEL\tTAXA\tROLLUP\tASSIGNED\t$header\n";
foreach my $rank ( sort {$major_level{$a}<=>$major_level{$b}} keys %$taxa )
{
	foreach my $name ( sort { $taxa->{$rank}->{$b}->{ROLLUP} <=> $taxa->{$rank}->{$a}->{ROLLUP} } keys %{$taxa->{$rank}}  ) {
		my $result = join "\t", @{$taxa->{$rank}->{$name}->{RESULT}} if defined $taxa->{$rank}->{$name}->{RESULT};
		printf "%s\t%s\t%s\t%s\t%s\n",
			$rank,
			$name,
			$taxa->{$rank}->{$name}->{ROLLUP},
			defined $taxa->{$rank}->{$name}->{READ_COUNT} ? $taxa->{$rank}->{$name}->{READ_COUNT} : "",
			defined $result ? $result : "";
	}
}
