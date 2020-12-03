#! /usr/bin/perl

# Copyright © 2011, Battelle National Biodefense Institute (BNBI);
# all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and
# Adam Phillippy
#
# See the LICENSE.txt file included with this software for license information.


use strict;
use Storable;

use lib (`ktGetLibPath`);
use KronaTools;

setOption('out', 'gottcha.krona.html');
setOption('name', 'Root');

my @options =
qw(
	out
	name
	include
	combine
	queryCol
	giCol
	scoreCol
	magCol
	depth
	hueBad
	hueGood
	local
	standalone
	url
);

getKronaOptions(@options);

my ($gCol, $qCol, $sCol, $mCol);
$gCol = 1 unless defined getOption('giCol');
$qCol = 1;  #unless defined getOption('queryCol');
$sCol = 5;  #unless defined getOption('scoreCol');
$mCol = 7 unless defined getOption('magCol');

if
(
	@ARGV < 1
)
{
	printUsage
	(
		'Creates a Krona chart based on taxonomy IDs and, optionally, magnitudes
and scores. Taxonomy IDs corresponding to a rank of "no rank" in the database
will be assigned to their parents to make the hierarchy less cluttered (e.g.
"Cellular organisms" will be assigned to "root").',
		'taxonomy',
		'Tab-delimited file with taxonomy IDs and (optionally) query IDs,
magnitudes and scores. By default, query IDs, taxonomy IDs and scores will be
taken from columns 1, 2 and 3, respectively (see -q, -t, -s, and -m). Lines
beginning with "#" will be ignored.',
		1,
		1,
		\@options
	);
	
	exit 0;
}

if ( optionsConflict('queryCol', 'taxCol', 'magCol', 'scoreCol') )
{
	ktWarn('Query column already in use; not reading query IDs.');
	setOption('queryCol', undef);
}

if ( optionsConflict('scoreCol', 'taxCol', 'magCol') )
{
	ktWarn('Score column already in use; not reading scores.');
	setOption('scoreCol', undef);
}

if ( optionsConflict('magCol', 'taxCol') )
{
	ktWarn('Magnitude column already in use; not reading magnitudes.');
	setOption('magCol', undef);
}

my $tree = newTree();

print "Loading taxonomy...\n";
loadTaxonomy();
#print "Loading species tree...\n";
#my $speciesTree = retrieve "/lato/traceyf/db/custom/2013-03/speciesTreeGI.dmp";
my $speciesTree;
#print "Loading genome tree...\n";
#my $genomeVitals = retrieve "/lato/traceyf/db/custom/2013-03/genomeVitals.dmp";
my $genomeVitals;

my $set = 0;
my @datasetNames;
my $useScore = 1; # is score column present?
my $eVal; # is score e-value (inferred by negative scores)?
my $useMag = getOption('magCol'); # using magnitude values?
my @header;

#my $scoreLineage = &parseScore();

foreach my $input (@ARGV)
{
	my ($file, $magFile, $name) = parseDataset($input);
	
	if ( !getOption('combine') ){
		push @datasetNames, $name;
	}
	
	print "Importing $file...\n";
	
	open IN, "<$file" or ktDie("Couldn't open \"$file\".");

	while ( <IN> )
	{
		next if ( /^#/ );
		chomp;
		if ( /^Accession/ ){
			@header = split /\t/, $_ ;
			next;
		}
		
		my @fields = split /\t/;
		my $acc = $fields[0];
		my $queryID;

		my $magnitude;
		
		$queryID   = $header[$mCol-1] if ($mCol);
		$magnitude = $fields[$mCol-1];
		
		#addByLineage($tree, $set, &getLineage($gi), undef, $magnitude);
		#my $gi = convName2Gi($name);
		
		my $taxID = getTaxIDFromAcc($acc);
		
		my @score;
		
		if ( $taxID ) {
			addByTaxID($tree, $set, $taxID, $queryID, $magnitude);
		}
		#else{
	#		my @lineage = &getLineage($gi);
	#		unless ( scalar @lineage ){
	#			print "NOT FOUND: $gi\t$fields[8]\n";
	#			next;
	#		}
	#		addByLineage($tree, $set, \@lineage, $queryID, $magnitude);
	#	}
	}

	if ( ! getOption('combine') )
	{
		$set++;
	}
	
	close IN;
}

sub getLineage {
	my ( $gi ) = @_;
	my $taxid;
	my $strain;
	my @lineage;

	#find lineage in species tree
	SPE: foreach my $tid (keys %$speciesTree) {
		$taxid = $tid;
		if( exists $speciesTree->{$taxid}->{GI} ) {
			foreach my $rep ( keys %{$speciesTree->{$taxid}->{GI}} ){
				foreach my $db ( keys %{$speciesTree->{$taxid}->{GI}->{$rep}} ){
					if( exists $speciesTree->{$taxid}->{GI}->{$rep}->{$db}->{$gi} ){
						$strain = $speciesTree->{$taxid}->{GI}->{$rep}->{$db}->{$gi};
						last SPE;
					}
				}
			}
    		}
	}

	#if found, add lineage into an array
	if ( defined $strain ) {
		my @ranks = ( 'SK', 'P', 'C', 'O', 'F', 'G', 'S' );
		foreach my $rank ( @ranks ){
			my $name = $speciesTree->{$taxid}->{$rank};
			if( $rank eq 'S' ){
				foreach my $species ( keys %{$speciesTree->{$taxid}->{$rank}} ) {
					$name = $species if $speciesTree->{$taxid}->{$rank}->{$species} eq 'scientific name';
				}
			}
			push @lineage, $name;
		}
		push @lineage, $strain;
	}

	return @lineage;
}

sub optionsConflict {
	my ($option, @others) = @_;
	
	if ( getOption($option) )
	{
		foreach my $other ( @others )
		{
			if ( getOption($option) == getOption($other) )
			{
				return 1;
			}
		}
	}
	
	return 0;
}

my @attributeNames =
(
        'magnitude',
);

my @attributeDisplayNames =
(
        $header[$mCol-1],
);

my @scoreArgs;

if ( $useScore )
{
	@scoreArgs =(
		$eVal ? getOption('hueGood') : getOption('hueBad'),
		$eVal ? getOption('hueBad') : getOption('hueGood')
	)
}

writeTree
(
	$tree,
	\@attributeNames,
	\@attributeDisplayNames,
	\@datasetNames
);

