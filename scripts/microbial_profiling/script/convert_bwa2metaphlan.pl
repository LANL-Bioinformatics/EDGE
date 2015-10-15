#!/usr/bin/perl -w
use Storable;
use strict;

$|=1;

my ($sqd_in) = @ARGV;

my $phylo_all;
my $node_phylo;

if( !$sqd_in){
	print STDERR "USAGE: $0 <bwa_idmapping_output>\n";
	exit 1;
}

my %dn_lvl = (
	k => 'p',
	p => 'c',
	c => 'o',
	o => 'f',
	f => 'g',
	g => 's'
);

print STDERR "Loading species tree...\n";
my $speciesTree = retrieve "/lato/traceyf/db/custom/2013-03/speciesTreeGI.dmp";
#print STDERR "Loading genome tree...\n";
#my $genomeVitals = retrieve "/lato/traceyf/db/custom/2013-03/genomeVitals.dmp";

my $count=1;

print STDERR "Parsing node2phylo output file...\n" ;
open NODES, $sqd_in || die "$!\n";
while (<NODES>) {
	chomp;
	next if /^GI/;
	my @temp = split "\t", $_;
	my $line = &getLineage($temp[0]);

	my $name = $temp[0];
	my $phylo = $line;
	my @ranks = split /\|/, $phylo;

	#initial ref
	my $ref = \%{$phylo_all};
	foreach my $rank (@ranks) {
		my $idx = $1 if $rank =~ /^(\w)__/;
		#last if $idx eq 's';

		die "\nERROR: Can't parse taxon - $name\n" unless $idx;
		$node_phylo->{$name}->{$idx}=$rank;
		$ref->{$idx}->{$rank}->{COUNT}=0;


		last if $idx eq 's';
		#creat an unclassified down-level
		my $next_idx = $dn_lvl{$idx};
		my $next_rank = $rank;
		$next_rank =~ s/^$idx\__/$next_idx\__/;
		$next_rank .= "_unclassified";
		$ref->{$idx}->{$rank}->{$next_idx}->{$next_rank}->{COUNT}=0;

		#reference
		$ref = \%{$ref->{$idx}->{$rank}};
	}
	print STDERR "$count\r";
	$count++;
}
close NODES;
print STDERR "Done\n";

print STDERR "Parsing bwa output to tree...";
open MCV, $sqd_in || die "$!\n";
while (<MCV>) {
	chomp;
	next if $_ =~ /^#/;
	next if $_ =~ /ID/;
	next if $_ =~ /^\s*$/;

	my @temp = split "\t", $_;
	my ($nid, $assigned_num) = ($temp[0], $temp[6]);

	my $ref = \%{$phylo_all};
	foreach my $idx (('k','p','c','o','f','g','s')) {
		last unless $node_phylo->{$nid}->{$idx};

		my $rank = $node_phylo->{$nid}->{$idx};
		$ref->{$idx}->{$rank}->{COUNT} += $assigned_num;
		
		my $next_idx = $dn_lvl{$idx};
		if ( defined $next_idx && !defined $node_phylo->{$nid}->{$next_idx}) {
			my $next_rank = $rank;
			$next_rank =~ s/^$idx\__/$next_idx\__/; 
			$next_rank .= "_unclassified";
			$ref->{$idx}->{$rank}->{$next_idx}->{$next_rank}->{COUNT} += $assigned_num;
		}
		
		$ref = \%{$ref->{$idx}->{$rank}};
	}
}
close MCV;
print STDERR "Done\n";

print STDERR "Output results...";
my $tol=0;
foreach my $idx ( keys %{$phylo_all->{k}} ){
	$tol += $phylo_all->{k}->{$idx}->{COUNT};
}

&travel( \%{$phylo_all}, "k", "");
print STDERR "Done\n";
print STDERR "Total reads assigned: $tol\n";

sub travel {
	my ($ref, $idx, $str) = @_;
	foreach my $rank ( keys %{$ref->{$idx}} ){
		my $count = $ref->{$idx}->{$rank}->{COUNT};
		if( $count > 0 && defined $dn_lvl{$idx} ){
			my $next_idx = $dn_lvl{$idx};
			my $next = \%{$ref->{$idx}->{$rank}};
			&travel( $next, $next_idx, "$str|$rank" );
		}
		my $out = "$str|$rank";
		$out =~ s/^\|//;
		#print $out, "\t", ($count/$tol)*100, "\n" if $count > 0 && $type eq 'pcnt';
		print $out, "\t", $count           , "\n" if $count > 0;
	}
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
			
			$rank =~ s/SK/K/;
			$name =~ s/ /_/g;
			push @lineage, lc($rank)."__$name";
		}
	}

	my $out = join "\|",  @lineage;
	return $out;
}

