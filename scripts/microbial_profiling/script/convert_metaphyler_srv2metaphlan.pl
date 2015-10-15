#!/usr/bin/perl -w
use Storable;
use strict;

my ($infile, $cutoff) = @ARGV;
my $node_phylo;
my $phylo_all;
$cutoff = 0.9 unless $cutoff;

if( !$infile ){
	print STDERR "USAGE: $0 <metaphyler_output> [<cutoff>]\n";
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

my @lvl = ('g','f','o','c','p');

print STDERR "Parsing metaphylerSRV output file..." ;
open NODES, $infile or die "$!\n";
while (<NODES>) {
	chomp;
	next if /@/;
	#"prefix".classify.tab
	#column 0: query sequence id
	#column 2: phylogenetic marker gene name
    #    column 2: best reference gene hit
    #    column 3: % similarity with best hit
    #    column 4: classification rule
	#column 5-9: taxonomic label at genus,family,order,class,phylum level

    #"prefix".genus|family|order|class|phylum.tab
    #    column 1: taxonomic clade name
    #    column 2: % relative abundances
    #    column 3: depth of coverage of genomes
    #    column 4: number of sequences binned to this clade
    #    column 5: similarity with reference genes (only available at the genus level)	

	my @temp = split /\t/, $_;
	my @taxas = @temp[5..9];
	my $name = $taxas[0];
	my @ranks;
	
	push @ranks, "k__Bacteria";
	for (my $i=4; $i>=0; $i--){
		last if $taxas[$i] =~ /\{\w+\}/;
		$taxas[$i] =~ s/ \(\w+\)//g;
		$taxas[$i] =~ s/ /_/g;
		push @ranks, "$lvl[$i]__$taxas[$i]";
	}

	#initial ref
	my $ref = \%{$phylo_all};
	foreach my $rank (@ranks) {
		my $idx = $1 if $rank =~ /^(\w)__/;
		last if $idx eq 's';

		die "\nERROR: Can't parse taxon - $name\n" unless $idx;
		$node_phylo->{$name}->{$idx}=$rank;
		$ref->{$idx}->{$rank}->{COUNT}=0;

		#creat an unclassified down-level
		my $next_idx = $dn_lvl{$idx};
		my $next_rank = $rank;
		$next_rank =~ s/^$idx\__/$next_idx\__/;
		$next_rank .= "_unclassified";
		$ref->{$idx}->{$rank}->{$next_idx}->{$next_rank}->{COUNT}=0;

		#reference
		$ref = \%{$ref->{$idx}->{$rank}};
	}
}
close NODES;
print STDERR "Done\n";

print STDERR "Parsing metaphyler classification file...";
open MPYL, $infile || die "$!\n";
while (<MPYL>) {
	chomp;
	next if $_ =~ /^#/;
	next if $_ =~ /^\s*$/;
	next if /@/;

	my @temp = split "\t", $_;

	my ($tid, $score) = ($temp[5],1);
	foreach my $tid_scr ( @temp[5..9] ){
		next if $tid_scr =~ /\{\w+\}/;

		my ($nid, $assigned_num) = ($tid, $score);
	
		my $ref = \%{$phylo_all};
		foreach my $idx (('k','p','c','o','f','g')) {
			last unless $node_phylo->{$nid}->{$idx};
	
			my $rank = $node_phylo->{$nid}->{$idx};
			$ref->{$idx}->{$rank}->{COUNT} += $assigned_num;
			
			my $next_idx = $dn_lvl{$idx};
			unless ( defined $node_phylo->{$nid}->{$next_idx} ) {
				my $next_rank = $rank;
				$next_rank =~ s/^$idx\__/$next_idx\__/; 
				$next_rank .= "_unclassified";
				$ref->{$idx}->{$rank}->{$next_idx}->{$next_rank}->{COUNT} += $assigned_num;
			}
			
			$ref = \%{$ref->{$idx}->{$rank}};
		}
		last;
	}
}
close MPYL;
print STDERR "Done\n";

print STDERR "Output results...\n";
my $tol=0;
foreach my $idx ( keys %{$phylo_all->{k}} ){
	$tol += $phylo_all->{k}->{$idx}->{COUNT};
}

&travel( \%{$phylo_all}, "k", "");
print STDERR "Done\n";
print STDERR "Total read assigned: $tol\n";

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
		#print $out, "\t", ($count/$tol)*100, "\n" if $count > 0;
		print $out, "\t", $count, "\n" if $count > 0;
	}
}

