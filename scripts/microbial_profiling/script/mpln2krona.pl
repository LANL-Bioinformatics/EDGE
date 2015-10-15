#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($mplnfile, $level);

GetOptions(
	'i=s', \$mplnfile,
	'l=s', \$level
);

my @levels = ( 'k', 'p', 'c', 'o', 'f', 'g', 's');

my %lvl_name = (
	'k' => 'kingdom',
	'p' => 'phylum',
	'c' => 'class',
	'o' => 'order',
	'f' => 'family',
	'g' => 'genus',
	's' => 'species'
);

my %up_lvl = (
	'p' => 'k',
	'c' => 'p',
	'o' => 'c',
	'f' => 'o',
	'g' => 'f',
	's' => 'g'
);

my $taxa;
my $count = 0;

open MPLN, $mplnfile || die "Can't open $mplnfile: $!\n";
while(<MPLN>){
	chomp;
	next if /LEVEL/;
	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;
	my ($lvl,$name) = $phylo =~ /(\w)__([^\|]+)$/;

	if( $name =~ s/_unclassified// ){
		$taxa->{$up_lvl{$lvl}}->{$name}->{ASSIGNED} = $count;
	}
	else {
		$taxa->{$lvl}->{$name}->{ROLLUP} = $count;
		$taxa->{$lvl}->{$name}->{ASSIGNED} = $count if $lvl eq $level;
		$taxa->{$lvl}->{$name}->{LINEAGE}  = $phylo;
	}
}
close MPLN;

foreach my $lvl ( @levels )
{
	foreach my $name ( keys %{$taxa->{$lvl}} )
	{
		if ( $taxa->{$lvl}->{$name}->{ASSIGNED} > 0 ){
			my $assigned = $taxa->{$lvl}->{$name}->{ASSIGNED};
			my $lineage = $taxa->{$lvl}->{$name}->{LINEAGE};
			$lineage =~ s/\|/\t/g;
			$lineage =~ s/\w__//g;
			$lineage =~ s/_/ /g;
			print "$assigned\t$lineage\n";
		}
	}
}
