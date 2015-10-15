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

open MPLN, $mplnfile || die "Can't open $mplnfile: $!\n";
while(<MPLN>){
	chomp;
	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;

	my ($lvl,$name) = $phylo =~ /(\w)__([^\|]+)$/;

	if( $name =~ s/_unclassified// ){
		$taxa->{$up_lvl{$lvl}}->{$name}->{ASSIGNED} = $count;
	}
	else {
		$taxa->{$lvl}->{$name}->{ROLLUP} = $count;
		$taxa->{$lvl}->{$name}->{ASSIGNED} = $count if $lvl eq $level;
	}
}
close MPLN;

print "LEVEL\tTAXA\tROLLUP\tASSIGNED\n";
foreach my $lvl ( @levels )
{
	if($lvl eq $level){
		print "";
	}
	foreach my $name ( sort { $taxa->{$lvl}->{$b}->{ROLLUP} <=> $taxa->{$lvl}->{$a}->{ROLLUP} } keys %{$taxa->{$lvl}} )
	{
		my $temp = $name;
		$temp =~ s/_/ /g;
		printf "%s\t%s\t%s\t%s\n",
			$lvl_name{$lvl},
			$temp,
			$taxa->{$lvl}->{$name}->{ROLLUP},
			defined $taxa->{$lvl}->{$name}->{ASSIGNED} ? $taxa->{$lvl}->{$name}->{ASSIGNED} : "";
	}
	last if $lvl eq $level;
}
