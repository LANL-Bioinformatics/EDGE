#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($mplnfile, $level);

GetOptions(
	'i=s', \$mplnfile,
	'l=s', \$level
);

my %lvl_name = (
	'k' => 'kingdom',
	'p' => 'phylum',
	'c' => 'class',
	'o' => 'order',
	'f' => 'family',
	'g' => 'genus',
	's' => 'species'
);

my @levels = ( 'k', 'p', 'c', 'o', 'f', 'g', 's');

my $taxa;
my $count = 0;

open MPLN, "$mplnfile" || die "Can't open $mplnfile: $!\n";
while(<MPLN>){
	chomp;
	my ($phylo,$count) = $_ =~ /^(\S+)\t(\S+)$/;
	my ($lvl,$name) = $phylo =~ /(\w)__([^\|]+)$/;
	$name =~ s/_/ /g;
	$taxa->{$lvl}->{$name} = $count;
}
close MPLN;

foreach my $lvl ( @levels )
{
	foreach my $name ( sort { $taxa->{$lvl}->{$b} <=> $taxa->{$lvl}->{$a} } keys %{$taxa->{$lvl}} )
	{
		#printf "%s\t%s\t%s\n",
		#	$lvl_name{$lvl},
		#	$name,
		#	$taxa->{$lvl}->{$name};
		my $value = $taxa->{$lvl}->{$name};
		$value = sprintf("%.0f",$value) ? sprintf("%.0f",$value) : 1;

		if ( $name =~ s/ unclassified// ){
			printf "%s\t%d\n", $name, $value;
			next;
		}
		printf "%s\t%d\n", $name, $value if $lvl eq $level;
	}
	last if $lvl eq $level;
}
