#!/usr/bin/env perl
use strict;
use File::Find;
use File::Basename;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use JSON;

my $basename_opt;
my $sort_by_size_opt;
GetOptions(
	'-basename' => \$basename_opt,
	'-sort_by_size'	=>\$sort_by_size_opt
);

my $list={};
my @dir = @ARGV;
usage() if scalar @dir < 1;

find(\&parseRefGb, @dir);
find(\&parseRefFa, @dir);
$list = sort_by_size($list) if ($sort_by_size_opt);
print to_json($list,{pretty=>1});
#print to_json($list);

sub parseRefGb {
	($basename_opt)?
	$File::Find::name =~ /([^\/]+)\/([^\/]+\.(gb|gbk|gbff|genbank))$/i && push @{$list->{$1}}, $_:
	$File::Find::name =~ /([^\/]+)\/([^\/]+\.(gb|gbk|gbff|genbank))$/i && push @{$list->{$1}}, $File::Find::name;
}

sub parseRefFa {
	($basename_opt)?
	$File::Find::name =~ /([^\/]+)\/([^\/]+\.(fa|fna|fasta))$/i && push @{$list->{$1}}, $_:
	$File::Find::name =~ /([^\/]+)\/([^\/]+\.(fa|fna|fasta))$/i && push @{$list->{$1}}, $File::Find::name;
}

sub sort_by_size {
        my $list = shift;
	my $new_list={};
   	foreach my $id  (keys %{$list})
	{
		my @tmp;
		if ($basename_opt)
		{
			map { my( $name,$dir,$sufix )=fileparse($_, qr/\.[^.]*/); push @tmp, $name;} sort { -s $b <=> -s $a } @{$list->{$id}};
		}
		else
		{
			@tmp = sort { -s $b <=> -s $a } @{$list->{$id}};
		}
		$id =~ s/uid\d+//;
		$id =~ s/\_+$//;
		@{$new_list->{$id}} = @tmp;
	}
	return $new_list;
}

sub usage {
	print STDERR "
USAGE: $0 <dir1>[,<dir2>,<dir3>...] > \$EDGE_HOME/edge_ui/data/ref_list.json
	-basename 	only print filename, no path
	-sort_by_size	sort file by its size (Descending)

";
	exit;
}
