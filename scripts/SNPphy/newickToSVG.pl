#!/usr/bin/env perl
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::Node;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";
use phylo_dot_plot;
use strict;


$ENV{PATH} = "$Bin:$Bin/../:$Bin/script/:$Bin/bin/:$ENV{PATH}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
            'input=s',
            'prefix=s',
            'scale=s',
            'mode=s',
            'help|?') || &usage();
&usage() if !-r $opt{input} || !defined $opt{prefix};

sub usage {
        print "
USAGE: $0 --input <TREE_FILE> --prefix <OUTPUT_PREFIX> [OPTIONS]
    
    --input  | -i <STRING>    A newick tree file
    
    --prefix | -p <STRING>    file prefix for output

OPTIONS:
    --scale | -s <STRING>     Scaling factor: max or total. Default is \"max\".
    --mode  | -m <STRING>     Scaling mode: normal or log. Default is \"log\".
";
exit;
}

my $intree  = $opt{input};
my $prefix = $opt{prefix};
my $scale  = $opt{scale};
my $mode   = $opt{mode};
$scale   ||= "max";
$mode    ||= "log";

die("[phylo_dot_plot] ERROR: Input file is empty.\n") unless( -s $opt{input} );

#init a tree
my $intre = Bio::TreeIO->new( -file=> $intree, -format => 'newick' );
my $tree = $intre->next_tree();

#draw tree
my $svg = new phylo_dot_plot(
        -tree    => $tree,
        -type    => 'PHYLO',
        -compact => 'compact',
        -branch  => 170,
        -scale   => $scale,
        -mode    => $mode,
        -size    => 16
);

open FILE, ">$prefix.svg";
my $content = $svg->print(-caveh => '100%', -cavew => '100%', -format => 'svg');
print FILE $content;
close FILE;
