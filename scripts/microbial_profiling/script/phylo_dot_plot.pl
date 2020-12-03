#!/usr/bin/perl -w
use Bio::TreeIO;
use Bio::Tree::TreeI;
use Bio::Tree::Node;
use Getopt::Long;
use IO::String;
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
            'score=s',
            'scale=s',
            'mode=s',
            'title=s',
            'help|?') || &usage();
&usage() if !-r $opt{input} || !defined $opt{prefix};

sub usage {
	print "
USAGE: $0 --input <TREE_TAB_FILE> --prefix <OUTPUT_PREFIX> [OPTIONS]
    
    --input  | -i <STRING>    A text file with read counts and taxonomic lineage using
                              TAB separator. The first column has to be read counts
                              following by the lineage from upper-rank (broader) to more
                              specific ranks.

    --score  | -l <STRING>    tsv with '[species name]<tab>[score]'
    --prefix | -p <STRING>    file prefix for output

OPTIONS:
    --scale | -s <STRING>     Scaling factor: max or total. Default is \"max\".
    --mode  | -m <STRING>     Scaling mode: normal or log. Default is \"log\".
    --title | -t <STRING>     Title
";
exit;
}

my $intab  = $opt{input};
my $prefix = $opt{prefix};
my $scale  = $opt{scale};
my $mode   = $opt{mode};
my $title  = $opt{title};
$scale   ||= "max";
$mode    ||= "log";

my %taxa;
my @cons_tree;
my %score_hash;
my $score_provided = 0;
$score_provided = 1 if -e $opt{score};

die("[phylo_dot_plot] ERROR: Input file is empty.\n") unless( -s $opt{input} );

#parsing tab file
open TAB, "$opt{input}" or die "Can't open input file: $!\n";
while(<TAB>){
	chomp;
	my @temp = split /\t/, $_;
	#keep lineage
	push @cons_tree, \@temp;
}
close TAB;

if($score_provided){
	open TAB, "$opt{score}" or die "Can't open score file: $!\n";
	while(<TAB>){
		chomp;
		my @temp = split /\t/, $_;
		$score_hash{$temp[0]} = $temp[1];
		$score_hash{$temp[0]} = "NA" if $temp[1] eq "";
	}
	close TAB;
}

#init a tree
my $init_tree = IO::String->new(";");
my $intre = Bio::TreeIO->new( -fh=> $init_tree, -format => 'newick' );
my $tree = $intre->next_tree();

#adding and taging tree nodes
foreach my $ref ( @cons_tree ){
	my @taxas = @$ref;
	my $cnt = shift @taxas;
	my $root = $tree->get_root_node;

	my $cur_node = $root;
    my $depth = 0;

	foreach my $taxa ( @taxas ){
		my @child_nodes = $cur_node->each_Descendent;
		$cur_node->add_tag_value("depth", $depth) unless $cur_node->has_tag("depth");
		my $taxa_node;

		if( scalar @child_nodes ){
			my @nodes = grep { defined $_->id && $_->id eq $taxa } @child_nodes;
			$taxa_node = shift @nodes;
		}
		
		unless($taxa_node){
			$taxa_node = new Bio::Tree::Node( -id => $taxa );
            $cur_node->add_Descendent($taxa_node);
		}

		$cur_node = $taxa_node;
		$depth++;
	}
	my $taxa = $cur_node->id;
	my $score = "NA";
	$score = $score_hash{$taxa} if defined $score_hash{$taxa};
	$cur_node->add_tag_value("count", $cnt);
	$cur_node->add_tag_value("score", $score);
}

#draw tree
my $svg = new phylo_dot_plot(
	-tree    => $tree,
	-type    => 'PHYLO',
	-compact => 'compact',
	-branch  => 170,
	-scale   => $scale,
	-mode    => $mode,
	-size    => 16,
	-title   => $title,
);

open FILE, ">$prefix.svg";
my $content = $svg->print(-caveh => '100%', -cavew => '100%', -format => 'svg');
print FILE $content;
close FILE;
