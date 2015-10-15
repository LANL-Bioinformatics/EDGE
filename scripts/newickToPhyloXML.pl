#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
## load Bioperl modules
use Bio::TreeIO;
use Bio::TreeIO::phyloxml;
use Bio::Tree::TreeI;
use Bio::Tree::NodeI;
use Bio::Tree::TreeFunctionsI;
## load Bio::Phylo modules
use lib "$RealBin/../lib";
use Bio::Phylo::IO qw (parse unparse);
use Bio::Phylo::Factory;
use Bio::Phylo::Util::CONSTANT qw':objecttypes :namespaces';

my $input;
my $annotation;
my $midpoint_flag;
my $outputDir="Output";
my $version = "0.1";

GetOptions(
   'i=s'      => \$input,
   'a=s'      => \$annotation,
   'o=s'      => \$outputDir,
   'm'        => \$midpoint_flag,
   'help|h'   => sub{usage()},
);

unless($input) {print "\nPlease provide newick tree file\n\n"; &usage(); exit;};

my ($input_name,$input_path,$input_suffix)=fileparse("$input",qr/\.[^.]*/);
my $output_newick = "$outputDir/$input_name$input_suffix.nwk";
my $output_phyloxml= "$outputDir/$input_name$input_suffix.xml";

# midpoint process using Bioperl Bio::TreeIO

my $treeio = Bio::TreeIO->new(-format => 'newick',
                              -file => "$input");
my $tree = $treeio->next_tree;

$tree = midpoint_reroot($tree) if ($midpoint_flag);

my $newick_out = Bio::TreeIO->new(-format => 'newick',
                                  -file   => ">$output_newick");
$newick_out->write_tree($tree);

# process using Bio::Phylo to add xml annotation

my $newick_out_string=$tree->as_text('newick');
my $proj = Bio::Phylo::IO->parse(
 '-file' => "$output_newick",
 '-format' => 'newick',
 '-as_project' =>1
);

my ($forest) = @{ $proj->get_items(_FOREST_) }; 
$tree = $forest->first; 
add_annotation($annotation,$tree) if ($annotation);

open (my $fh ,">$output_phyloxml") or die "Cannot write $output_phyloxml\n";
print  $fh unparse('-format' => 'phyloxml', '-phylo' => $proj);
close  $fh ;


sub add_annotation
{
    my $annotaion = shift;
    my $tree=shift;
    my $annotation_r=&read_annotation($annotation);
    my $fac = Bio::Phylo::Factory->new;
    foreach my $node ( @{ $tree->get_entities } ) {
        my $node_name=$node->get_name;
        my $desc = $annotation_r->{$node_name}->{desc};
        my $url = $annotation_r->{$node_name}->{url};
        if (defined $annotation_r->{$node_name}->{id})
        {
            my $arch = _create_dummy_architecture($desc,$url); 
            $node->add_meta( 
                $fac->create_meta( 
                        '-namespaces' => { 'pxml' => _NS_PHYLOXML_ }, 
                        '-triple' => { 'pxml:annotation' => $arch }, 
                ) 
            );
        }
    }
}


sub _create_dummy_architecture { 
    my $desc=shift;
    my $url=shift;
return 
"
	<desc>$desc</desc> 
	<uri>$url</uri> 
"; 
} 


sub add_annotation_not_used
{
    my $annotaion = shift;
    my $tree=shift;
    my $annotation_r=&read_annotation($annotation);
    my @names = keys %{$annotation_r};

    foreach my $name(@names)
    {
        my ($obj) = $tree->find_node("$name");  
        my $desc = $annotation_r->{$name}->{desc};
        my $url = $annotation_r->{$name}->{url};
        $treeio->add_phyloXML_annotation(
            -obj => $obj,  
            -xml => "
            <annotation> 
                  <desc>$desc</desc>
                  <url>$url</url>
            </annotation>    
                "
            );   
    }
}




sub midpoint_reroot 
{
    my $tree = shift;
    #get leaf
    my @leaves = $tree->get_leaf_nodes;

    #find midpoint
    my ($dist,$node1,$node2);
    $dist=0;

    for( my $i=0; $i<=$#leaves; $i++ ){
	    for( my $j=0; $j<=$#leaves; $j++ ){ 
		    my $curdist = $tree->distance($leaves[$i],$leaves[$j]);
		    if($curdist > $dist){
			    $dist = $curdist;
			    ($node1,$node2) = ($leaves[$i],$leaves[$j]);
		    }
	    }
    }

    my $midpoint = $dist/2;
    my $curdist = 0;
    my $node = $node1;
    while( $midpoint){
    	    my $anc = $node->ancestor;
            if (!$anc)
            { 
                $node=$node2;
                $anc = $node->ancestor;
            }
	    $curdist +=  $tree->distance($anc,$node);
                
	    last if $curdist > $midpoint;
	    $node = $anc;
    }

    $tree->reroot_at_midpoint($node);
    return $tree;
}

sub read_annotation
{
    my $file=shift;
    my %annotation;
    open (my $fh, "$file") or die "Cannot read $file\n";
    #Name Description URL
    my $header=<$fh>;
    while(<$fh>)
    {
        chomp;
        my ($id,$desc,$url)=split /\t/,$_;
        # Percent-encoding reserved characters
        $url =~ s/ /%20/g;
        $url =~ s/\&/%26/g;
        $url =~ s/\#/%23/g;
        $url =~ s/\(/%28/g;
        $url =~ s/\)/%29/g;
        $url =~ s/\"/%22/g;
        $annotation{$id}->{desc}=$desc;
        $annotation{$id}->{url}=$url;
        $annotation{$id}->{id}=$id;
    }
    close $fh;
    return (\%annotation);
}

sub usage
{
print <<"END";
     Usage: perl $0 [options] -i input.tre -o out_directory 
     Version $version
     Input:
            -i            Newick tree file
            
     Output:
            -o            Output Directory (default: Output)
                          will output .nwk and .xml
     
     Options:
            -m            Reroot with Midpoint
 
            -a            annotation file
                          A three tab-delimited columns table with a header
                          Name  Description  URL
            -version      Print verison

END
exit;
}
