# phylo_dot_plot.pm
# ver 0.1
# 2013/12/16
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
#
# This BioPerl module is used to generate phylogenetic tree with 
# scaled dot (MEGAN-like) plot in SVG format.
#
# Some of the codes were taken and modified from
# Bio::Tree::Draw::Cladogram
# Cladogram.pm,v 1.8 2005/09/04 07:35:05 valiente Exp
# Cared for by Gabriel Valiente <valiente@lsi.upc.edu>
# Copyright Gabriel Valiente
# You may distribute this module under the same terms as Perl itself

package phylo_dot_plot;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

# The following private package variables are set by the new method
# and used by the print method.

my %xx;           # horizontal coordinate for each node
my %yy;           # vertical coordinate for each node
my $t1;           # first Bio::Tree::Tree object
my $font;         # font name
my $size;         # font size
my $width;        # total drawing width
my $height;       # total drawing height
my $xstep;        # branch length in drawing
my $tip;          # extra space between tip and label
my $tipwidth1;    # width of longest label among $t1 taxa
my $compact;      # whether or not to ignore branch lengths
my $ratio;        # horizontal to vertical ratio
my $type;         # clado or phylo
my $title;        # title

my $std_size   = 16; # standard font size for initial
my $std_branch = 20; # standard font size for initial
my $xoffset_from_root;

my $content = "";
my $cnt_tol = 0; #init
my $cnt_max = 1; #init
my $cnt_scale;
my $cnt_scale_mode;
my $cnt_circle_max_r = 20;

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Draw::Cladogram();
 Function: Builds a new Bio::Tree::Draw::Cladogram object 
 Returns : Bio::Tree::Draw::Cladogram
 Args    : -tree => Bio::Tree::Tree object
           -font => font name [string] (optional)
           -size => font size [integer] (optional)
           -top => top margin [integer] (optional)
           -bottom => bottom margin [integer] (optional)
           -left => left margin [integer] (optional)
           -right => right margin [integer] (optional)
           -tip => extra tip space [integer] (optional)
           -compact => ignore branch lengths [boolean] (optional)
           -ratio => horizontal to vertical ratio [integer] (optional)
           -branch => stem length [integer] (optinal)
           -title => title of the plot [string] (optional)

=cut

sub new {
	my ( $class, @args ) = @_;

	my $self = $class->SUPER::new(@args);
	(
		$t1,      $type,     $font, $size,    my $top, my $bottom,
		my $left, my $right, $tip,  $compact, $ratio,  my $branch, my $scale, $cnt_scale_mode, $title
	  )
	  = $self->_rearrange(
		[
			qw(TREE TYPE FONT SIZE TOP BOTTOM 
			    LEFT RIGHT TIP COMPACT RATIO BRANCH SCALE MODE TITLE)
		],
		@args
	  );
	$font    ||= "Arial,Verdana,Helvetica,Sans-serif";
	$type    ||= "PHYLO";
	$size    ||= 16; #px
	$top     ||= 30;
	$bottom  ||= 30;
	$left    ||= 30;
	$right   ||= 10;
	$tip     ||= 15;
	$compact ||= 1;
	$ratio   ||= 1;
	$branch  ||= 20; #branch length
	$scale   ||= "max";
	$cnt_scale_mode ||= "normal";
	$xoffset_from_root = 1;
	$title	 ||= "";

	# Roughly, a cladogram is set according to the following parameters.

	#################################
	# title                     # T #   $top (T, top margin)
	#        +---------+ XXX    #   #   $bottom (B, bottom margin)
	#        |                  #   #   $left (L, left margin)
	#        |                  #   #   $right (R, right margin)
	#   +----+                  #   #   $tip (X, extra tip space)
	#        |    +----+ XXXX   #   #   $width (total drawing width)
	#        |    |             #   #   $height (total drawing height)
	#        +----+             # Y #   $xstep (S, stem length)
	#             |             #   #   $ystep (Y, space between taxa)
	#             +----+ XX     #   #   $tiplen (string length of longest name)
	#                           # B #   $tipwidth (N, size of longest name)
	#################################
	# L         S       X  N  R #
	#############################

	my @taxa = $t1->get_leaf_nodes;

	$tipwidth1 = 0;

	foreach my $taxon (@taxa) {
		my $w = length( $taxon->id );
		$tipwidth1 = $w if ( $w > $tipwidth1 );
	}
	$tipwidth1 = $tipwidth1*$size*0.55;    #max label width

	my $ystep = $std_size + $std_size / 2.5;
	$height = $bottom + $ystep * ( scalar @taxa ) + $top;
	
	my $y = $top+5;
	for my $taxon ( @taxa ) {
		$yy{$taxon} = $y;
		$y += $ystep;
	}

	my @stack;
	my @queue; # postorder traversal
	push @stack, $t1->get_root_node;
	
	while (@stack) {
		my $node = pop @stack;
		push @queue, $node;
		foreach my $child ( $node->each_Descendent( -sortby => 'internal_id' ) ){
			push @stack, $child;
		}
	}
	@queue = reverse @queue;

	for my $node (@queue) {
		if( $node->has_tag("count") ){
			my ($cnt) = $node->get_tag_values("count");
			$cnt_tol += $cnt;
			$cnt_max = $cnt if $cnt > $cnt_max;
		}

		if ( !$node->is_Leaf ) {
			my @children = $node->each_Descendent;
			my $child    = shift @children;
			my $ymin     = my $ymax = $yy{$child};
			foreach $child (@children) {
				$ymax = $yy{$child} if $yy{$child} > $ymax;
				$ymin = $yy{$child} if $yy{$child} < $ymin;
			}
			$yy{$node} = ( $ymin + $ymax ) / 2;
		}
	}

	#scale factor
	$cnt_scale = $cnt_max;
	$cnt_scale = $cnt_tol if $scale eq 'total';

	if ($compact) {    # ragged right, ignoring branch lengths

		my @preorder = $t1->get_nodes( -order => 'depth' );
		$width = 0;
		
		my $root = shift @preorder;    # skip root
		$xx{$root} = $left;
		
		$xstep = $branch;
		
		for my $node (@preorder) {
			$xx{$node} = $xx{ $node->ancestor } + $xstep;
			$width = $xx{$node} if $xx{$node} > $width;
		}
		$width += $tip + $tipwidth1 + $right;
	}
	else {    # set to aspect ratio and use branch lengths if available

		my $total_height = ( scalar( $t1->get_leaf_nodes ) - 1 ) * $ystep;
		my $scale_factor = $total_height * $ratio / $t1->get_root_node->height;

		$width = $t1->get_root_node->height * $scale_factor;
		$width += $left + $xstep;
		$width += $tip + $tipwidth1 + $right;

		my @preorder = $t1->get_nodes( -order => 'depth' );
		my $root = shift @preorder;    # skip root
		$xx{$root} = $left;
		for my $node (@preorder) {
			my $bl = $node->branch_length;
			$bl = 1 unless ( defined $bl && $bl =~ /^-?\d+(\.\d+)?$/ );
			$xx{$node} = $xx{ $node->ancestor } + $bl * $scale_factor * $branch/$std_branch + $xoffset_from_root;
		}
	}

	if ( $type =~ /clado/i ) {
		@taxa = $t1->get_leaf_nodes;
		for my $taxa (@taxa) {
			$xx{$taxa} = $width - $right - $tipwidth1 - $tip;
		}
	}

	return $self;
}

=head2 print

 Title   : print
 Usage   : $obj->print();
 Function: Outputs $obj in Encapsulated PostScript (EPS) format 
 Returns : 
 Args    : 

=cut

sub print {
	my ( $self, @args ) = @_;

	# original $width and $height are use to viewBox attribute
	my ( $format, $caveh, $cavew ) = $self->_rearrange( [qw(FORMAT CAVEH CAVEW)], @args );
	$format ||= "svg";

	# original $width and $height are use to viewBox attribute
	$caveh ||= $height;    # svg height
	$cavew ||= $width;     # svg weight
	
	my $root1 = $t1->get_root_node;

	if ( $format =~ /svg/i ) {

		# print out the svg header and inline css
		$content .=
		  '<?xml version="1.0" encoding="UTF-8" standalone="no"?>' . "\n";

		#$content .= '<svg width="'.$cavew.'"  height="'.$caveh.'"
		$content .= '<svg width="' . $width . 'px"  height="' . $height . 'px"
			id="treesvg"
			xmlns="http://www.w3.org/2000/svg"
			xmlns:xlink="http://www.w3.org/1999/xlink"
			preserveAspectRatio="xMidYMin meet">
<style type="text/css">

<![CDATA[
	.tree polyline {
		fill: none;
		stroke: black;
	}
	
	.tree text {
		stroke: none;
		fill: black;
        font-family: '.$font.';
        font-szie: '.$size.'px;
	}
	
	.tree circle {
		fill-opacity: 1;
		fill: lightgrey;
	}
	
	.text-id-tag {
		fill-opacity: 1;
	}
	
	]]>
</style>
<g id="viewport" class="tree">';

		# add title
		$content .= "<text x='5' y='20'>$title</text>" if ($title);

	  # call $root1->each_Descendent to avoid root node
	  # because $root do not have ancestor, it will make wannings in __svghelper
		foreach my $child ( $root1->each_Descendent ) {
			&_svghelper( $child, $root1 );
		}
		$content .= '</g></svg>';
	}    #end elsif
	return $content;
}

=head2 __svghelper 

Title   : __svghelper
Usage   : &__svghelper(node obj)
Function: private method to print the subtree in SVG format
Returns :
Args    : none

=cut

sub _svghelper {

	my ( $node, $root1 ) = @_;

	#draw the line from node to ancestor node
	$content .= sprintf "<g id='%s'>", $node->internal_id;

	$content .=
	  sprintf "<polyline id='polyline%d' points='%f %f, %f %f, %f %f'/>\n",
	  $node->internal_id,
	  $xx{$node},
	  $yy{$node},
	  $xx{ $node->ancestor },
	  $yy{$node},
	  $xx{ $node->ancestor },
	  $yy{ $node->ancestor };

    #prepare to print
	my $cnt = 0;
	($cnt) = $node->get_tag_values("count") if $node->has_tag("count");
	my $cnt_scale_factor = 0;
	if($cnt>0){
		$cnt_scale_factor = $cnt/$cnt_scale;
		$cnt_scale_factor = _log10($cnt)/_log10($cnt_scale) if $cnt_scale_mode eq "log" && $cnt_scale > 1;
	}

	#score color
	my $score;
	my $style_for_score = ""; #style='fill:rgb(0,0,255)';
	my @max_color = (196,39,49);
	($score) = $node->get_tag_values("score") if $node->has_tag("score");
	if( defined $score && $score ne "NA" ){
		my @rgb;
		foreach my $color (@max_color){
			$color = 255 - int($score*(255-$color));
			push @rgb, $color;
		}
		$style_for_score = "style='fill:rgb($rgb[0],$rgb[1],$rgb[2]);'";
	}


	my $cnt_text = $cnt;
    $cnt_text = sprintf("%.2f", $cnt) if $cnt =~ /\./;
	$cnt_text .= ", $score" if defined $score && $score ne "NA";
    $cnt_text = $cnt>0 ? " ($cnt_text)" : "";
    
    my ($depth) = $node->get_tag_values("depth") if $node->has_tag("depth");

	if ( $node->is_Leaf ) {    #if node is leaf one, output the taxaon label
		$content .= sprintf "<circle id='circle%s' cx='%d' cy='%d' r='%d' %s stroke='grey' stroke-width='1'  />\n",
			$node->internal_id,
			$xx{$node},
			$yy{$node},
			$cnt_circle_max_r*$cnt_scale_factor,
			$style_for_score;

		$content .= sprintf '<text class="text-taxon" id="%s" x="%s" y="%s">%s</text>\n"',
			$node->id,
			( $xx{$node} + $tip ),
			( $yy{$node} + $size / 3 ),
			$node->id.$cnt_text;
	}
	else
	{  # The group tag format <g id="t.internal_id"> and draw a circle for click
		# travel sub tree
		foreach my $child ( $node->each_Descendent ) {
			&_svghelper($child);
		}
		
		$content .= sprintf "<circle id='circle%s' cx='%d' cy='%d' r='%d' stroke='grey' stroke-width='1'/>\n",
		        $node->internal_id,
		        $xx{$node},
		        $yy{$node},
		        $cnt_circle_max_r*($cnt_scale_factor);

		if( $cnt || $depth==1 || $depth==2 ){
			$content .= sprintf "<text id='idtext%s' class='text-id-tag' x='%s' y='%s'>%s</text>\n",
				$node->internal_id,
				$xx{$node} - length( $node->id.$cnt_text )*$size*0.58,
				$yy{$node} + ( $depth%2==0 ? $size*1.25 : -$size*0.75 ),
				$node->id.$cnt_text;
		}
	}
	$content .= "</g>\n";
}

sub _log10 {
	my $n = shift;
	return log($n)/log(10);
}

1;
