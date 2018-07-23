# gi2lineage: This is a Perl package to convert gi/NCBI taxonomy ID
# to taxonomy id/name/lineage. Some codes are taken from the package
# from Krona (Brian Ondov, Nicholas Bergman, and Adam Phillippy
# , Battelle National Biodefense Institute (BNBI))
#
# Po-E (Paul) Li
# B-11, LANL
# 04/2013 v0.1
#
# Changes log:
# 04/07/2017 - Add getAccFromSeqID, getTaxIDFromAcc acc2taxID acc2name acc2rank acc2lineage
# 04/04/2014 - Add taxid2lineage method
#
# 10/29/2013 - fixed bug retrieving GI>536870909 (2GB limitation)
#              NOTE: THIS PACKAGE REQUIRED PERL v5.16 or ABOVE.
#
# 10/28/2013 - able to map old GI to updated GI by providing a mapping
#              table using loadUpdatedGI().
#              (eg: /users/218817/scratch/opt/src/krona/taxonomy/gi_updated.tab).
#
package gi2lineage;

use strict;
#use Storable;
use Getopt::Long;
use File::Basename;
use File::Path;
#use threads;
#use threads::shared;

use base 'Exporter';
use Cwd 'abs_path';

#############
# Exported #
############
our @EXPORT = qw
(
    acc2taxID
    acc2name
    acc2rank
    acc2lineage
    gi2taxID
    gi2name
    gi2rank
    gi2rank_taxid
    taxid2rank
    taxid2rank_taxid
    taxid2name
    taxid2lineage
    gi2lineage
    getTaxRank
    getTaxName
    getTaxDepth
    getTaxIDFromGI
    getTaxIDFromName
    getTaxIDFromAcc
    getAccFromSeqID
    getTaxParent
    name2taxID
    loadTaxonomy
    loadUpdatedGI
);

####################
# Global constants #
####################

my $libPath = `ktGetLibPath`;
my $taxonomyDir = "$libPath/../taxonomy";
my $fileTaxByAcc = 'all.accession2taxid.sorted';
my $DEBUG=0;

#################
# Lookup tables #
#################
my @taxDepths   ;
my @taxParents  ;
my @taxRanks    ;
my @taxNames    ;
my %nameTax     ;
my %updatedGI   ;
my %taxIDByAcc  ;
my %taxIDByGI   ;
my $taxIDByGIStr;
my %invalidAccs ;
my %missingAccs ;

#######################################################
# sub: name2taxID
# input scientifc name and return its taxonomy ID
# argument: taxonomy scientific name
# return: <taxonomy id>
#######################################################
sub name2taxID{
    checkTaxonomy();
    my ($name) = @_;
    return $nameTax{$name};

}
#################################################################
# sub: acc2taxID
# function: input a Accession number and return its taxomomy ID.
# argument: <accession>
# return: <taxomomy id>
#################################################################
sub acc2taxID {
    checkTaxonomy();
    my ($acc) = @_;
	return getTaxIDFromAcc($acc);
}
#################################################################
# sub: acc2name
# function: input a Accession number and return its taxomomy.
# argument: <accession>
# return: <taxomomy>
#################################################################
sub acc2name {
	checkTaxonomy();
	my ($acc) = @_;
	my $taxID = getTaxIDFromAcc($acc);
	my $name = getTaxName($taxID);
	return $name;
}
#################################################################
# sub: acc2rank
# function: convert gi to a specific rank
# argument: <accession>, <RANK>
# return: <taxomomy>
# example:
#      ac2rank("NC_000964","genus")
#      return: Bacillus
#################################################################
sub acc2rank {
	checkTaxonomy();
	my ($acc,$r) = @_;

	my @rank;
	my $taxID = getTaxIDFromAcc($acc);
	my $rank = getTaxRank($taxID);
	my $name = getTaxName($taxID);

	return "root" if $taxID == 1;
	return "root" if $r eq "root";
	return $name if $r =~ /^strain$/i;

	while( $taxID ){
		if( uc($rank) eq uc($r) ){
			return $name;
		}
    		last if $name eq 'root';

		$taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}
	return "";
}
#################################################################
# sub: gi2taxID
# function: input a GI number and return its taxomomy ID.
# argument: <GI>
# return: <taxomomy id>
#################################################################
sub gi2taxID {
    checkTaxonomy();
    my ($gi) = @_;
	return getTaxIDFromGI($gi);
}

#################################################################
# sub: gi2name
# function: input a GI number and return its taxomomy.
# argument: <GI>
# return: <taxomomy>
#################################################################
sub gi2name {
    checkTaxonomy();
	my ($gi) = @_;
	my $taxID = getTaxIDFromGI($gi);
    my $name = getTaxName($taxID);
	return $name;
}

sub taxid2name {
    checkTaxonomy();
	my ($taxID) = @_;
    my $name = getTaxName($taxID);
	return $name;
}

#################################################################
# sub: gi2rank
# function: convert gi to a specific rank
# argument: <GI>, <RANK>
# return: <taxomomy>
# example:
#      gi2rank("47552137","genus")
#      return: Bacillus
#################################################################
sub gi2rank {
    checkTaxonomy();
	my ($gi,$r) = @_;

	my @rank;
	my $taxID = getTaxIDFromGI($gi);
    my $rank = getTaxRank($taxID);
    my $name = getTaxName($taxID);

    return "root" if $taxID == 1;
    return "root" if $r eq "root";
    return $name if $r =~ /^strain$/i;

    while( $taxID )
    {
		if( uc($rank) eq uc($r) ){
			return $name;
		}
    	last if $name eq 'root';

        $taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}
	return "";
}

sub taxid2rank {
    checkTaxonomy();
	my ($id,$r) = @_;

    return "root" if $id == 1;
    return "root" if $r eq "root";
	
    my @rank;
	my $taxID = $id;
    my $rank = getTaxRank($taxID) || 'no rank';
    my $name = getTaxName($taxID) || '';;

    return $name if $r =~ /^strain$/i;
    
    while( $taxID )
    {
		if( uc($rank) eq uc($r) ){
			return $name;
		}
    	last if $name eq 'root';

        $taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}
	return "";
}
#################################################################
# sub: gi2rank_taxid
# function: convert gi to a specific rank
# argument: <taxid>, <RANK>
# return: <taxomomy id>
# example:
#      gi2rank_taxid("47552137","genus")
#      return: Bacillus
#################################################################
sub gi2rank_taxid {
    checkTaxonomy();
	my ($gi,$r) = @_;

	my @rank;
	my $taxID = getTaxIDFromGI($gi);
    my $rank = getTaxRank($taxID);
    my $name = getTaxName($taxID);

    return $taxID if $rank eq 'strain';
    return 1 if $rank eq "root";
    
    while( $taxID )
    {
		if( uc($rank) eq uc($r) ){
			return $taxID;
		}
    	last if $name eq 'root';

        $taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}
	return "";
}

sub taxid2rank_taxid {
    checkTaxonomy();
	my ($id,$r) = @_;

	my @rank;
	my $taxID = $id;
    my $rank = getTaxRank($taxID);
    my $name = getTaxName($taxID);

    return $taxID if $r eq $rank || ( $r eq 'strain' && $rank eq 'no rank');
    return 1 if $r eq "root";
    
    while( $taxID )
    {
		if( uc($rank) eq uc($r) ){
			return $taxID;
		}
    	last if $name eq 'root';

        $taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}
	return "";
}


#################################################################
# sub: gi2lineage(<GI>[,<PRINT_ALL_RANK>,<PRINT_STRAIN>,<REMOVE_SPACE>])
# function: input a GI number and return a full lineage following by Kingdom (k), Phylum (p), Class 
#           (c), Order (o), Family (f), Genus (g), Species (s) and/or Strain(n)
# argument:
#	<GI>       - gi number
#
#	<PRINT_ALL_RANK> (1|0) - (optional, default is 1)
#             1: display 7 major ranks even they are not available (display as no_rank).
#             0: display available ranks only.
#
#	<PRINT_STRAIN> (1|0) - (optional, default is 0)
#             1: display guessed strain name (n).
#             0: do not display strain name.
#	<REMOVE_SPACE> (1|0) - (optional, default is 1)
#             1: remove space from organism name to "_"
#             0: do not remove space from organism name
# return:
#   The full lineage string is joined of the ranks with "|" separator. Each rank represents by it's
#   initial as an abbreviation following by two underscores ("__") and the taxomomy.
#
# example:
#   gi2lineage("47552137")
#   return: k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_anthracis
#################################################################
sub taxid2lineage {
    my ($id,$print_all_rank,$print_strain,$replace_space2underscore) = @_;
    &id2lineage($id,$print_all_rank,$print_strain,$replace_space2underscore,"taxid");
}

sub gi2lineage {
    my ($id,$print_all_rank,$print_strain,$replace_space2underscore) = @_;
    &id2lineage($id,$print_all_rank,$print_strain,$replace_space2underscore,"gi");
}

sub acc2lineage {
    my ($id,$print_all_rank,$print_strain,$replace_space2underscore) = @_;
    &id2lineage($id,$print_all_rank,$print_strain,$replace_space2underscore,"acc");
}

sub id2lineage {
    checkTaxonomy();
	my ($id,$print_all_rank,$print_strain,$replace_space2underscore,$input_type) = @_;
	$print_all_rank           = 1 unless defined $print_all_rank;
	$print_strain             = 0 unless defined $print_strain;
    $replace_space2underscore = 1 unless defined $replace_space2underscore;
	my @rank;

	my %major_level = ( 
	    'superkingdom' => 'k', 
	    'phylum'       => 'p', 
	    'class'        => 'c', 
	    'order'        => 'o', 
	    'family'       => 'f', 
	    'genus'        => 'g', 
	    'species'      => 's' 
	);

	my %level = ( 
	    'k' => '', 
	    'p' => '', 
	    'c' => '', 
	    'o' => '', 
	    'f' => '', 
	    'g' => '', 
	    's' => ''
	);

    my $taxID;

    if( $input_type eq 'taxid' ){
        $taxID = $id;
    }
    elsif( $input_type eq 'gi' ){
	    $taxID = getTaxIDFromGI($id);
    }
    elsif ( $input_type eq 'acc'){
            $taxID = getTaxIDFromAcc($id);
    }

    my $rank = getTaxRank($taxID);
    my $name = getTaxName($taxID);
    my $str_name = $name;
    $str_name =~ s/ /_/g if $replace_space2underscore;

    while( $taxID )
    {
		if( defined $major_level{$rank} ){
			$name =~ s/ /_/g if $replace_space2underscore;
			$level{$major_level{$rank}} = $name;
		}
    	last if $name eq 'root';

        $taxID = getTaxParent($taxID);
		$rank = getTaxRank($taxID);
		$name = getTaxName($taxID);
	}

    my $last=$str_name;
	foreach my $lvl ( ('s','g','f','o','c','p','k') ){
		if( $print_all_rank == 0 ){
			next unless $level{$lvl};
		}
		$level{$lvl} = "$last - no_${lvl}_rank" unless $level{$lvl};
        $last=$level{$lvl};
		push @rank, "${lvl}__$level{$lvl}";
	}

    @rank = reverse @rank;
    push @rank, "n__$str_name" if $print_strain;

	return join "\|", @rank;
}

sub getTaxDepth
{
	my ($taxID) = @_;
	return $taxDepths[$taxID];
}

sub getTaxName
{
	my ($taxID) = @_;
	return $taxNames[$taxID]
}

sub getTaxParent
{
	my ($taxID) = @_;
	
	do
	{
		$taxID = $taxParents[$taxID] || 1;
	}
	while ( $taxID > 1 && $taxRanks[$taxID] eq 'no rank');
	
	return $taxID;
}

sub getTaxRank
{
	my ($taxID) = @_;
	return $taxRanks[$taxID];
}

sub getTaxIDFromGI
{
	my ($gi) = @_;
    $gi = $updatedGI{$gi} if defined $updatedGI{$gi};
	my $taxID;

    if( defined $taxIDByGI{$gi} )
    {
        return $taxIDByGI{$gi};
    }
    elsif( defined $taxIDByGIStr )
    {
        my $pos = $gi*4;
        eval {
            #the last one x2147483636L (2GB limitation)
            #$taxID = unpack "x${pos}L", $taxIDByGIStr;
            my $data = substr $taxIDByGIStr, $pos, 4;
            $taxID = unpack "L", $data;
        };

        die "ERROR: GI to TaxID table is out of date.\n" if $@;
    }
	elsif ( ! defined $taxIDByGI{$gi} )
	{
		open GI, "<$taxonomyDir/gi_taxid.dat" or die "ERROR: GI to TaxID data not found.  Was updateTaxonomy.sh run?\n";
		
		seek GI, $gi * 4, 0;
	
        my $data;
	
        read GI, $data, 4;
	
        $taxID = unpack "L", $data;
		
        close GI;
	}

    #make sure $taxID is in database
    if ( defined $taxRanks[$taxID] ){
    	$taxIDByGI{$gi} = $taxID;
    }

	return $taxIDByGI{$gi};
}

sub getTaxIDFromAcc
{
	my ($acc) = @_;
	
	if ( $acc =~ /^\d+$/ )
	{
		return $acc;
	}
	
	$acc =~ s/\.\d+$//;
	
	if ( defined $taxIDByAcc{$acc} )
	{
		return $taxIDByAcc{$acc};
	}
	
	my $size = -s "$taxonomyDir/$fileTaxByAcc";
	my $accCur="";
	my $taxID;
	
	if ( ! open ACC, "<$taxonomyDir/$fileTaxByAcc" )
	{
		print "ERROR: Sorted accession to taxID list not found. Was updateAccessions.sh run?\n";
		exit 1;
	}
	
	my $min = 0;
	my $max = $size;
	while ( $acc ne $accCur && $min < $max )
	{
		my $posNew = int(($min + $max) / 2);
		
		seek ACC, $posNew, 0;
		
		if ( $posNew != $min )
		{
			<ACC>; # eat up to newline
		}
		
		my $line = <ACC>;
		
		my $accNew;
		($accNew, $taxID) = split /\t/, $line;
		
		if ( $acc gt $accNew && $accCur ne $accNew && $accNew )
		{
			if ( $accNew )
			{
				$posNew = tell ACC;
			}
			
			$min = $posNew;
			
			if ( $min >= $max )
			{
				$max = $min + 1;
			}
		}
		else
		{
			$max = $posNew;
		}
		
		$accCur = $accNew;
	}
	
	close ACC;

	chomp $taxID;
	
	if ( $accCur ne $acc )
	{
		$missingAccs{$acc} = 1;
		$taxID = 0;
	}
	
	$taxIDByAcc{$acc} = $taxID;
	
	return $taxIDByAcc{$acc};
}

sub getAccFromSeqID
{
	my ($seqID) = @_;
	
	$seqID =~ /^>?(\S+)/;
	
	my $acc = $1;
	
	if ( $acc =~ /\|/ )
	{
		$acc = (split /\|/, $acc)[3];
	}
	
	if ( $acc !~ /^\d+$/ && $acc !~ /^[A-Z]+_?[A-Z]*\d+(\.\d+)?$/ )
	{
		$invalidAccs{$acc} = 1;
		return undef;
	}
	
	return $acc;
}



sub loadUpdatedGI
{
    my $giUpdatedFile = shift;

    open INFO, $giUpdatedFile or die "ERROR: Can't load updated GI table.\n";

    foreach(<INFO>)
	{
	    chomp;
	    my ($oldgi, $newgi) = split /\t/, $_;
		$updatedGI{$oldgi} = $newgi;
	}
}

sub loadTaxonomy
{
    my ($giTaxFile) = @_;

    #@taxParents = retrieve( "$taxonomyDir/taxParents.dat" );
    #@taxDepths  = retrieve( "$taxonomyDir/taxDepths.dat" );
    #@taxRanks   = retrieve( "$taxonomyDir/taxRanks.dat" );
    #@taxNames   = retrieve( "$taxonomyDir/taxNames.dat" );

    open INFO, "<$taxonomyDir/taxonomy.tab" or die
		"Taxonomy not found.  Was updateTaxonomy.sh run?";

    foreach( <INFO> )
	{
        my $line = $_; 
		chomp $line;
		my ($id, $depth, $parent, $rank, $name) = split /\t/, $line;
		
		$taxParents[$id] = $parent;
		$taxDepths[$id] = $depth;
		$taxRanks[$id] = $rank;
		$taxNames[$id] = $name;
		$nameTax{$name}=$id;
		$name =~ s/\W/ /g;
		$name =~ s/_/ /g;
		$nameTax{$name}=$id;
	}

	print STDERR "Done parsing taxonomy.tab\n" if $DEBUG;
	
    if ( $taxParents[2] == 1 )
	{
		ktDie( "Local taxonomy database is out of date. Update using updateTaxonomy.sh." );
	}
	close INFO;

    if( defined $giTaxFile && -e $giTaxFile ){
        my $ref = retrieve( $giTaxFile );
        %taxIDByGI = %$ref;
    }
    elsif( defined $giTaxFile && !-e $giTaxFile && $giTaxFile =~ /load/i ){
        if ( ! open GI, "<$taxonomyDir/gi_taxid.dat" ){
	    	print "ERROR: GI to TaxID data not found.  Was updateTaxonomy.sh run?\n";
	    	exit 1;
	    }
        #local $/ = undef;
	    #binmode GI;
        #$taxIDByGIStr = <GI>;
	    
        my $size = -s "$taxonomyDir/gi_taxid.dat";
        print STDERR "size of $taxonomyDir/gi_taxid.dat: $size\n" if $DEBUG;
        read(GI, $taxIDByGIStr, $size);
        print STDERR "size of taxIDByGIStr: ".length($taxIDByGIStr)."\n" if $DEBUG;
        print STDERR "Done loading gi_taxid.dat\n" if $DEBUG;
	    close GI;
    }
}

sub taxContains
{
	# determines if $parent is an ancestor of (or equal to) $child
	
	my ($parent, $child) = @_;
	
	my $depthParent = $taxDepths[$parent];
	
	while ( $taxDepths[$child] > $taxDepths[$parent] )
	{
		$child = $taxParents[$child];
	}
	
	return $parent == $child;
}

################
# Not exported #
################
sub ktDie
{
    my ($error) = @_;
    
    *STDOUT = *STDERR;
    printColumns('[ ERROR ]', $error);
    exit 1;
}

sub checkTaxonomy
{
	if ( ! @taxParents )
	{
		die 'Taxonomy not loaded. "loadTaxonomy()" must be called first.';
	}
}

sub taxonLink
{
	my ($taxID) = @_;
	
	return $taxID;
}

1;
