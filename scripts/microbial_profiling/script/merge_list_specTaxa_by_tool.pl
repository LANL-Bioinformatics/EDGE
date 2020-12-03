#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib $RealBin;
use gi2lineage; 

my %opt;
my $res=GetOptions(\%opt,
                   'level|l=s',
                   'output|o=s',
                   'otu=s',
                   'prefix|p=s',
		   'taxalist=s',
                   'display_read_count',
                   'top|t=i',
                   'help|h|?') || &usage();

&usage() if ( defined $opt{help} || scalar @ARGV < 1 );

sub usage {
	print STDERR "
USAGE: merge_list_specTaxa_by_tool.pl (--level [rank] --taxalist [list] --display_read_count --top (INT) --h) *.out.list
\n";
	exit;
}

my @listfile = @ARGV;
my $matrix;
my $count;
my $order;
my $data;
my @datasets;
my $tool;
#
# parsing list files
#
my ($dataset,$tool);
print STDERR "Parsing taxa result lists...";
foreach my $file ( @listfile ){
	($dataset,$tool) = $file =~ /\d+_([^\/]+)\/([^\/]+)\/([^\/]+)\.[^\.]+$/;
	push @datasets, $dataset;

	open LIST, $file || die "Can't open $file\n";
	while(<LIST>){
		chomp;
		my @fields = split /\t/, $_;
		next if $fields[0] ne $opt{level};

		$matrix->{$fields[1]}->{$dataset}->{$tool} = $fields[2];
		$count->{$dataset}->{$tool} += $fields[2];
		$order->{BY_TAXALIST}->{$fields[1]} = 1000;
	}
	close LIST;
}
print STDERR "...done.\n";

#
# parsing taxa files
#
my $taxalist;
if( defined  $opt{taxalist} ){
	print STDERR "Parsing specific taxa lists...";
	open LIST, $opt{taxalist} || die "Can't open $opt{taxalist}\n";
	my $specific_taxa_count=-100000;
	while(<LIST>){
		chomp;
		my @fields = split /\t/, $_;
		my ($lvl) = $opt{level} =~ /^(\w)/;
		my ($taxa) = $fields[1] =~ /$lvl\__([^\|]+)/;
		$taxa =~ s/_/ /g;
		$taxalist->{$taxa} = $specific_taxa_count++;
		$order->{BY_TAXALIST}->{$taxa} = $specific_taxa_count++;
	}
	close LIST;
	print STDERR "...done.\n";
}

#
# Create order list
#
foreach my $dataset ( @datasets ){
	my $order_cnt=0;
	foreach my $taxa ( sort { $matrix->{$b}->{$dataset}->{$tool} <=> $matrix->{$a}->{$dataset}->{$tool}  } keys %$matrix ) {
		$order->{BY_DATASET}->{$dataset}->{$taxa} = $order_cnt++;
	}
}

#my $order_cnt=0;
#if( $taxalist ){
#	foreach my $taxa ( sort { $taxalist->{$a} <=> $taxalist->{$b} } keys %$taxalist ){
#		$order->{BY_TAXALIST}->{$taxa} = $order_cnt++;
#	}
#}
#else{
#	foreach my $taxa ( keys %$taxalist ){
#		$order->{BY_TAXALIST}->{$taxa} = 1;
#	}
#}

my $order_cnt=0;
foreach my $taxa ( keys %$matrix ){
	my $sum=0;
	foreach my $dataset ( @datasets ){
		my $number = $matrix->{$taxa}->{$dataset}->{$tool} ? $matrix->{$taxa}->{$dataset}->{$tool} : 0;
		$sum += $opt{display_read_count} ? $number : ($count->{$dataset}->{$tool} ? $number/$count->{$dataset}->{$tool}*100 : 0 ) ;
	}
	$order->{BY_TAXA_SUM}->{$taxa} = -$sum;
}

if ( defined $opt{otu}){
	loadTaxonomy();
	open OTU, ">", "$opt{otu}";
	print OTU "#OTU ID\t$tool\ttaxonomy\n";
}
print STDERR "Generating matrix...";
my $fh;
if (defined $opt{output}){
	open $fh, ">" , "$opt{output}";
}else{
	$fh = *STDOUT;
}
print $fh "ID\t",join("\t",@datasets),"\n";
my $cnt=0;
my $unknow=1;
foreach my $taxa ( sort { $order->{BY_TAXALIST}->{$a} <=> $order->{BY_TAXALIST}->{$b} ||
			  $order->{BY_TAXA_SUM}->{$a} <=> $order->{BY_TAXA_SUM}->{$b}
			} keys %$matrix ) {
	my $value=0;
	print $fh $taxa,"\t" if( defined $opt{top} && $cnt < $opt{top});
	my @values=();
	foreach my $dataset ( @datasets ){
		my $number = $matrix->{$taxa}->{$dataset}->{$tool} ? $matrix->{$taxa}->{$dataset}->{$tool} : 0;
		if( $count->{$dataset}->{$tool} > 0 ){
			$value = defined $opt{display_read_count} ? $number : $number/$count->{$dataset}->{$tool}*100;
		}else{$value=0;}
		push @values, $value;
	}
	print $fh join("\t",@values),"\n" if( defined $opt{top} && $cnt < $opt{top});

	if ( defined $opt{otu}){
		(my $taxId, my $lineage, $unknow) = get_lineage($taxa,$unknow);
		print OTU $taxId,"\t", join("\t",@values),"\t",$lineage,"\n";
	}
	++$cnt;
	
}

close $fh;
close OTU if (defined $opt{otu});
print STDERR "...done\n";

sub get_lineage{
	my $tax_name = shift;
	my $unknow_count = shift;
	my @rank=();
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
	my $taxId = name2taxID($tax_name);
	my $inTaxID = $taxId;
	my $rank = getTaxRank($taxId);
	my $lineage = "k__; p__; c__; o__; f__; g__; s__$tax_name";
	if (!$inTaxID){
		$inTaxID = "unknow.$unknow_count";
		$unknow_count++;
		return ($inTaxID, $lineage,$unknow_count);
	}

	while ($taxId){
		if( defined $major_level{$rank} ){
			$level{$major_level{$rank}} = $tax_name;
		}
		last if $tax_name eq 'root';
		$taxId = getTaxParent($taxId);
		$rank = getTaxRank($taxId);
		$tax_name = getTaxName($taxId);	
		$tax_name =~ s/ /_/g;
	}
	foreach my $lvl ( ('s','g','f','o','c','p','k') ){
		unshift @rank,"${lvl}__$level{$lvl}";
	}
	$lineage=join "; ", @rank;
	return ($inTaxID, $lineage,$unknow_count);
}
