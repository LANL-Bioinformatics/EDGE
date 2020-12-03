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
		   'tplist|t=s',
		   'filter|f=s',
                   'display_read_count',
		   'show_all_tp',
                   'top|t=i',
                   'help|h|?') || &usage();

&usage() if ( defined $opt{help} || scalar @ARGV < 1 );

my @listfile = @ARGV;
my $matrix;
my $count;
my $order;
my $data;
my $tools;
my $filters;
my %tool_order=(
	'gottcha2-genDB-v' => 1,
	'gottcha2-speDB-v' => 2,
	'gottcha2-strDB-v' => 3,
	'gottcha2-speDB-b' => 4,
	'gottcha-genDB-v'  => 5,
	'gottcha-speDB-v'  => 6,
	'gottcha-strDB-v'  => 7,
	'blastn'           => 8,
	'bwa-target'       => 9,
	'kraken'           => 10,
	'kraken2'          => 11,
	'kraken_mini'      => 12,
	'sequedex-opg'     => 13,
	'sequedex-tol'     => 14,
	'metaphlan'        => 15,
	'metaphlan2'       => 16,
	'metacv'           => 17,
	'bwa'              => 18,
	'pangia'           => 19,
	'centrifuge'       => 20,
	'diamond'          => 21 
);

#parsing filter
if( defined $opt{filter} && $opt{filter} ){
	my @terms = split /,/, $opt{filter};
	foreach my $term ( @terms ){
		my ($t, $cond) = $term =~ /^([\w\d-]+)(.*)$/;
		$filters->{$t}=$cond if( $t );
	}
}

############################
# parsing list files
############################

print STDERR "Parsing taxa result lists...\n";
my ($dataset,$tool);

print STDERR "\tTOOL\tFILTERED\tCLASS_TOL\tTOL_READS\n";

foreach my $file ( @listfile ){
	($dataset,$tool) = $file =~ /\d+_([^\/]+)\/([^\/]+)\/([^\/]+)\.[^\.]+$/;

	if( $opt{level} =~ /species/i ){
		next if $tool =~ /gottcha\d?-genDB/;
	}
	elsif( $opt{level} =~ /strain/i ){
		next if $tool =~ /gottcha\d?-speDB/;
		next if $tool =~ /gottcha\d?-genDB/;
	}

	$tools->{$tool}=1;
	my $filtered_taxa=0;

	open LIST, $file || die "Can't open $file\n";
	while(<LIST>){
		chomp;
		my @fields = split /\t/, $_;
		next if $fields[0] ne $opt{level};

		#apply filter to certain tool
		if( defined $filters->{$tool} ){
			my $cond = $filters->{$tool};
			if ( eval("$fields[2] $cond") ){
				$filtered_taxa++;
				next;
			}
		}

		$matrix->{$fields[1]}->{$dataset}->{$tool} = $fields[2];

		# counting total reads/abundance
		$count->{$dataset}->{$tool} = 0 unless defined $count->{$dataset}->{$tool};
		$count->{$dataset}->{$tool} += $fields[2];

		# counting number of classifications
		$count->{CLASS}->{$tool} = 0 unless defined $count->{CLASS}->{$tool};
		$count->{CLASS}->{$tool}++;
	}
	close LIST;

	print STDERR "\t$tool\t$filtered_taxa\t$count->{CLASS}->{$tool}\t$count->{$dataset}->{$tool}\n";
}

#############################
# Parsing true positive list
#############################

my $tplist;
if( defined $opt{tplist} ){
	print STDERR "Parsing specific taxa lists...";
	open LIST, $opt{tplist} || die "Can't open $opt{tplist}\n";
	my $specific_taxa_count=-10000000;
	while(<LIST>){
		chomp;
		my @fields = split /\t/, $_;
		my ($lvl) = $opt{level} =~ /^(\w)/;

		my $taxa;
		if( defined $opt{level} && $opt{level} eq 'strain' ){
			$taxa = $fields[0];
		}
		else{
			($taxa) = $fields[1] =~ /$lvl\__([^\|]+)/;
			$taxa =~ s/_/ /g;
		}

		$tplist->{$taxa} = $specific_taxa_count++;

		if( defined $opt{show_all_tp} ){
			$matrix->{$taxa}->{$dataset}->{show_all_tp} = 0 unless defined $matrix->{$taxa}->{$dataset};
		}

		$order->{BY_TPLIST}->{$taxa} = $specific_taxa_count++;
	}
	close LIST;
	print STDERR "...done.\n";
}

#############################
# Create orders
#############################

foreach my $tool ( keys %$tools ){
	my $order_cnt=0;
	foreach my $taxa ( sort { $matrix->{$b}->{$dataset}->{$tool} <=> $matrix->{$a}->{$dataset}->{$tool}  } keys %$matrix ) {
		$order->{BY_TOOL}->{$tool}->{$taxa} = $order_cnt++;
	}
}

my $order_cnt=0;
foreach my $taxa ( keys %$matrix ){
	my $sum=0;
	foreach my $tool ( keys %{$matrix->{$taxa}->{$dataset}} ){
		my $number = $matrix->{$taxa}->{$dataset}->{$tool} ? $matrix->{$taxa}->{$dataset}->{$tool} : 0;
		next unless $count->{$dataset}->{$tool};
		$sum += defined $opt{display_read_count} ? $number : $number/$count->{$dataset}->{$tool}*100;
	}
	$order->{BY_TAXA_SUM}->{$taxa} = -$sum;
}

#############################
# Generating matrix
#############################

if ( defined $opt{otu}){
	loadTaxonomy();
	open OTU, ">", "$opt{otu}";
	print OTU "#OTU ID\t",join("\t",sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools),"\ttaxonomy\n";
}

print STDERR "Generating matrix...\n";

my $fh;
if (defined $opt{output}){
	open $fh, ">" , "$opt{output}";
}else{
	$fh = *STDOUT ;
}
print $fh "ID\t",join("\t",sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools),"\n";
my $cnt=0;
my $unknow=1;
foreach my $taxa ( sort { $order->{BY_TPLIST}->{$a} <=> $order->{BY_TPLIST}->{$b} ||
		          $order->{BY_TAXA_SUM}->{$a} <=> $order->{BY_TAXA_SUM}->{$b}
			} keys %$matrix ) {
	my $value=0;
	print $fh $taxa,"\t" if( defined $opt{top} && $cnt < $opt{top});
	my @values=();
	foreach my $tool ( sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools ){
		my $number = $matrix->{$taxa}->{$dataset}->{$tool} ? $matrix->{$taxa}->{$dataset}->{$tool} : 0;
		if( $count->{$dataset}->{$tool} > 0 ){
			$value = defined $opt{display_read_count} ? $number : $number/$count->{$dataset}->{$tool}*100;
		}
		else{
			$value=0;
		}
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

sub usage {
	print STDERR '
USAGE: merge_list_specTaxa_by_dataset.pl <DATASET>/<TOOL>/<DATASET>.out.list [OPTIONS]

EXAMPLE: 
    merge_list_specTaxa_by_dataset.pl 1_CDC_10_clean/bwa/*.out.list 1_CDC_10_clean/metaphlan/*.out.list
    merge_list_specTaxa_by_dataset.pl \
		--level genus \
		--top 60 \
		--show_all_tp \
		--tplist /users/218817/bin/users/218817/bin/HMPRP_sT1-Mock_lineage.txt \
		--filter "bwa<10,blastn<10" \
		1_454_EVEN/*/*.out.list

OPTIONS:

  --level <RANK>

      Default is \'species\'. Any ranks from \'kingdom\' to \'strain\' are allowable.

  --tplist|-t <LIST>

      A list contains organism and lineage as true positive. The organisms on the list 
	  will be put on the top of the matrix in order.

	  Available lists:
	      HMP : /users/218817/bin/users/218817/bin/CDC_spiked_lineage.txt
		  CDC : /users/218817/bin/users/218817/bin/HMPRP_sT1-Mock_lineage.txt
		  300M NO-ERROR:
                /users/218817/scratch/projects/20130526_HC100_test/100genus_no_error_lineage.txt
          300M ERROR:
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus1-100species.gilist_lineage.txt
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus1-200species.gilist_lineage.txt
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus1-300species.gilist_lineage.txt
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus2-100species.gilist_lineage.txt
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus2-200species.gilist_lineage.txt
                /users/218817/scratch/projects/20130625_HC150/calc_fp_tp/100genus2-300species.gilist_lineage.txt

  --show_all_tp
      
      Show all organisms in tplist even not being reported by any tools.

  --display_read_count

	  The matrix displays the real value instead of relative abundance.

  --filter <EXP1>[,<EXP2>,<EXP3>...]

      One filter can be applied to a tool. The expression starts with the tool
	  name following by a operator and a value. Note that the input value is a
	  raw data not a normalized value.
	  E.g.:
	  	bwa<10
		bwa<5,blastn<5

  --top (INT)

      Only display a certain number of classifications.

  --order <ORDER1>[,<ORDER2>,<ORDER3>...]
      
      **** [not available yet] ****
      Default sorting order is TPLIST following by the sum of a organism.
 
  --output   matrix output filename
  --otu      output name for otu table format text file      
  --help

';
	exit;
}
