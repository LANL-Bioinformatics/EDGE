#!/usr/bin/perl
use strict;
use Getopt::Long;

my %opt;
my $res=GetOptions(\%opt,
                   'level|l=s',
                   'outdir|o=s',
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
	'gottcha-genDB-v'  => 4,
	'gottcha-speDB-v'  => 5,
	'gottcha-strDB-v'  => 6,
	'blastn'           => 7,
	'bwa-target'       => 8,
	'kraken'           => 9,
	'kraken_mini'      => 10,
	'sequedex-opg'     => 11,
	'sequedex-tol'     => 12,
	'metaphlan'        => 13,
	'metacv'           => 14,
	'bwa'              => 15
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

print STDERR "Generating matrix...\n";

print "ID\t",join("\t",sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools),"\n";
my $cnt=0;
foreach my $taxa ( sort { $order->{BY_TPLIST}->{$a} <=> $order->{BY_TPLIST}->{$b} ||
					      $order->{BY_TAXA_SUM}->{$a} <=> $order->{BY_TAXA_SUM}->{$b}
																					} keys %$matrix ) {
	print $taxa;
	foreach my $tool ( sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools ){
		my $number = $matrix->{$taxa}->{$dataset}->{$tool} ? $matrix->{$taxa}->{$dataset}->{$tool} : 0;
		if( $count->{$dataset}->{$tool} > 0 ){
            my $value = defined $opt{display_read_count} ? $number : $number/$count->{$dataset}->{$tool}*100;
			print "\t",$value;
		}
		else{
			print "\t0";
		}
	}
	print "\n";

	if( defined $opt{top} ){
		last if ++$cnt == $opt{top};
	}
}
print STDERR "...done\n";

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
      
  --help

';
	exit;
}
