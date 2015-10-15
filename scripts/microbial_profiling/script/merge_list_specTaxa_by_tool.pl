#!/usr/bin/perl
use strict;
use Getopt::Long;

my %opt;
my $res=GetOptions(\%opt,
                   'level|l=s',
                   'outdir|o=s',
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

print STDERR "Generating matrix...";
print "ID\t",join("\t",@datasets),"\n";
my $cnt=0;
foreach my $taxa ( sort { $order->{BY_TAXALIST}->{$a} <=> $order->{BY_TAXALIST}->{$b} ||
					      $order->{BY_TAXA_SUM}->{$a} <=> $order->{BY_TAXA_SUM}->{$b}
																					} keys %$matrix ) {
	print $taxa;
	foreach my $dataset ( @datasets ){
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
