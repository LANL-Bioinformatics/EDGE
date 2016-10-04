#!/usr/bin/env perl
use strict;
use Getopt::Long;

my $fasta;
my $taxa;
my $csv;
my $rank = "strain";

GetOptions(
	"fasta=s"	=> \$fasta,
	"taxa=s"	=> \$taxa,
	"csv=s"		=> \$csv,
	"rank=s"	=> \$rank,
	"help|?"       => sub{Usage()}
);

sub Usage{
print <<"END";
Usage: perl $0 -fasta <contig> -csv <ctg_class.top.csv> -taxa NameOftaxa -rank TaxRank

	e.g: perl $0 -fasta contigs.fasta -csv output.ctg_class.top.csv -taxa "Zaire ebolavirus" -rank strain

	-rank    phylum, class, order, family, genus, species, or strain (default: strain)
END
exit;
}

&Usage unless ($fasta && $taxa && $csv);

##SEQ	RANK	ORGANISM	TAX_ID	PARENT	LENGTH	NUM_HIT	TOL_HIT_LEN	TOL_MISM	AVG_IDT	LINEAR_LEN	RANK_LINEAR_LEN	COV	SCALED_COV	ACC_COV_RGN	ACC_COV_LEN
#Ebola.plasma.IDBA_00000	superkingdom	Viruses	10239	root	14637	1	14637	0	1.0000	14637	14637	1.0000	1.0000	1..14637	14637
my %hit_id;
open (my $fh, $csv) or die "Cannot read $csv file \n";
while(<$fh>){
	chomp;
	my @tmp = split/\t/,$_;
	if ($tmp[1] =~ /$rank/ && $tmp[2] =~ /$taxa/){
		$hit_id{$tmp[0]}= $tmp[2];
	}
}
close $fh;

open (my $cfh, $fasta) or die "Cannot read $fasta file\n";
$/ = ">";
while (my $line=<$cfh>){
	$line =~ s/\>//g;
	my ($head, @seq) = split /\n/, $line;
	next if (!$head);
	my ($id,$desc) = $head =~ /^(\S+)(.*)/;
	if ( $hit_id{$id} ){ 
		$id = $id."$desc $hit_id{$id}";
		my $seq = join "", @seq;
		$seq =~ s/(.{70})/$1\n/g; 
		chomp $seq;
		print ">$id\n$seq\n";
	}
}    
$/="\n";
close $cfh;

	 
