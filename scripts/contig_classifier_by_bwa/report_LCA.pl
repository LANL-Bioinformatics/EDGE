#!/usr/bin/perl -w
use strict;

my $acc;
my $name = "";

my %rank = (
	"unclassified" => 0,
	"superkingdom" => 1,
	"phylum"       => 2,
	"class"        => 3,
	"order"        => 4,
	"family"       => 5,
	"genus"        => 6,
	"species"      => 7,
	"strain"       => 8
);

print "##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tNUM_MERGED\tACC_COV_LEN\n";

while(<STDIN>){
    next if /^#/;
	my $line = $_;
    chomp;
	next if /\t0$/;
	my @temp = split /\t/, $_;
	#SEQ    RANK    ORGANISM    TAX_ID  PARENT  LENGTH  NUM_HIT TOL_HIT_LEN TOL_MISM    AVG_IDT LINEAR_LEN  RANK_LINEAR_LEN COV SCALED_COV  ACC_COV_RGN  ACC_COV_LEN
    # 0       1        2           3      4       5        6         7         8           9        10              11       12     13           14           15
	
	if( $temp[0] ne $name )
	{
		foreach my $rank ( sort {$rank{$b}<=>$rank{$a}} keys %rank ){
			next unless defined $acc->{$rank};
			if( scalar(@{$acc->{$rank}}) == 1 ){
				 print $acc->{$rank}[0];
				 last;
			}
		}
		$acc = {};
		$name = $temp[0];
	}
	
	push @{$acc->{$temp[1]}}, $line;
}
