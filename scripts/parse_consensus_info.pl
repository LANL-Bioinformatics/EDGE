#!/usr/bin/env perl
## Usage: parse_consensus_info.pl consensus_fasta changelog

my $fasta = $ARGV[0];  # consensus_fasta
my $log = $ARGV[1]; # changelog 
open (my $fh, "<", $fasta) or die "Cannot open $fasta";
my ($id, $seq, %info);
$/ = ">";
while(my $line=<$fh>){
		$line =~ s/\>//g;
		my ($id, @seq) = split /\n/, $line;
		next if (!$id);
		my $seq = join "", @seq;
		my ($lead_N) = $seq =~ /^(N*)/i;
		my ($tail_N) = $seq =~ /(N*)$/i;
		my $total_cap_N = $seq=~ tr/N/N/;
		my $total_small_n = $seq=~ tr/n/n/;
		my $total_N = $total_cap_N + $total_small_n;
		my $total_G = $seq=~ tr/G/G/;
		my $total_C = $seq=~ tr/C/C/;
		my $lead_cap_N = $lead_N =~ tr/N/N/;
		my $lead_small_n = $lead_N =~ tr/n/n/;
		my $tail_cap_N = $tail_N =~ tr/N/N/;
		my $tail_small_n = $tail_N =~ tr/n/n/;
		my $gap_num = $seq =~ s/(n{1,})/$1/g;
		my $ambigous_num = $seq =~ s/(N{1,})/$1/g;
		print "Genome Length:", length($seq),"\n";
		print "Total n/N count:", $total_N , "\n"; 
		print "Total N count:", $total_cap_N,"\n";
		print "Total n count:", $total_small_n,"\n";
		print "5\' n/N count:", length($lead_N) , "\n"; 
		print "5\' N count:", $lead_cap_N,"\n";
		print "5\' n count:", $lead_small_n,"\n";
		print "3\' n/N count:", length($tail_N) , "\n"; 
		print "3\' N count:", $tail_cap_N,"\n";
		print "3\' n count:", $tail_small_n,"\n";
		$gap_num = 0 if (!$gap_num);
		print "Ambiguous:",$ambigous_num,"\n";
		print "Gaps:", $gap_num,"\n";
		print "GC content:", sprintf("%.2f\n", ($total_G + $total_C) / length($seq) * 100 );
}
$/="\n";
close $fh;

my %info2;
open (my $fh2, "<", $log) or die "Cannot open $log";
while(<$fh2>){
		chomp;
		next if /RATIO_OF_REF/;
		my @field=split("\t",$_);
		if ($field[2] eq '-' or $field[3] eq '-'){
				$info2{$field[0]}->{indel_num}++;
		}elsif($field[3] eq 'n'){
                        $info2{$field[0]}->{gaps_num}++;
		}elsif($field[3] eq 'N'){
                        $info2{$field[0]}->{ambiguousN_num}++;
		}else{
				$info2{$field[0]}->{snps_num}++;
		}
}
close $fh2;
foreach my $key( keys %info2 ){
print "INDELS count:", $info2{$key}->{indel_num} , "\n";
print "SNPs count:" , $info2{$key}->{snps_num},"\n";
}



