#! /usr/bin/perl
use Getopt::Long;
use strict;

$|=1;
my ($inFastq, $maxReads);

GetOptions(
	'i=s', \$inFastq,
	'm=i', \$maxReads
);
&usage() unless -e $inFastq;

my $bases=0;
my $reads=0;
my ($header,$min,$max) = ("",164,0);
my $tol_avg_score=0;
my $tol_read_num=0;

if( $maxReads > 0 ){
	my $cmdout = `wc -l $inFastq`;
	($tol_read_num) = $cmdout =~ /^(\d+)/;
	
	if( $tol_read_num%4 ){
		die "ERROR: total number of lines($tol_read_num)%4 != 0. FASTQ format error!\n";
	}
	else{
		$tol_read_num = $tol_read_num/4;
	}
}

open IN, $inFastq or die "Can't open $inFastq: $!\n";
while ( <IN> )
{
	chomp;
	$header = $_;
	
	my $seq=<IN>;
	chomp $seq;
	
	#$seq =~ s/(\n|\r)$//;
	$bases += length($seq);
	$reads++;
	print STDERR "processing $reads reads in $inFastq\r" unless $reads%5000;
	
	my $tmp=<IN>; #+

	my $qual_seq=<IN>;
	chomp $qual_seq;
	
	#score
	my @qual = split //, $qual_seq;
	my $score = 0;
	foreach my $ascii ( @qual ){
		$score +=  ord($ascii);
		$min = ord($ascii) if ord($ascii) < $min;
		$max = ord($ascii) if ord($ascii) > $max;
	}
	
	$tol_avg_score += $score/length($seq);
	last if defined $maxReads && $reads >= $maxReads;
}
$tol_read_num = $reads unless $maxReads;
close IN;

my $platform = &guessPlatform($header);
my ($format,$offset) = &guessOffset($min,$max);

printf "%14s%14s%16s%15s%14s%11s%8s%18s\n","TOL_READS","PROCESSED","PLATFORM","TOL_BASES","AVG_LENGTH","AVG_SCORE","OFFSET","FASTQ_FMT";
printf "----------------------------------------------------------------------------------------------------------------\n";
printf "%14s%14s%16s%15s%14.2f%11.2f%8d%18s\n",
	$tol_read_num,
	$reads,
	$platform,
	$bases,
	($bases/$reads),
	($tol_avg_score/$reads-$offset),
	$offset,
	$format;

sub guessPlatform {
	my $header = shift;
	$header =~ s/^@//;

	my @guess_ilu = split /:/, $header;
	my @guess_pac = split /_/, $header;

	my @guess;

	if( $header =~ /MiSeq/i ){
		push @guess, "illumina";
	}
	elsif( $header =~ /HiSeq/i ){
		push @guess, "illumina";
	}
	elsif( scalar @guess_ilu >= 5 ){
		push @guess, "Illumina";
	}
	elsif( scalar @guess_pac >= 5 && $header =~ /^m.*_s\d+_p\d+/i ){
		push @guess, "pacbio";
	}
	elsif ( $header =~ /^[\w\d]{14}\s*/ ){
		push @guess, "ionTorrent";
	}
	else{
		push @guess, "unknown";
	}

#	if( scalar @guess_ilu == 10 ){
#		push @guess, "Casava 1.8";
#	}

	return join " ", @guess;
}


sub guessOffset {
	my ($min,$max) = @_;
	if( $min < 59 ){
		return ( "Illumina 1.8+", 33) if $max == 74;
		return ( "Sanger", 33);
	}
	elsif( $max > 74){
		return ( "Illumina 1.5+", 64) if $min >= 67;
		return ( "Illumina 1.3+", 64) if $min >= 64;
		return ( "Solexa", 64) if $min >= 59;
	}
	else{
		return ( "Unknown", "NA");
	}
}

#usage
sub usage {
print<<'end';
Usage: $0 -i <fastq file> [-m (PROCESS_NUMBER_OF__READS)]
end
exit;
}

1;
