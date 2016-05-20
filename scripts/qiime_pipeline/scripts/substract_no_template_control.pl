#!/usr/bin/env perl
# chienchi at lanl.gov
# 20150908

use strict;
use Getopt::Long;

my $input;
my $output;
my $ntc_ids;

GetOptions(
        "input|i=s"	=> \$input,
        "output|o=s"	=> \$output,
	"ntc|n=s"	=> \$ntc_ids,
        "help|?"	=> sub{Usage()}
);

sub Usage{
print <<"END";
Usage: perl $0 -i <otu_table> -ntc No_template_control_ids

	-input    Tab-delimited otu table file. e.g. otu_table_tabseparated.txt
	-ntc      No_template_control Sample IDs, sepearted by comma
	-output   Output to file

END
exit;
}

&Usage unless ($input && $ntc_ids);
my $ofh;
if ($output){
	open($ofh, ">$output") or die "Cannot write to $output file \n";
}else{
	open($ofh, ">&STDOUT") or die "Cannot wirte to STDOUT \n";
}

open (my $fh, $input) or die "Cannot read $input file \n";
my @header;
my $otu; # otu maxtrix
my $max_NTC;
	#$otu->{sampleid}->{otuID} = otu_count;
while(<$fh>){
        chomp;
	@header = split /\t/,$_ if (/OTU ID/);
	next if (/^#/);
        my @tmp = split /\t/,$_;
	my $max_ntc_count=0;
	for my $i (1..$#tmp){
		if ($ntc_ids =~ /$header[$i]/){
			$max_ntc_count = $tmp[$i] if ($max_ntc_count < $tmp[$i]);
			$max_NTC->{$tmp[0]}=$max_ntc_count;
		}else{
			$otu->{$header[$i]}->{$tmp[0]}=$tmp[$i];
		}
	}
	
}
close $fh;

my ( @no_NTC_header_index ) = grep {  $ntc_ids !~ /$header[$_]/ } 0..$#header;
my @no_NTC_header;
map { push @no_NTC_header, $header[$_] } @no_NTC_header_index;
print $ofh join("\t",@no_NTC_header),"\n";
 
foreach my $otu_id ( sort { $otu->{$header[1]}->{$b} <=> $otu->{$header[1]}->{$a} }  keys %{$otu->{$header[1]}}){
	my $count=0;
	my @sample_otu_count;
	my $taxa;
	foreach my $i (@no_NTC_header_index){
		next if ($i == 0); # otu id 
		if ($header[$i] eq "taxonomy")
		{
			$taxa= $otu->{$header[$i]}->{$otu_id};
		}else{
			my $substraction = $otu->{$header[$i]}->{$otu_id} - $max_NTC->{$otu_id};
			$otu->{$header[$i]}->{$otu_id} = ($substraction > 0)? $substraction : 0;
			$count = $count + $otu->{$header[$i]}->{$otu_id};
			push @sample_otu_count, sprintf ("%.1f",$otu->{$header[$i]}->{$otu_id});
		}
	}
	if ($count > 0 ){
		my $out_string = ($taxa)? join("\t",$otu_id,@sample_otu_count,$taxa): join("\t",$otu_id,@sample_otu_count); 
		print $ofh $out_string,"\n";
	}
}
close $ofh;

