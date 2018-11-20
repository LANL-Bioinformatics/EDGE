#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;
use JSON;
use FindBin qw($RealBin);
use POSIX qw/strftime/;
use lib $RealBin;
use gi2lineage;

my %opt;
my $workingDir = Cwd::getcwd();

GetOptions( \%opt,
            'prefix=s',
            'gottcha=s',
            'gottcha2=s',
            'bwa=s',
            'kraken=s',
            'metaphlan2=s',
            'pangia=s',
            'diamond=s',
	    'centrifuge=s',
            'taxfile=s',
            'help|?' 
);
if ( $opt{help} || scalar keys %opt == 0 ) { &Usage(); }

my $prefix = $opt{prefix} || "_";
my ($gottcha, $gottcha2, $bwa, $kraken, $metaphlan, $pangia, $diamond);
my @tools = ("gottcha","gottcha2","bwa","kraken","metaphlan2","pangia","diamond","centrifuge");
my $time = time;

my $html_file = "${prefix}_taxonomyDBtable.html";
my $tsv_file = "${prefix}_taxonomyDBtable.tsv";
my $acc_h;
my @acc_h2;
my @input; 
my @tsv_head;
my $th;
#&print_timeInterval($time, "Load Tanonomy info ... ");
#
my $time_string= strftime('%Y-%m-%d',localtime);

loadTaxonomy();

if ($opt{taxfile} && -e $opt{taxfile}){
	open (my $fh, "<", $opt{taxfile});
	my $header = <$fh>;
	#TaxID	NAME	GOTTCHA	GOTTCHA2	BWA	KRAKEN	METAPHLAN2	PANGIA	DIAMOND	
	chomp $header;
	@tsv_head = split("\t",$header);
	for my $i (2..$#tsv_head){
		my $tool = $tsv_head[$i];
		$th .= "<th>".uc($tool)."</th>\n";
	}
	
	while(<$fh>){
		chomp;
		my @array = split("\t",$_);
		my $taxId = $array[0];
		my $tax_name = $array[1];
		#print $taxId."\t".getTaxRank($taxId)."\n";	
		#print $_,"\taa $array[8]\n" if ( getTaxRank($taxId) eq "family"); 
		my $rank = getTaxRank($taxId);
		next if $rank and $rank =~ /family|class|order|phylum|kingdom|tribe/;
		my $tax_id_link = "<a href=\'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=$taxId\' target='_blank' >$taxId</a>";
		my $ncbi_name_link = ($tax_name)?"<a href=\'http://www.ncbi.nlm.nih.gov/genome/?term=\"$tax_name\"\' target='_blank' >$tax_name</a>":"";
		for my $j (2..$#tsv_head){
			$acc_h->{$taxId}->[0] = $tax_id_link;
			$acc_h->{$taxId}->[1] = $ncbi_name_link;
			$acc_h->{$taxId}->[$j]=($array[$j])?$array[$j]:"";
        	}
	}
	#exit;
}else{

	for my $tool (@tools){
		if ( $opt{$tool} && -e $opt{$tool}){
			push @input, $tool; 
		}
	}

	@tsv_head = ("TaxID","NAME");
	for my $i (1..scalar(@input)){
		my $tool = $input[$i-1];
		&get_info($opt{$tool},$tool,$i,scalar(@input));
	
		$th .= "<th>".uc($tool)."</th>\n";
		push @tsv_head, uc($tool);
	}
}
#exit;


open (my $tsv_out, ">" , $tsv_file);
print $tsv_out join("\t",@tsv_head,"\n");
#&print_timeInterval($time, "All Done.\n");
foreach my $tax (keys %$acc_h){
	my @content = @{$acc_h->{$tax}};
	push @acc_h2,[@content];
	my $string = join ("\t",@content);
	$string =~ s|<.+?>||g;
	print $tsv_out $string."\n";
}
close $tsv_out;

my $dataset = encode_json(\@acc_h2);
#print join("\n",split(/\[/,$dataset)),"\n";

open( my $ofh, ">", $html_file);

print $ofh <<"HTML";
<!DOCTYPE html>
<html>
<head>
	<meta charset='utf-8'>
	<meta name="description" content="EDGE : EDGE bioinformatics">
	<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/jquery.dataTables.css">
	<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js"></script>
	<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.js"></script>
</head>
<body>
<script>
var dataset = $dataset;

\$(document).ready(function() {
    \$('#taxonomy').DataTable( {
	"pageLength": 20,
        "data": dataset,
	"deferRender": true,
	"columnDefs": [
                {"className": "dt-left", "targets": [0]},
                {"className": "dt-center", "targets": "_all"}
        ]
    } );
} );
</script>
<div style='background:#50a253;'><h2 style='position:inherit; padding-left:20px;'>EDGE bioinformatics<span style='float:right; font-size: small; background:#50a253;'>$time_string</span></h2></div>
<table id="taxonomy" class="display" style="width:100%">
        <thead>
            <tr>
		<th>TaxID</th>
		<th>NAME</th>
                $th
            </tr>
        </thead>
        <tfoot>
            <tr>
		<th>TaxID</th>
		<th>NAME</th>
                $th
            </tr>
        </tfoot>
    </table>
</body></html>
HTML

close $ofh;

sub Usage {
print <<END;
Usage: $0  -gottcha db.ann ...

Options:
-prefix         Output file prefix (default: _)
-gottcha        ann file
-gottcha2       stats file
-bwa            ann file
-pangia         ann file
-metaphlan2      mpa_v20_m200_marker_info.txt
-kraken   	seqid2taxid.map
-diamond	seq.headers
-centrifuge     centrifuge-inspect -n > dbinfo.txt
-taxfile        existing taxonomy tsv file
-help           Show this usage

END

exit;
}

sub get_info{
	my $file=shift;
	my $tool=shift;
	my $index=shift;
	my $num_tool=shift;
	my %de_repliate=();
	my $test_num=0;
	open (IN, "<",$file);
	while(<IN>)
	{
  		chomp;
  		next if ($_ !~ /\|/ && $tool =~ /gottcha|pangia/);
		my $line=$_;
		my @array=split /\s+/,$_;
		my $id_string=$array[1];
		my $tax_name="";
		my $taxId;
		my $acc='';
		if ($tool eq "centrifuge"){
			$id_string = $line;
			$taxId = $1 if $line =~ /cid\|(\d+)/;
		}
		if ($tool eq "gottcha"){
			$id_string=$array[1];
		}
		if ($tool eq "pangia"){
			my @array2 = split /\|/,$array[1];
			$id_string=$array2[0]; 
		}
		if ($tool eq "gottcha2" ){
			$id_string=$array[2];
			$taxId = $id_string;
		}
		if ($tool eq "bwa"){
			next if $id_string !~ /\./;
		}
		if ($tool =~ /metaphlan/){
			$id_string=$array[0];
			if ($id_string =~ /^GeneID/){
				$tax_name = $1 if $line =~ /clade': '\w__(\w+)'/;
				$tax_name =~ s/_/ /g;
				$taxId = name2taxID($tax_name);
			}
			
		}
		if ($tool eq 'kraken'){
			$id_string=$array[0];
			if ($id_string =~ /taxid\|(\d+)/){
				$taxId = $1;
			}
		}
		if ($tool eq 'diamond'){
			$id_string=$array[0];
		}

		if (!$taxId){
			$acc = getAccFromSeqID($id_string) || '';
			$taxId = acc2taxID($acc);
		}
		
		if (!$taxId){
			print "Cannot find taxId for accession $acc $tax_name $id_string from $tool database\n";
			next;
		}

		$tax_name = getTaxName($taxId);
		next if $de_repliate{$taxId};
		$de_repliate{$taxId}=1;
		my $tax_id_link = "<a href=\'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=$taxId\' target='_blank' >$taxId</a>";
		my $ncbi_name_link = ($tax_name)?"<a href=\'http://www.ncbi.nlm.nih.gov/genome/?term=\"$tax_name\"\' target='_blank' >$tax_name</a>":"";
		for my $j (2..($num_tool+1)){
			$acc_h->{$taxId}->[0] = $tax_id_link;
			$acc_h->{$taxId}->[1] = $ncbi_name_link;
			$acc_h->{$taxId}->[$index+1]="v" if ($j == $index+1);
			$acc_h->{$taxId}->[$j]="" if (!$acc_h->{$taxId}->[$j]);
		}
		#last if ($test_num++ > 1000);
		#push @{$acc_h->{data}},["$ncbi_name_link","v"];
		#push @acc_h2,["$ncbi_name_link","v"];
		# print $acc_to_name,"\n";
	}
	close IN;
}

sub saveListToJason {
        my ($list, $file) = @_;
        open JSON, ">", "$file" or die "Can't write to file: $file\n";
        my $json = encode_json($list);
        print JSON $json;
        close JSON;
}


sub print_timeInterval{
    my $now = shift;
    my $msg = shift;
    $now = time - $now;
    my $string=sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
    print STDERR "[$string]  $msg\n";
}
