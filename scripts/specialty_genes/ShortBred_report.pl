#!/usr/bin/env perl

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use JSON;
use Data::Dumper;

my $sb_result;
my $meta_data;
my $category;
my $project;
my $out_dir;
my $card_json="$RealBin/../../database/ShortBRED/AR/card.json";
print $card_json;
my $ARDB_tab_PATH="$RealBin/../../database/ShortBRED/AR/ARDB";

GetOptions( "SBresult=s"  => \$sb_result,
            "metadata=s"        => \$meta_data,
            "project=s"     => \$project,
            "category=s"        => \$category,
            "CARD=s" => \$card_json,
            "ARDB=s" => \$ARDB_tab_PATH,
            "out=s"     => \$out_dir,
            "help|?"    => sub {Usage()}
          );

sub Usage{
	print "perl $0 -SBresult ShortBred_result -metadata metaDataTable -category SpcealityGene -project ProjectNAME -out OutDir\n";
	print "  -category   AR or VF\n";
	print "  -CARD      card.json from https://card.mcmaster.ca/download\n";
	print "  -ARDB      PATH to ARDB tab files ftp://ftp.cbcb.umd.edu/pub/data/ARDB/ARDBflatFiles.tar.gz\n";
	exit;
}

if (!$sb_result || !$meta_data ) {&Usage;}

$project ||= "";
$category ||= "AR";
$out_dir ||= "Output";

`mkdir -p $out_dir`;

my $info;
my $aro_r;
my $drug_r;
if ($category eq "AR"){
	$drug_r = &load_ARDB_drugs($ARDB_tab_PATH);
	$aro_r = &load_card_json($card_json);
}

my %sb_result_T;

open (my $fh,$sb_result) or die "Cannot read $sb_result\n";
my $sb_result_header = <$fh>;
$sb_result_header=~ s/\r\n//;
$sb_result_header=~ s/\n//;
while(<$fh>){
	s/\r\n//;
	s/\n//;
	my @array = split /\s+/,$_;
	my $hit  = (scalar(@array)>2)? $array[2]:$array[1];
	@{$sb_result_T{$array[0]}->{row}} = @array if ( $hit > 0);
	$sb_result_T{$array[0]}->{count} = $hit if ( $hit > 0);
}
close $fh;

my $text_output = "$out_dir/${project}_${category}_genes_ShortBRED_table.txt";
my $html_output = "$out_dir/${project}_${category}_genes_ShortBRED_table.html";
my $json_output = "$out_dir/${project}_${category}_genes_ShortBRED_table.json";
open (my $ofh,">$text_output") or "die cannot write to $text_output\n";
open (my $ofhtml,">$html_output") or "die cannot write to $html_output\n";
open (my $ofhjson,">$json_output") or "die cannot write to $json_output\n";

open (my $fh2,$meta_data) or die "Cannot read $meta_data\n";
my $header = <$fh2>;
$header =~ s/\r\n//;
$header =~ s/\n//;
$header =~ s/\"//g;
$header =~ s/\./_/g;
$header =~ s/Family\t//;
$header =~ s/Source ID\t//;
$header =~ s/PMID\tAssertion//;
$header =~ s/Locus Tag\t//;
$header =~ s/Gene ID\t//;
$header = $sb_result_header."\t".$header;
my $sb_result_column_num = scalar(split/\t/,$sb_result_header);
my @headers = split /\t/, $header;
splice @headers, $sb_result_column_num, 0, "Resistant_To"  if ($category eq "AR");


print $ofh join("\t",@headers),"\n";
print $ofhtml "<!DOCTYPE html>
<html>
<head>
<script type=\"text/javascript\" src=\"http://code.jquery.com/jquery-1.12.0.min.js\"></script>
<script type=\"text/javascript\" src=\"https://cdn.datatables.net/1.10.11/js/jquery.dataTables.min.js\"></script>
<link rel=\"stylesheet\" type=\"text/css\" href=\"https://cdn.datatables.net/1.10.11/css/jquery.dataTables.min.css\">
<script>
\$(document).ready(function() 
    { 
        \$(\"#myTable\").DataTable({
		\"scrollX\": true,
		\"order\": [[ 1, \"desc\" ],[0,\"desc\"]]
	}); 
    } 
);
</script>
</head>
<body>
<h2>$project $category Genes ShortBRED result table</h2>
<table id=\"myTable\" class=\"stripe compact hover nowrap\" cellspacing=\"0\" width=\"100%\">
	<thead>
	<tr>
";
for my $i  (0..$#headers){
	print $ofhtml "		<th>$headers[$i]</th>\n";
	my $hash;
	$hash->{title}= $headers[$i];
	$hash->{data}= $headers[$i];
	push @{$info->{columns}},$hash;
}
print $ofhtml "	</tr>\n	</thead>\n";
 
my $column_num=0;
while(<$fh2>){
	s/\r\n//;
	s/\n//;
	s/\\\"//g;
	s/\"//g;
	my @array = split /\t/,$_;
	if ($category eq "AR" && ($sb_result_T{$array[1]}->{count} > 0) ){
		my $resistant_drugs = ($drug_r->{$array[1]})? join("; ",@{$drug_r->{$array[1]}}):"";
		push @{$sb_result_T{$array[1]}->{row}} , $resistant_drugs, $array[0],@array[2..$#array];
		$column_num = (scalar(@{$sb_result_T{$array[1]}->{row}}) > $column_num)? scalar(@{$sb_result_T{$array[1]}->{row}}):$column_num;
	}elsif($category eq "VF" && ($sb_result_T{$array[2]}->{count}>0) ){
		push @{$sb_result_T{$array[2]}->{row}} , @array[0..1],$array[3],@array[6..($#array)];
		$column_num = (scalar(@{$sb_result_T{$array[2]}->{row}}) > $column_num)? scalar(@{$sb_result_T{$array[2]}->{row}}):$column_num;
	}
}
close $fh2;


foreach my $id (sort {$sb_result_T{$b}->{count}<=>$sb_result_T{$a}->{count}|| $a cmp $b } keys %sb_result_T){
	my $data;
	if ($sb_result_T{$id}->{count}>0){
		my @tmp=@{$sb_result_T{$id}->{row}};
		# print raw table
		print $ofh join("\t",@tmp),"\n";
		# print html/json
		$tmp[0] = "<a href='http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID=$tmp[0]' target='_blank'>$tmp[0]</a>" if ($category eq "VF");
		if ($category eq "AR"){
			if ($tmp[-2] eq "ARDB"){
				if ($tmp[-6]){# drug field
					my @tmp2 = split /; /,$tmp[-6];
					foreach my $i (0..$#tmp2){
						$tmp2[$i] ="<a href='http://ardb.cbcb.umd.edu/cgi/search.cgi?db=B&and0=O&term=$tmp2[$i]' target='_blank'>$tmp2[$i]</a>";
					}
					$tmp[-6] = join("; ",@tmp2);
				}
				$tmp[0] = "<a href='http://ardb.cbcb.umd.edu/cgi/search.cgi?db=R&term=$tmp[0]' target='_blank'>$tmp[0]</a>";
			}else{ # WUSTL need ARO
				my $desc = $tmp[-1];
				my ($aro) = $desc  =~ /ARO:(\d+)/;
				if ( $aro_r->{$aro}){
					my $aro_id = $aro_r->{$aro};
					$tmp[0] = "<a href='https://card.mcmaster.ca/ontology/$aro_id' target='_blank'>$tmp[0]</a>";
				}
			}
		}
		$tmp[1] = sprintf("%.2f",$tmp[1]);
		print $ofhtml "	<tr>\n";
		for my $i (0..($column_num-1)){
			$tmp[$i] ||= "";
			# print raw table
			print $ofhtml "		<td>$tmp[$i]</td>\n";
			$data->{$headers[$i]}=$tmp[$i];
		}
		push @{$info->{data}},$data;
		print $ofhtml "	</tr>\n";
	}
}

print $ofhtml "</table></body></html>\n";

my $json = to_json($info);
print $ofhjson $json;

close $ofh;
close $ofhtml;
close $ofhjson;




sub load_card_json{
	my $card_json = shift;
	open( my $fh, '<', $card_json );
	my $json_text = <$fh>;
	$json_text =~ s/\r//;
	my $card_r = decode_json($json_text);
	my %aro;
	foreach my $key ( keys %$card_r){
		next if ($key =~ /^_/); # comment version ...
		my $aro_accession = $card_r->{$key}->{ARO_accession};
		my $aro_id =  $card_r->{$key}->{ARO_id};
		$aro{$aro_accession}=$aro_id;
		if (ref($card_r->{$key}->{ARO_category}) eq "HASH"){
			foreach my $key2 (keys %{$card_r->{$key}->{ARO_category}}){
				$aro_accession = $card_r->{$key}->{ARO_category}->{$key2}->{category_aro_accession};
				$aro_id = $key2;
				$aro{$aro_accession}=$aro_id;
			}
		}
	}
	close $fh;
	return \%aro
}

sub load_ARDB_drugs{
	my $ardb_path=shift;
	my $resistance_profile = "$ardb_path/resistance_profile.tab";
	my $ar_genes = "$ardb_path/ar_genes.tab";
	my %ardb;
	open( my $fh1, '<',$resistance_profile);
	while(<$fh1>){
		chomp;
		my ($id,$drug) = split /\t/,$_;
		push @{$ardb{$id}},$drug;
	}
	close $fh1;

	open (my $fh, '<' , $ar_genes);
	while (<$fh>){
		chomp;
		my @array = split /\t/,$_;
		# protein id map to ardb id;
		@{$ardb{$array[0]}}=@{$ardb{$array[1]}} if ($ardb{$array[1]});
	}
	close $fh;
	return \%ardb;
}
