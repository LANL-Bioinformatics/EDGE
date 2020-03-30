#!/usr/bin/env perl

my $file1=$ARGV[0]; # readsToRef.SNP_report.txt
my $file2=$ARGV[1]; # contigsToRef.SNP_report.txt

if (!@ARGV){print "$0 readsToRef.SNP_report.txt contigsToRef.SNP_report.txt\n";exit;}
my %hash;
open (my $fh,$file1);
my $header=<$fh>;
while(<$fh>){
	chomp;
	my @arrays=split /\s+/,$_;
	my $unique_id=$arrays[0].'_'.$arrays[1];
	$hash{$unique_id}=$_;
}
close $fh;

my %common_id;
my %hash2;
open (my $fh2,$file2);
$header=<$fh2>;
while(<$fh2>){
	chomp;
	@arrays2=split /\s+/,$_;
	my $unique_id2=$arrays2[0].'_'.$arrays2[1];
	if (! $hash{$unique_id2}){
		$hash2{$unique_id2}=$_;
	}else{
		$common_id{$unique_id2}=1;
	}
}
close $fh2;


foreach (keys %hash){
	if ($common_id{$_}){
		delete $hash{$_};
	}
}
my $output_gff;
if ($file1 =~ /SNP/){
	$output_gff =&convertCtg2refSnps(\%hash2,$header,"unique_contigsToRef");
	$output_gff.=&convertCtg2refSnps(\%hash,$header,"unique_readsToRef");
}else{
	$output_gff =&convertCtg2refIndels(\%hash2,$header,"unique_contigsToRef");
	$output_gff.=&convertCtg2refIndels(\%hash,$header,"unique_readsToRef");
}
print $output_gff;

sub convertCtg2refSnps {
	my $content = shift;
	my $header = shift;
	my $source = shift;
	my $out_text = "";
	my $id=0;
	chomp $header;
	my @header = split /\t/,$header;
	foreach my $key(keys %$content){
		my $record = $content->{$key};
		# [0]         [1]             [2]         [3]         [4]     [5]     [6]         [7]     [8]         [9]     [10]
		# Chromosome  SNP_position    Ref_codon   Sub_codon   aa_Ref  aa_Sub  Synonymous  Product CDS_start   CDS_end CDS_strand
	
		next if $record =~ /Merged with SNP /;
	
		my @temp = split /\t/, $record;
		#next if scalar @temp < 11;
	
		# id generater
		$id++;
	
		# fill null
		$temp[4] ||= "-";
		$temp[5] ||= "-";
		$temp[6] ||= "-";
		$temp[8] ||= "-";
		$temp[9] ||= "-";
		$temp[10] = ($temp[10] == -1) ? "-":"+";
		
		my $attr;

		for( my $i=0; $i<=$#temp; $i++){
			push @{$attr->{$header[$i]}}, $temp[$i];
		}
		
		#add ID attribute
		push @{$attr->{ID}}, "SNP_$id";
		push @{$attr->{Name}}, "SNP_$id";
		push @{$attr->{Description}}, "$temp[2]($temp[4])->$temp[3]($temp[5])";

		my $attr_str = &gff_attr_string($attr);
	
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n",
				$temp[0],
				$source,,
				"SNPs",
				$temp[1],
				$temp[1]+length($temp[2])-1,
				$temp[10],
				$attr_str
		;
	
	}

	close IN;
	return $out_text;
}

sub convertCtg2refIndels {
	my $content = shift;
	my $header = shift;
	my $source = shift;
	my $out_text = "";
	my $id=0;
	chomp $header;
	my @header = split /\t/,$header;
	foreach my $key ( keys %$content){
		my $record = $content->{$key};

		# [0] Chromosome [1] INDEL_position [2] Sequence [3] Length [4] Type [5] Product [6] CDS_start [7] CDS_end
	
		my @temp = split /\t/, $record;
		#next if scalar @temp < 8;
	
		# id generater
		$id++;
	
		my $attr;

		for( my $i=0; $i<=$#temp; $i++){
			push @{$attr->{$header[$i]}}, $temp[$i];
		}
		
		#add ID attribute
		push @{$attr->{ID}}, "indels_$id";
		push @{$attr->{Name}}, "indels_$id";
		push @{$attr->{Description}}, "$temp[4]:$temp[2] ($temp[3]bp);";

		my $attr_str = &gff_attr_string($attr);
	
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n",
				$temp[0],
				$source,
				$temp[4],
				$temp[1],
				($temp[4] eq "Insertion") ? $temp[1] : ($temp[1]+$temp[3]-1),
				"+",
				$attr_str
		;
	
	}

	close IN;
	return $out_text;
}

sub gff_attr_string {
	my $attr = shift;
	my $has_id = 0;
	my $gff_str = "";

	foreach my $attrn ( keys %$attr ){
		my @attrv = map { &encode($_) } @{$attr->{$attrn}};
		my $temp = join ",", @attrv;

		$attrn =~ s/[<>#]//g;
		if( $attrn =~ /^id$/i ){
			$attrn = uc($attrn); 
			$gff_str = "$attrn=$temp;$gff_str";
		}
		else{
			$gff_str .= "$attrn=$temp;";
		}
	}
	
	return $gff_str;
}

sub encode {
	my $str = shift;
	$str =~ s/,/%2c/g;
	$str =~ s/;/%3b/g;
	return $str;
}
