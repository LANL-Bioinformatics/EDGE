#!/usr/bin/env perl

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use JSON;
use Data::Dumper;

my $sb_hit;
my $sb_metadata;
my $faa;
my $contig_gff;
my $category;
my $outdir;
my $cov_table;
my $card_json="$RealBin/../database/CARD/card.json";

GetOptions( "SBhits=s"	=> \$sb_hit,
            "metadata=s"	=> \$sb_metadata,
            "faa=s"	=> \$faa,
            "gff=s"	=> \$contig_gff,
            "covTable=s" => \$cov_table,
            "CARD=s" => \$card_json,
            "category=s"	=> \$category,
            "out=s"	=> \$outdir,
            "help|?"	=> sub {Usage()}
          );

sub Usage{
	print "perl $0 -SBhits ShortBRED_SBhits -metadata metaDataTable -faa FAA -gff contig_gff -category SpcealityGene -out outdir\n";
	print "  -category  AR or VF\n";
	print "  -covTable  readsToContigs_coverage.table\n";
	print "  -CARD      card.json from https://card.mcmaster.ca/download\n";
	exit;
}

if (!$sb_hit || !$sb_metadata || !$faa || !$contig_gff){ &Usage; }

$category ||= "AR";
$outdir ||= "Output"; 
$outdir =~ s/\/$//;

`mkdir -p $outdir`;


my $metadata = &load_metadata($sb_metadata,$category);
my $annotation = &load_gff($contig_gff);
my $hit = &parse_sb_hit($sb_hit);

my $faaOutput = "$outdir/With_${category}_markers_hit.faa";
&get_faa_with_hit($faa,$hit,$faaOutput);

my $gffOutput = "$outdir/With_${category}_markers_hit.gff";
&convert_hit_gff($hit,$annotation,$metadata,$category,$gffOutput);

if ($cov_table){
	my $aro_r = ($category eq "AR")? load_card_json($card_json):"";
	my $contigTableOutput = "$outdir/With_${category}_contigs_coverage.json";
	&get_hit_contig_table($hit,$annotation,$metadata,$aro_r,$cov_table,$contigTableOutput);
}

# END MAIN
sub get_hit_contig_table{
	my $hit=shift;
	my $annotation=shift;
	my $metadata=shift;
	my $aro_r=shift;
	my $covTable_file=shift;
	my $output=shift;
	my %contig;
	my $info;
	foreach my $pid (sort {$a cmp $b} keys %{$hit}){
		foreach my $hit (@{$hit->{$pid}}){
			my ($markerID,$faastart,$faaend,$identity,$bitscore)=split /:/,$hit;
			my $cid = $annotation->{$pid}->{cid};
			my ($familyID) = $markerID =~ /(\w+)_[TQJ]M_#\d+/;
			$contig{$cid}{$familyID}=1;
		}
	}
	my $contigsizelimit=700; ## bigger 700 with JBrowse track
	open(my $fh, $covTable_file) or die "Cannot read $covTable_file\n";
	my $header = <$fh>;
	chomp $header;
	my @headers = split /\t/,$header;
	splice @headers, 1, 0,  "NCBI BLAST";
	splice @headers, 2, 0,  "$category Marker Family";
	$headers[0]= "CONTIG_ID";
	my $length_index;
	foreach my $i (0..$#headers){
		my $hash;
		$hash->{title}= $headers[$i];
		$hash->{data}= $headers[$i];
		push @{$info->{columns}},$hash;
		$length_index = $i if ($headers[$i] =~ /length/i);
	}
	my @tmp = split/\//,$output;
	my $projname = $tmp[-4];
	while(<$fh>){
		chomp;
		my $data;
		my @array=split(/\t/,$_);
		if ($contig{$array[0]}){ # with AR/VF hits
			my @familyID = keys %{$contig{$array[0]}};
			if ($category eq "AR"){
				foreach my $i (0..$#familyID){
					if ($metadata->{$familyID[$i]}->{source} eq "ARDB") {
						$familyID[$i]="<a href='http://ardb.cbcb.umd.edu/cgi/search.cgi?db=R&term=$familyID[$i]' target='_blank'>$familyID[$i]</a>";
					}else{ # WUSTL need ARO
						my $desc = $metadata->{$familyID[$i]}->{desc};
						my ($aro) = $desc  =~ /ARO:(\d+)/;
						if ( $aro_r->{$aro}){
							my $aro_id = $aro_r->{$aro};
							$familyID[$i]="<a href='https://card.mcmaster.ca/ontology/$aro_id' target='_blank'>$familyID[$i]</a>";
						}
					}
				}
			}else{ # VF
				for my $i(0..$#familyID){
					$familyID[$i]="<a href='http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID=$familyID[$i]' target='_blank'>$familyID[$i]</a>";
				}
			}
			splice @array, 1, 0,  "<a href='#' class='edge-contigBlast' >BLAST</a>";
			splice @array, 2, 0, join("; ",@familyID);
			if ($array[$length_index] >= $contigsizelimit){
				my $end = ($length_index)? $array[$length_index] : $contigsizelimit;
				my $SpecialtyGeneTrack = ($category eq "AR")? "%2CARmarker":"%2CVFmarker";
				$array[0]="<a href='JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA${SpecialtyGeneTrack}&loc=$array[0]%3A1..$end' target='_blank'>$array[0]</a>" if ($projname);
			}
			foreach my $i (0..$#array){
				$data->{$headers[$i]}=$array[$i];
			}
			push @{$info->{data}},$data;
		}
	}

	open (my $ofh, ">$output") or die "Cannot write to $output\n";
	my $json = to_json($info);
	print $ofh $json;
	close $ofh;
}

sub convert_hit_gff{
	my $hit=shift;
	my $annotation=shift;
	my $metadata=shift;
	my $category=shift;
	my $output=shift;
	my $num=1;
	open (my $ofh, ">$output") or die "Cannot write to $output\n";
	#test_AR_VF_081  barrnap:0.4.2   rRNA    2       108     .       +       .       ID=test_AR_VF_03205;locus_tag=test_AR_VF_03205;product=5S ribosomal RNA
	foreach my $pid (sort {$a cmp $b} keys %{$hit}){
		foreach my $hit (@{$hit->{$pid}}){
			my ($markerID,$faastart,$faaend,$identity,$bitscore)=split /:/,$hit;
			my $cid = $annotation->{$pid}->{cid};
			my $start = $annotation->{$pid}->{start} + $faastart * 3 - 3 ;
			my $end = $annotation->{$pid}->{start} + $faaend * 3 - 1 ; 
			my $strand = $annotation->{$pid}->{strand};
			my ($familyID) = $markerID =~ /(\w+)_[TQJ]M_#\d+/;
			my $note = $metadata->{$familyID}->{classification};
			print $ofh $cid."\tShortBRED\t${category}marker\t$start\t$end\t\.\t$strand\t0\t";
			print $ofh "ID=${category}_$num;NAME=$markerID;Note=$note;\n";
			$num++;
		}
	}
	close $ofh;
}

sub parse_sb_hit{
	my $sb_hit=shift;
# marker             query faa  id  alnlen  mismatch  gapopen qstart qend sstart send evalue bitscore
#NP_418505_TM_#04        test_AR_VF_03049        100.0   9       0       0       1       9       27      35      -0.2    24.3
	my %uniq_aa;  # hash of array =>   $uniq_aa = [ $markerID:$faastart:$faaend:$identity:$bitscore ; 
	      #					$markerID:$faastart:$faaend:$identity:$bitscore ]
	open (my $fh, $sb_hit) or die "Cannot read $sb_hit\n";
	while(<$fh>){
		s/\r\n//;
		s/\n//;
		my ($markerID,$faaId,$identity,$alnLen,$mismatch,$gapopen,$mstart,$mend,$faastart,$faaend,$evalue,$bitscore) = split /\t/,$_;
		push @{$uniq_aa{$faaId}}, "$markerID:$faastart:$faaend:$identity:$bitscore";
	}
	close $fh;
	return \%uniq_aa;
}

sub get_faa_with_hit{
	my $inFasta = shift;
	my $list= shift;
	my $outputfile = shift;
	open (OUT, ">$outputfile") or "Cannot write to $outputfile\n";
	$/ = ">";
	open (FASTA, "<$inFasta") or "Cannot read $inFasta\n";
	my $junk = (<FASTA>);
	while (my $frecord = <FASTA>) {
        	chomp $frecord;
        	my ($fdef, @seqLines) = split /\n/, $frecord;
		my ($pid) = $fdef =~ /^(\S+)/;
        	my $seq = join '', @seqLines;
		$seq =~ s/(.{70})/$1\n/g;
		chomp $seq; 
        	if ($list->{$pid}){
			print OUT ">$fdef ".join(";",@{$list->{$pid}})."\n$seq\n"; 
		}
	}
	close FASTA;
	close OUT;
	$/ = "\n";
}

sub load_gff{
	my $file=shift;
	my %annotation;
	open (my $fh, $file) or die "Cannot read $file";
	while(<$fh>){
		chomp;
		my @array = split /\t/,$_;
		next if ($array[2] ne "CDS");
		last if (/^>/);
		my ($pid) = $array[8] =~ /ID=([^;]*)/;
		$annotation{$pid}->{start} = $array[3];
		$annotation{$pid}->{end} = $array[4];
		$annotation{$pid}->{cid} = $array[0];
		$annotation{$pid}->{strand} = $array[6];
		$annotation{$pid}->{strand} = $array[6];
	}
	close $fh;
	return \%annotation;
}

sub load_metadata{
	my $file=shift;
	my $category=shift;
	my %metadata;
	open (my $fh, $file) or die "Cannot read $file";
	my $header = <$fh>;
	while(<$fh>){
		s/\r\n//;
	        s/\n//;
        	s/\\\"//g;
        	s/\"//g;
		my @array = split /\t/,$_;
		if ($category eq "AR"){
			$metadata{$array[1]}->{classification}=$array[2];
			$metadata{$array[1]}->{source}=$array[4];
			$metadata{$array[1]}->{desc}=$array[5];
		}elsif($category eq "VF"){
			$metadata{$array[2]}->{classification}=$array[10];
			$metadata{$array[2]}->{desc}=$array[9];
			$metadata{$array[2]}->{product}=$array[8];
		}
	}
	close $fh;
	return \%metadata;
}

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
	return \%aro;
}
