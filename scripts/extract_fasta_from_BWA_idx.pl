#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use JSON;

my $project_name;
my $outputDir="Reference";
my $ref_id_map;
my $genomes;
my $genomesFiles;
my $ref_json_file;
my $bwa_index_id_map;
my $NCBI_bwa_genome_index;

GetOptions(
   'o=s'      => \$outputDir,
   'map=s'    => \$ref_id_map,
   'ref_json=s'	=> \$ref_json_file,
   'bwa_id_map=s' => \$bwa_index_id_map,
   'bwa_genome_index=s' =>\$NCBI_bwa_genome_index,
   'genomesList=s'	=> \$genomes,
   'help|h'   => sub{usage()},
);

$ref_json_file ||= "$RealBin/../edge_ui/data/Ref_list.json";
$bwa_index_id_map ||= "$RealBin/../database/bwa_index/id_mapping.txt";
$NCBI_bwa_genome_index ||= "$RealBin/../database/bwa_index/NCBI-Bacteria-Virus.fna";

#########################################################
`mkdir -p $outputDir`;

my $outputDir_abs_path = Cwd::abs_path("$outputDir");
my $reffile;
my $gff_file;
my %ref_id;
my $refdir;
my $annotation="$outputDir/annotation.txt";
my @ref_abs;

my $list = &readListFromJson($ref_json_file); 
$refdir = "$outputDir_abs_path/reffiles";
`mkdir -p $refdir`;
if ($genomes){
	open (my $a_fh, ">$annotation") or die "Cannot write to $annotation\n";
	print $a_fh  "Name\tDescription\tURL\n";
	## extract ref genome from bwa index and fasta file
	if ( -f $genomes) # a list file
	{
		open(my $fh, $genomes) or die "Cannot read $genomes\n";
		while(<$fh>)
		{
			chomp;
			&extract_ref_seq($_,$list,$refdir,$a_fh);
			print "$refdir/$_.fna\n";
		}
		close $fh;
	}
	else
	{
		my @names = split /,/,$genomes;
		foreach (@names)
		{
			&extract_ref_seq($_,$list,$refdir,$a_fh);
			print "$refdir/$_.fna\n";
		}
	}	
	close $a_fh;
}

sub extract_ref_seq {
	my $name=shift;
	my $list=shift;
	my $dir=shift;
	my $annotion_fh=shift;
	my $out_file = "$dir/$name.fna";
	open (my $g_fh, ">$out_file") or die "Cannot write $out_file";
	my $s;
	my $count=0;
	foreach my $acc (@{$list->{$_}})
	{
		my $get = `grep $acc $bwa_index_id_map`;
		my @ref= split /\s+/,$get;
		my $extract_id = shift @ref;
		my ($gi) = $extract_id =~ /gi\|(\d+)/;
		my @seq = `samtools faidx $NCBI_bwa_genome_index \"$extract_id\"`;
		$seq[0] =~ s/\n/ $name\n/m;
		$s .= join("",@seq);
		print $annotion_fh "$name\t".join(" ",@ref)."\thttp://www.ncbi.nlm.nih.gov/nuccore/$gi\n" if ($count<1);
		$count++;
	}
	print $g_fh $s."\n";
	close $g_fh;
}

sub readListFromJson {
	my $json = shift;
	my $list = {};
	if( -r $json ){
		open JSON, $json;
		flock(JSON, 1);
  		local $/ = undef;
  		$list = decode_json(<JSON>);
  		close JSON;
	}
	return $list;
}
