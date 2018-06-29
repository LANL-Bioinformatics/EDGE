#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use JSON;

my $contig;
my @pair_read;
my $single_end_read;
my $project_name;
my $outputDir="SNP_phylogeny";
my $ref_id_map;
my $SNPdbName;
my $genomes;
my $genomesFiles;
my $reference;
my $ref_json_file;
my $bwa_index_id_map;
my $treeMaker="FastTree";
my $kingdom = "Bacteria";
my $bootstrap=0;
my $bootstrap_n=100;
my $PosSelect;
my $NCBI_bwa_genome_index;
my $modeltest=0;
my $numCPU = 4;
my $nanopore=0;
my $aligner="bowtie";
my $version = "0.3";

my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/..";

GetOptions(
   'c=s'      => \$contig,
   'p=s{2}'   => \@pair_read,
   's=s'      => \$single_end_read,
   'n=s'      => \$project_name,
   'o=s'      => \$outputDir,
   'map=s'    => \$ref_id_map,
   'db=s'     => \$SNPdbName,
   'ref_json=s'	=> \$ref_json_file,
   'bwa_id_map=s' => \$bwa_index_id_map,
   'bwa_genome_index=s' =>\$NCBI_bwa_genome_index,
   'tree=s'	=>\$treeMaker,
   'kingdom=s' => \$kingdom,
   'genomesList=s'	=> \$genomes,
   'genomesFiles=s' => \$genomesFiles,
   'modeltest',  => \$modeltest,
   'bootstrap'	=>\$bootstrap,
   'bootstrap_n=i' =>\$bootstrap_n,
   'PosSelect=s'  => \$PosSelect,
   'reference=s'  => \$reference,
   'nanopore'  => \$nanopore,
   'cpu=i'    => \$numCPU,
   'help|h'   => sub{usage()},
);

$project_name ||= "${SNPdbName}_SNP_phylogeny";
$ref_id_map ||= "$RealBin/../database/SNPdb/reference.txt";
$ref_json_file ||= "$RealBin/../edge_ui/data/Ref_list.json";
$bwa_index_id_map ||= "$RealBin/../database/bwa_index/id_mapping.txt";
$NCBI_bwa_genome_index ||= "$RealBin/../database/bwa_index/NCBI-Bacteria-Virus.fna";

unless($contig || @pair_read || $single_end_read) {print "\nPlease provide contig or reads files\n\n"; &usage(); exit;};

`mkdir -p $outputDir`;

my $contig_abs_path = Cwd::abs_path("$contig");
my $outputDir_abs_path = Cwd::abs_path("$outputDir");
my $random_ref=0;
my $reffile;
my $gff_file;
my $updateSNP=2;
my $data_type;
my $read_type;
my $cdsSNPS=0;
my $treeMethod = ($treeMaker =~ /FastTree/i)? 1:2;
my $molecular_analysis=0;
my $kingdom_code=0; #Bacteria
my %ref_id;
my $refdir;
my $annotation="$outputDir/annotation.txt";
#Name	Description	URL
#Ecoli_042	Escherichia coli 042, complete genome	http://www.ncbi.nlm.nih.gov/nuccore/387605479

## prepare reference SNPdb
if ($genomes || $genomesFiles)
{
    $updateSNP=1;
	#my $list = &readListFromJson($ref_json_file); 
	#my @ref_list = keys %$list;
	$refdir = "$outputDir_abs_path/reffiles";
	`mkdir -p $refdir`;
	if ($genomes){
		if ( ! -e $bwa_index_id_map || ! -s $NCBI_bwa_genome_index){
			print "Cannot Find $bwa_index_id_map \n or \n Cannot find $NCBI_bwa_genome_index\n";
			exit;
		}
		open (my $a_fh, ">$annotation") or die "Cannot write to $annotation\n";
		print $a_fh  "Name\tDescription\tURL\n";
		## extract ref genome from bwa index and fasta file
		if ( -f $genomes) # a list file
		{
			open(my $fh, $genomes) or die "Cannot read $genomes\n";
			while(my $line=<$fh>)
			{
				chomp $line;
				#my @ref_name = grep { /$line/ } @ref_list;
				#&extract_ref_seq($ref_name[0],$list,$refdir,$a_fh);
				&extract_ref_seq($line,$refdir,$a_fh);
			}
			close $fh;
		}
		else
		{
			my @names = split /,/,$genomes;
			foreach my $name (@names)
			{
				#my @ref_name = grep { /$name/ } @ref_list;
				#&extract_ref_seq($ref_name[0],$list,$refdir,$a_fh);
				&extract_ref_seq($name,$refdir,$a_fh);
			}
		}	
		close $a_fh;
	}
	if($genomesFiles){
		map{
			my ($name,$path,$suffix)=fileparse("$_",qr/\.[^.]*/);
			open (my $fh, $_) or die "Cannot read $_\n";
			open (my $ofh, ">$refdir/$name.fna") or die "Cannot write $refdir/$name.fna\n";
			print $ofh ">$name\n";
			while(my $line=<$fh>){
				next if($line =~ />/);
				if($line =~ /\w+/){
					print $ofh $line;
				}
			}
			close $fh;
			close $ofh;
		} split /,/,$genomesFiles;
	}
	if($reference){
		my @tmpl = `ls -d $EDGE_HOME/database/NCBI_genomes/$reference*`;
		chomp @tmpl;
                my @gfiles = `ls -S $tmpl[0]/*gbk $tmpl[0]/*gbff 2>/dev/null`;
		my @gfffiles;
		foreach my $gbk (@gfiles){
			chomp $gbk;
			my $gbk_basename=basename($gbk);
			system("$RealBin/genbank2gff3.pl -e 3 --outdir stdout $gbk > $refdir/$gbk_basename.gff");
			push @gfffiles, "$refdir/$gbk_basename.gff";
		}
		my $cat_cmd="$RealBin/cat_gff.pl -i ". join(" ",@gfffiles) . "> $refdir/$reference.gff";
		system($cat_cmd);
		unlink @gfffiles;
		$reffile = "$reference.fna";
		$gff_file = "$reference.gff";
		$cdsSNPS=1;
		$random_ref=1;
	}
}
else{
	## precompute SNPdb
	## Check species in SNPdb
	$random_ref=1;
	$cdsSNPS=0;
	$refdir = "$outputDir_abs_path/files";
	open (my $mapfh,$ref_id_map) or die "Cannot open $ref_id_map\n";
	while(<$mapfh>)
	{
	    chomp;
            next if (/^#/);
	    my($id,$ref) = split /\s+/,$_;
	    $ref_id{$id}=$ref;
	}
	close $mapfh;
	$reffile= $ref_id{$SNPdbName}.".fna";
	$gff_file= $ref_id{$SNPdbName}.".gff";
	my ($name,$path,$suffix)=fileparse("$ref_id_map",qr/\.[^.]*/);
	$cdsSNPS=1 if ( -e "$path/$SNPdbName/files/$gff_file");
	my $current_db = join(", ",keys(%ref_id));
	unless($SNPdbName) {print "\nPlease specify a db_Name in the SNPdb\nCurrent available db are $current_db.\n\n"; &usage(); exit;}

	if (!$ref_id{$SNPdbName}) 
	{
    	print "\nThe SNPdbName=$SNPdbName SNPdb is not available.\nCurrent available db are $current_db.\n\n"; 
    	exit;
	}
}

## Prepare contig and reads fastq

my $control_file= "$outputDir_abs_path/phame.ctrl";

if (@pair_read)
{
	$read_type += 2;
    my $R1_abs_path = Cwd::abs_path("$pair_read[0]");
    my $R2_abs_path = Cwd::abs_path("$pair_read[1]");
    system("ln -sf $R1_abs_path $outputDir_abs_path/${project_name}_R1.fastq");
    system("ln -sf $R2_abs_path $outputDir_abs_path/${project_name}_R2.fastq");
}
if ($single_end_read)
{
    $read_type += 1;
    my $S_abs_path = Cwd::abs_path("$single_end_read");
    system("ln -sf $S_abs_path $outputDir_abs_path/${project_name}.fastq");
}
if ($contig)
{
    system("ln -sf $contig_abs_path $outputDir_abs_path/${project_name}.contig");
} 

$data_type = 1 if ( !(@pair_read || $single_end_read) && $contig && (!$genomes||!$genomesFiles)) ;
$data_type = 2 if ( (@pair_read || $single_end_read) && !$contig && (!$genomes||!$genomesFiles)) ;
$data_type = 3 if ( !(@pair_read || $single_end_read) && $contig && ($genomes||$genomesFiles)) ;
$data_type = 4 if ((@pair_read || $single_end_read) && !$contig && ($genomes||$genomesFiles)) ;
$data_type = 5 if ((@pair_read || $single_end_read) && $contig && (!$genomes||!$genomesFiles)) ;
$data_type = 6 if ((@pair_read || $single_end_read) && $contig && ($genomes||$genomesFiles)) ;
$modeltest = 0 if ($treeMaker =~ /FastTree/i);
$molecular_analysis = 1 if ($PosSelect =~ /PAML/i);
$molecular_analysis = 2 if ($PosSelect =~ /HyPhy/i);
$kingdom_code = 1 if ($kingdom =~ /Virus/i || $SNPdbName =~ /Virus/i);
$kingdom_code = 2 if ($kingdom =~ /Eukaryote/i);
$aligner="minimap2" if ($nanopore);


## Prepare Reference and control file
    system("cp -R $RealBin/../database/SNPdb/${SNPdbName}/* $outputDir_abs_path/.") if (! $genomes && $SNPdbName);
    
    open (my $fh, ">$control_file") or die "Cannot write $control_file\n";
    print $fh <<"CONTRL";
       refdir = $refdir  # directory where reference files are located
      workdir = $outputDir_abs_path  # directory where contigs/reads files are located and output is stored
    reference = $random_ref  # 0:pick a random reference; 1:use given reference
      reffile = $reffile  # reference species to use
      project = $project_name
      cdsSNPS = $cdsSNPS  # 0:no cds SNPS; 1:cds SNPs
    FirstTime = $updateSNP  # 1:yes; 2:update existing SNP alignment
         data = $data_type  # *See below 0:only complete(F); 1:only contig(C); 2:only reads(R); 
                   # 3:combination F+C; 4:combination F+R; 5:combination C+R; 
                   # 6:combination F+C+R; 7:realignment  *See below 
        reads = $read_type  # 1: single reads; 2: paired reads; 3: both types present;
      aligner = $aligner # support bowtie/bwa/minimap2
         tree = $treeMethod  # 0:no tree; 1:use FastTree; 2:use RAxML; 3:use both;
    modelTest = $modeltest  # 0:no; 1:yes; # Only used when building a tree using RAxML
    bootstrap = $bootstrap  # 0:no; 1:yes;  # Run bootstrapping  *See below
            N = $bootstrap_n  # Number of bootstraps to run *See below   
            
    PosSelect = $molecular_analysis  # 0:No; 1:use PAML; 2:use HyPhy; 3:use both *See below
         code = $kingdom_code  # 0:Bacteria; 1:Virus; 2:Eukaryote
        clean = 1  # 0:no clean; 1:clean
      threads = $numCPU  # Number of threads to use
       cutoff = 0  # Mismatch cutoff - ignore SNPs within cutoff length of each other.
CONTRL


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

sub extract_ref_seq {
	my $name=shift;
	#my $list=shift;
	my $dir=shift;
	my $annotion_fh=shift;
	my $out_file = "$dir/$name.fna";
	my @tmpl = `ls -d $EDGE_HOME/database/NCBI_genomes/$name*`;
	chomp @tmpl;
	my @gfiles = `ls -S $tmpl[0]/*gbk $tmpl[0]/*gbff 2>/dev/null`;
	foreach my $gbk (@gfiles){
		chomp $gbk;
		my ($gbk_basename,$gbk_path,$gbk_suffix)=fileparse($gbk,qr/\.[^.]*/);
		my $acc = ($gbk_basename =~ /^(GC[FA]_\d+\.\d)/)? $1 : $gbk_basename;
		system("$EDGE_HOME/scripts/genbank2fasta.pl $gbk >> $out_file");
		my $get = `grep ">" $out_file | head -n 1 | sed -e 's/>//'`;
		my @ref= split /\s+/,$get;
		my $ref_name= join(" ",@ref);
		print $annotion_fh "$name\t".$ref_name."\thttp://www.ncbi.nlm.nih.gov/nuccore/$acc\n" if ($ref_name !~ /plasmid/);
	}
}

sub usage {
     print <<"END";
     Usage: perl $0 [options] -c contig.fasta -p 'reads1.fastq reads2.fastq' -o out_directory [-db Ecoli | -genomes file]
     Version $version
     Input File:
            -c            Contig fasta file
            
            -p            Paired reads in two fastq files and separate by space
            
            -s            Single end reads

            -o            Output Directory (default: SNP_phylogeny)
     
            -db	          Available choices are Ecoli, Yersinia, Francisella, Brucella, Bacillus.
                 OR 
            -genomesList  A comma separated NCBI RefSeq genome names or list in a file (one name per line)
            -genomesFiles A comma separated genome Files
	    -reference    A reference genome name for reads/contigs mapping to 
  
     Options:
            -map          SNPdb name text file. (default: $RealBin/../database/SNPdb/reference.txt) 
            
            -ref_json     SNP reference genome list json file (default: $RealBin/../edge_ui/data/SNP_ref_list.json)
     
            -bwa_id_map   BWA index id map file (default: $RealBin/../database/bwa_index/id_mapping.txt)
            
            -bwa_genome_index  (default: $RealBin/../database/bwa_index/NCBI-Bacteria-Virus.fna)
            
            -kingdom      Bacteria or Virus or Eukaryote (default: Bacteria). This option will affect SNP calling
            
            -tree         Tree Maker: FastTree or RAxML  (default: FastTree)

	    -modeltest    [bool] perform modeltest when building a tree using RAxML
            
            -bootstrap    [bool] Run bootstrapping
            
            -bootstrap_n  Number of bootstraps to run (default 100)
            
            -PosSelect    Molecular evolutionary analysis: PAML or HyPhy (default: no ME analysis)

            -nanopore     [bool] use minimap2 as aligner 
            
            -n            Name of the project

            -cpu          number of CPUs (default: 4)
 
            -version      Print verison

END
exit;
}
