#!/usr/bin/env perl

use strict;
use warnings;
use Time::Piece;
use Cwd;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
my $workingDir = Cwd::getcwd();

$ENV{PATH} = "$RealBin:$RealBin/../bin/:$RealBin/../thirdParty/Anaconda2/bin/:$ENV{PATH}";
$ENV{EDGE_HOME} = $RealBin;

my $version="0.1";

my $starttime = localtime;
my $OPSYS = $^O;

msg("** Local time is $starttime");
msg("** You are", $ENV{USER} || 'not telling me who you are!');
msg("** Operating system is $OPSYS");

&check_tools;

sub check_tools{

	my %tools = (
		'perl' => {
			GETVER  => "perl -v | grep 'version'",
			REGEXP  => qr/\(v(\d+\.\d+)/,
			LICENSE => "GPL"
		},
		'R' => {
			GETVER  => "R --version | grep 'version'",
			REGEXP  => qr/R version\s+(\d+\.\d+\.\d+)/,
			CITATION  => "R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/",
			LICENSE => "GPLv2"
		},
		'samtools' => {
			GETVER  => "samtools 2>\&1| grep 'Version'",
			REGEXP  => qr/Version:\s+(\d+\.\d+\.?\d*)/,
			CITATION  => "Li, H., et al. (2009) The Sequence Alignment/Map format and SAMtools, Bioinformatics, 25, 2078-2079.",
			LICENSE => "MIT"
		},
		'bcftools' => {
			GETVER  => "bcftools 2>\&1| grep 'Version'",
			REGEXP  => qr/Version:\s+(\d+\.\d+\.?\d*)/,
			CITATION  => "Li H (2011), A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data, Bioinformatics 27(21) 2987-93.",
			LICENSE => "MIT"
		},
		'BWA' => {
			GETVER  => "bwa 2>\&1| grep 'Version'",
			REGEXP  => qr/Version:\s+(\d+\.\d+\.\d+)/,
			CITATION  => "Li, H. and Durbin, R. (2009) Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics, 25, 1754-1760.",
			LICENSE => "GPLv3"
                },
		'IDBA-UD' => {
			GETVER => "idba_ud",
			REGEXP => "",
			CITATION => "Peng, Y., et al. (2012) IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth, Bioinformatics, 28, 1420-1428.",
			LICENSE => "GPLv2",
			VERSION => "1.1.1"
		},
		'SPAdes' => {
			GETVER => "spades.py -v 2>\&1",
			REGEXP => qr/SPAdes\s+v(\d+\.\d+\.\d+)/,
			CITATION => "Nurk, Bankevich et al. (2013) Assembling single-cell genomes and mini-metagenomes from chimeric MDA products. J Comput Biol. 2013 Oct;20(10):714-37",
			LICENSE => "GPLv2"
		},
		'MEGAHIT' => {
			GETVER => "megahit -v 2>\&1",
			REGEXP => qr/MEGAHIT\s+v(\d+\.\d+\.\d+)/,
			CITATION => "Li D. et al. (2015) MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics. 2015 May 15;31(10):1674-6",
			LICENSE => "GPLv3"
		},
		'RATT' => {
			GETVER => "../thirdParty/RATT/start.ratt.sh",
			REGEXP => "",
			CITATION => "Otto, T.D., et al. (2011) RATT: Rapid Annotation Transfer Tool, Nucleic acids research, 39, e57",
			LICENSE => ""
		},
		'Prokka' => {
			GETVER => "prokka -v 2>\&1",
			REGEXP => qr/prokka\s+(\d+\.\d+)/,
			CITATION => "Seemann, T. (2014) Prokka: rapid prokaryotic genome annotation, Bioinformatics, 30,2068-2069.",
			LICENSE => "GPLv2"
		},
		'tRNAscan' => {
			GETVER => "tRNAscan-SE 2>\&1|grep  tRNAscan-SE ",
			REGEXP => qr/tRNAscan-SE\s+(\d+\.\d+\.\d+)/,
			CITATION => "Lowe, T.M. and Eddy, S.R. (1997) tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence, Nucleic acids research, 25, 955-964.",
			LICENSE => "GPLv2"
		},
		'Barrnap' => {
			GETVER => "barrnap -v 2>\&1",
			REGEXP => qr/barrnap\s+(\d+\.\d+\.\d+)/,
			CITATION => "",
			LICENSE => "GPLv3"
		},
		'BLAST+' => {
			GETVER => "blastn -version 2>\&1 | grep Package",
			REGEXP => qr/blast\s+(\d+\.\d+\.\d+)/,
			CITATION => "Camacho, C., et al. (2009) BLAST+: architecture and applications, BMC bioinformatics, 10, 421.",
			LICENSE => "Public domain"
		},
		'blastall' => {
			GETVER => "blastall | grep blastall",
			REGEXP => qr/blastall\s+(\d+\.\d+\.\d+)/,
			CITATION => "Altschul, S.F., et al. (1990) Basic local alignment search tool, Journal of molecular biology, 215, 403-410.",
			LICENSE => "Public domain"
		},
		'Phage_Finder' => {
			GETVER => "perl ../thirdParty/phage_finder_v2.1/bin/Phage_Finder_v2.1.pl -V 2>\&1",
			REGEXP => qr/(\d+\.\d+)/,
			CITATION => "Fouts, D.E. (2006) Phage_Finder: automated identification and classification of prophage regions in complete bacterial genome sequences, Nucleic acids research, 34, 5839-5851.",
			LICENSE => "GPLv3"
		},
		'Glimmer' => {
			GETVER => "glimmer3",
			REGEXP => "",
			CITATION => "Delcher, A.L., et al. (2007) Identifying bacterial genes and endosymbiont DNA with Glimmer, Bioinformatics, 23, 673-679.",
			LICENSE => "Artistic License",
			VERSION => "302b"
		},
		'ARAGORN' => {
			GETVER => "aragorn -h | grep ARAGORN",
			REGEXP => qr/ARAGORN\s+v(\d+\.\d+\.\d+)/,
			CITATION => "Laslett, D. and Canback, B. (2004) ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences, Nucleic acids research, 32, 11-16.",
			LICENSE => ""
		},
		'Prodigal' => {
			GETVER => "prodigal -v 2>\&1",
			REGEXP => qr/Prodigal\s+V(\d+\.\d+)/,
			CITATION => "Hyatt, D., et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification, BMC bioinformatics, 11, 119.",
			LICENSE => "GPLv3"
		},
		'tbl2asn' => {
			GETVER => "tbl2asn.orig --help | grep tbl2asn",
			REGEXP => qr/tbl2asn\s+(\d+\.\d+)/,
			CITATION => "",
			LICENSE => ""
		},
		'HMMER3' => {
			GETVER => "hmmalign -h | grep HMMER",
			REGEXP => qr/HMMER\s+(\d+\.\S+)/,
			CITATION => "Eddy, S.R. (2011) Accelerated Profile HMM Searches, PLoS computational biology, 7, e1002195",
			LICENSE => "GPLv3"
		},
		'Infernal' => {
			GETVER => "cmalign -h | grep INFERNAL",
			REGEXP => qr/INFERNAL\s+(\d+\.\S+)/,
			CITATION => "Nawrocki, E.P. and Eddy, S.R. (2013) Infernal 1.1: 100-fold faster RNA homology searches, Bioinformatics, 29, 2933-2935.",
			LICENSE => "GPLv3"
		},
		'Bowtie2' => {
			GETVER => "bowtie2 --version | grep bowtie",
			REGEXP => qr/version\s+(\d+\.\d+\.\d+)/,
			CITATION => "Langmead, B. and Salzberg, S.L. (2012) Fast gapped-read alignment with Bowtie 2, Nature methods, 9, 357-359.",
			LICENSE => "GPLv3"
		},
		'MUMmer3' => {
			GETVER => "nucmer -V",
			REGEXP => "",
			CITATION => "Kurtz, S., et al. (2004) Versatile and open software for comparing large genomes, Genome biology, 5, R12.",
			LICENSE => "GPLv3",
			VERSION => "3.23"
		},
		'RAPSearch2' => {
			GETVER => "rapsearch2 2>\&1| grep rapsearch",
			REGEXP => qr/rapsearch2\s+v(\d+\.\d+)/,
			CITATION => "Zhao et al. (2012) RAPSearch2: a fast and memory-efficient protein similarity search tool for next-generation sequencing data. Bioinformatics. 2012 Jan 1;28(1):125-6",
			LICENSE => "GPL"
		},
		'Kraken' => {
			GETVER => "kraken -v 2>\&1",
			REGEXP => qr/Kraken version\s+(\S+)/,
			CITATION => "Wood, D.E. and Salzberg, S.L. (2014) Kraken: ultrafast metagenomic sequence classification using exact alignments, Genome biology, 15, R46.",
			LICENSE => "GPLv3"
		},
		'Metaphlan' => {
			GETVER => "metaphlan2.py -v 2>\&1",
			REGEXP => qr/version\s+(\d+\.\d+\.\d+)/,
			CITATION => "Segata, N., et al. (2012) Metagenomic microbial community profiling using unique clade-specific marker genes, Nature methods, 9, 811-814.",
			LICENSE => "Artistic Liense"
		},
		'GOTTCHA' => {
			GETVER => "gottcha.pl -h | grep VERSION",
			REGEXP => qr/VERSION:\s+(\S+)/,
			CITATION => "Tracey Allen K. Freitas, Po-E Li, Matthew B. Scholz, Patrick S. G. Chain (2015) Accurate Metagenome characterization using a hierarchical suite of unique signatures. Nucleic Acids Research (DOI: 10.1093/nar/gkv180)",
			LICENSE => "GPLv3"
		},
		'GOTTCHA2' => {
			GETVER => "../thirdParty/gottcha2/gottcha.py -h | grep VERSION",
			REGEXP => qr/VERSION:\s+(\S+)/,
			CITATION => "",
			LICENSE => "GPL"
		},
		'FastTree' => {
			GETVER => "FastTreeMP  2>\&1| grep version",
			REGEXP => qr/version\s+(\d+\.\d+\.\d+)/,
			CITATION => "Morgan N. Price, Paramvir S. Dehal, and Adam P. Arkin. 2009. FastTree: Computing Large Minimum Evolution Trees with Profiles instead of a Distance Matrix. Mol Biol Evol (2009) 26 (7): 1641-1650",
			LICENSE => "GPLv2"
		},
		'RAxML' => {
			GETVER => "raxmlHPC-PTHREADS -v | grep version",
			REGEXP => qr/version\s+(\d+\.\d+\.\d+)/,
			CITATION => "Stamatakis,A. 2014. RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30:1312-1313",
			LICENSE => "GPLv2"
		},
		'PhaME' => {
			GETVER => "./PhaME/src/runPhaME.pl",
			REGEXP => "",
			CITATION => "Sanaa Afroz Ahmed, Chien-Chi Lo, Po-E Li, Karen W Davenport, Patrick S.G. Chain. From raw reads to trees: Whole genome SNP phylogenetics across the tree of life. bioRxiv doi: http://dx.doi.org/10.1101/032250",
			LICENSE => "GPLv3",
			VERSION => "1.0"
		},
		'ShortBRED' => {
			GETVER => "../thirdParty/ShortBRED-0.9.4M/shortbred_quantify.py --version",
			REGEXP => "",
			CITATION => "Kaminski J, et al. (2015) High-specificity targeted functional profiling in microbial communities with ShortBRED. PLoS Comput Biol.18;11(12):e1004557.",
			LICENSE => "MIT",
			VERSION => "0.9.4M"
		},
		'RGI (Resistance Gene Identifier)' => {
			GETVER => "rgi -V 2>\&1",
			REGEXP => qr/rgi-(\d+\.\d+\.\d+)/,
			CITATION => "McArthur & Wright. (2015) Bioinformatics of antimicrobial resistance in the age of molecular epidemiology. Current Opinion in Microbiology, 27, 45-50.",
			LICENSE => "Apache Software License"
		},
		'JBrowse' => {
			GETVER => "../edge_ui/JBrowse/bin/prepare-refseqs.pl",
			REGEXP => "",
			CITATION => "Skinner, M.E., et al. (2009) JBrowse: a next-generation genome browser, Genome research, 19, 1630-1638.",
			LICENSE => "Artistic License 2.0/LGPLv.1",
			VERSION => "1.11.6"
		},
		'Bio::Phylo' => {
			GETVER => "",
			REGEXP => "",
			CITATION => "Rutger A Vos, Jason Caravas, Klaas Hartmann, Mark A Jensen and Chase Miller, (2011). Bio::Phylo - phyloinformatic analysis using Perl. BMC Bioinformatics 12:63.",
			LICENSE => "GPLv3",
			VERSION => "0.58"
		},
		'jsPhyloSVG' => {
			GETVER => "",
			REGEXP => "",
			CITATION => "Smits SA, Ouverney CC, (2010) jsPhyloSVG: A Javascript Library for Visualizing Interactive and Vector-Based Phylogenetic Trees on the Web. PLoS ONE 5(8): e12267.",
			LICENSE => "GPL",
			VERSION => "1.55"
		},
		'KronaTools' => {
			GETVER => "ktImportTaxonomy | grep KronaTools ",
			REGEXP => qr/KronaTools\s+(\d+\.\d+)/,
			CITATION => "Ondov, B.D., Bergman, N.H. and Phillippy, A.M. (2011) Interactive metagenomic visualization in a Web browser, BMC bioinformatics, 12, 385.",
			LICENSE => "BSD"
		},
		'BEDTools' => {
			GETVER => "bedtools --version 2>\&1",
			REGEXP => qr/bedtools\s+v(\d+\.\d+\.\d+)/,
			CITATION => "Quinlan, A.R. and Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, 26, 841-842.",
			LICENSE => "GPLv2"
		},
		'MetaComp' => {
			GETVER => "",
			REGEXP => "",
			CITATION => "",
			LICENSE => "",
			VERSION => "1.0.1"
		},
		'GNU_parallel' => {
			GETVER => "parallel --version",
			REGEXP => qr/parallel\s+(\d+)/,
			CITATION => "O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47",
			LICENSE => "GPLv3"
		},
		'tabix' => {
			GETVER => "tabix 2>&1| grep Version",
			REGEXP => qr/Version:\s+(\d+\.\d+\.\d+)/,
			CITATION => "",
			LICENSE => ""
		},
		'Primer3' => {
			GETVER => "primer3_core  -h 2>&1| grep release",
			REGEXP => qr/release\s+(\d+\.\d+\.\d+)/,
			CITATION => "Untergasser, A., et al. (2012) Primer3–new capabilities and interfaces, Nucleic acids research, 40, e115.",
			LICENSE => "GPLv2"
		},
		'FaQCs' => {
			GETVER => "illumina_fastq_QC.pl -v",
			REGEXP => qr/Version:\s+(\d+\.\d+)/,
			CITATION => "Chienchi Lo, PatrickS.G. Chain (2014) Rapid evaluation and Quality Control of Next Generation Sequencing Data with FaQCs. BMC Bioinformatics. 2014 Nov 19;15",
			LICENSE => "GPLv3"
		},
		'wigToBigWig' => {
			GETVER => "wigToBigWig 2>&1 | grep wigToBigWig",
			REGEXP => qr/wigToBigWig\s+v\s+(\d+)/,
			CITATION => "Kent, W.J., et al. (2010) BigWig and BigBed: enabling browsing of large distributed datasets, Bioinformatics, 26, 2204-2207.",
			LICENSE => ""
		},
		'NCBI-sratoolkit' => {
			GETVER => "fastq-dump -version",
			REGEXP => qr/fastq-dump\s+:\s+(\d+\.\d+\.\d+)/,
			CITATION => "",
			LICENSE => ""
		},
		'ea-utils' => {
			GETVER => "fastq-join",
			REGEXP => "",
			CITATION => "Erik Aronesty (2011) ea-utils : “Command-line tools for processing biological sequencing data”",
			LICENSE => "MIT",
			VERSION => "1.1.2-537"
		},
		'Anaconda2' => {
			GETVER => "",
			REGEXP => "",
			CITATION => "",
			LICENSE => "3-clause BSD",
			VERSION => "4.1.1"
		},
		'Anaconda3' => {
			GETVER => "",
			REGEXP => "",
			CITATION => "",
			LICENSE => "3-clause BSD",
			VERSION => "4.1.1"
		},
		'QIIME' => {
			GETVER => "split_libraries_fastq.py --version",
			REGEXP => qr/(\d+\.\d+\.\d+)/,
			CITATION => "Caporaso et al. (2010) QIIME allows analysis of high-throughput community sequencing data. Nat Methods. 2010 May;7(5):335-6",
			LICENSE => "GPLv2"
		},
		'minimap2' => {
			GETVER  => "minimap2 --version",
			REGEXP  => qr/(\d+\.\d+\.?\d*)-/,
			CITATION  => "Li, H. (2017). Minimap2: fast pairwise alignment for long nucleotide sequences. arXiv:1708.01492",
			LICENSE => "MIT"
		},
		'TargetedNGS' => {
			GETVER => "../thirdParty/targetngs/targetedNGS -v",
			REGEXP => qr/(\d+\.\d+\.\d+)/,
			CITATION => "",
			LICENSE => "GPL"
		},
	);
	msg("\nEDGE TOOLS:");
	foreach my $toolname (sort keys %tools){
		my $tool = $tools{$toolname};
		my ($exe) = $tool->{GETVER} =~ /^(\S+)/;
		my $fp = find_exe($exe) if $tool->{GETVER};
		msg("* $toolname:");
		msg("*** Can't find '$toolname' in your \$PATH ***") if !$fp && $tool->{GETVER};
		if ($fp){
			if($tool->{REGEXP}) {
				my($s) = qx($tool->{GETVER});
				if (defined $s) {
					$s =~ $tool->{REGEXP};
					$version = $1 if defined $1;
					msg("  Version: $version");
				}else{
					msg("  Version: Could not determine version");
				}
			}else{
				#msg("  Version: Could not determine version");
			}
		}
		if($tool->{VERSION}){
			msg("  Version: $tool->{VERSION}");
		}
		if(defined $tool->{CITATION}){
			msg("  Citation: $tool->{CITATION}");
		}
		if(defined $tool->{LICENSE}){
			msg("  License: $tool->{LICENSE}");
		}
	}
}

sub msg {
        my $t = localtime;
        my $line = "@_\n";
        print STDERR $line;
}

sub find_exe {
        my($bin) = shift;
        for my $dir (File::Spec->path) {
                my $exe = File::Spec->catfile($dir, $bin);
                return $exe if -x $exe;
        }
        return;
}
