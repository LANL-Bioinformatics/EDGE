#!/usr/bin/perl
# sra2fastq.pl
# ver 0.1
# 2014/12/19
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# ver 0.2 (2015/01/06)
# - Input SRA accessions support studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions (SRA*/ERA*/DRA*)
# - Provide "--platform-restrict" option to limit the platform of SRAs
# - Provide "--concat" option to concatenate multiple FASTQ files into a singal (single-end) or two (paired-end) files
# - Remove dependency of File::Which

use strict;
use Getopt::Long;
use FindBin qw($Bin);
#use File::Which;

$|=1;
$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";

checkRequiredExec();

my $local_mode = 0;
my $OUTDIR = ".";
my $CONCAT;
my $PLAT_R;
my $CLEAN;
my $CONCAT_PREFIX = "output";
my $gzip;

my $res=GetOptions(
    'local|l'                  => \$local_mode,
    'outdir|d=s'               => \$OUTDIR,
    'concat|c'                 => \$CONCAT,
    'concat-prefix|p=s'        => \$CONCAT_PREFIX,
    'platform-restrict|r=s'    => \$PLAT_R,
    'clean|n'                => \$CLEAN,
    'gzip|z'		       => \$gzip,
    'help|?'        => sub{usage()}
) || &usage();

if ( !scalar @ARGV ) { &usage(); }

## init temp directory
`mkdir -p $OUTDIR/sra2fastq_temp`;
`rm -rf $OUTDIR/sra2fastq_temp/*`;

## Main
my $base = "ftp://ftp-trace.ncbi.nih.gov";
my @srr;

if( $local_mode ){
	@srr = @ARGV;
}
else{
	foreach my $sra ( @ARGV ){
		my ($sra_acc_first6, $sra_acc_prefix) = $sra =~ /^((\w{3})\d{3})/;
		my $sra_type = "ByRun";
		if( $sra_acc_prefix =~ /(SRR|ERR|DRR)/ ){
			push @srr, $sra;
		}
		else{
			$sra_type = "ByExp"    if $sra_acc_prefix =~ /(SRX|ERX|DRX)/;
			$sra_type = "BySample" if $sra_acc_prefix =~ /(SRS|ERS|DRS)/;
			$sra_type = "ByStudy"  if $sra_acc_prefix =~ /(SRP|ERP|DRP)/;
			my $path = "$base/sra/sra-instant/reads/$sra_type/sra/$sra_acc_prefix/$sra_acc_first6/$sra/";
			
			print STDERR "Checking $sra...\n";
			my $list = `/usr/bin/curl $path 2>/dev/null`;
			#my @sra_srr_list = $list =~ /dr-xr-xr-x \d+ ftp\s+anonymous\s+\d+\s+\w+\s+\d+\s+\d+\s+(\w{3}\d+)/g;
			my @sra_srr_list = $list =~ /ftp\s+anonymous.*(SRR\d+)\n/g;
		
			if( @sra_srr_list ){
				print STDERR "Found SRA runs $sra_type: $sra\n";
				print STDERR "  - ".join("\n  - ",@sra_srr_list)."\n";
				push @srr, @sra_srr_list;
			}
			else{
				print STDERR "SRA runs not found.\n\n";
			}
		}
	}
}

my $flag;
foreach my $sra ( @srr ){
	print STDERR "\nProcessing $sra...\n";
	my $sra_file = $sra;
	$sra_file = getSRA($sra) unless( $local_mode );
	die "ERROR: SRA file '$sra_file' doesn't exist.\n" unless -e $sra_file;

	my $info = getSRA_info($sra_file);
	print "\n$sra Metadata\n";
	print "==========================================================\n";
	print "$info\n";

	if( $PLAT_R && $info !~ /$PLAT_R/i ){
		print "\nWARNING: Platform doesn't match. Skipping $sra.\n\n";
		next;
	}

	print STDERR "Converting $sra to FASTQ...\n";
	SRA_to_FASTQ($sra_file);
	print STDERR "Done converting $sra.\n";
	$flag=1;# has something coverted..
}

if( $flag && $CONCAT ){
	print STDERR "\nConcatenating FASTQ files...";
	my $suffix = ($gzip)?".gz":"";
	`cat $OUTDIR/sra2fastq_temp/*_1.fastq$suffix > $OUTDIR/${CONCAT_PREFIX}_1.fastq$suffix 2>/dev/null`;
	`cat $OUTDIR/sra2fastq_temp/*_2.fastq$suffix > $OUTDIR/${CONCAT_PREFIX}_2.fastq$suffix 2>/dev/null`;
	`cat $OUTDIR/sra2fastq_temp/*[:digit:][:digit:].fastq$suffix > $OUTDIR/${CONCAT_PREFIX}.fastq$suffix 2>/dev/null`;
	`rm -rf $OUTDIR/${CONCAT_PREFIX}_1.fastq 2>/dev/null` unless -s "$OUTDIR/${CONCAT_PREFIX}_1.fastq";
	`rm -rf $OUTDIR/${CONCAT_PREFIX}_2.fastq 2>/dev/null` unless -s "$OUTDIR/${CONCAT_PREFIX}_2.fastq";
	`rm -rf $OUTDIR/${CONCAT_PREFIX}.fastq 2>/dev/null` unless -s "$OUTDIR/${CONCAT_PREFIX}.fastq";
	print STDERR "Done.\n\n";
}
elsif($flag){
	($gzip)?`mv $OUTDIR/sra2fastq_temp/*.gz $OUTDIR`:
		`mv $OUTDIR/sra2fastq_temp/*.fastq $OUTDIR`;
}

if( $CLEAN ){
	print STDERR "\nCleaning up...";
	`rm -rf $OUTDIR/sra2fastq_temp`;
	print STDERR "Done.\n";
}

## Subroutines ########################################################################

sub getSRA {
	my ($sra_acc) = @_;
	my ($sra_acc_first6, $sra_acc_prefix) = $sra_acc =~ /^((\w{3})\d{3})/;
	my $sra_type;
	$sra_type = "ByRun"    if $sra_acc_prefix =~ /(SRR|ERR|DRR)/;
	die "ERROR: Invalid SRA prefix. Please input accessions in SRA runs.\n" unless $sra_type;

	print STDERR "Downloading $sra_acc...\n";
	my $url = "$base/sra/sra-instant/reads/$sra_type/sra/$sra_acc_prefix/$sra_acc_first6/$sra_acc/$sra_acc.sra";
	system("/usr/bin/curl -o $OUTDIR/sra2fastq_temp/$sra_acc.sra $url");
	die "ERROR: Failed to download SRA file from ftp://ftp-trace.ncbi.nih.gov. Please check http://www.ncbi.nlm.nih.gov/sra/?term=$sra_acc\n" unless -s "$OUTDIR/sra2fastq_temp/$sra_acc.sra";

	return "$OUTDIR/sra2fastq_temp/$sra_acc.sra";
}

sub getSRA_info {
	my ($sra_file) = @_;
	my $info = `vdb-dump --info $sra_file`;
	return $info;
}

sub SRA_to_FASTQ {
	my ($sra_file) = @_;
	my $gzip_flag = ($gzip)? "-gzip":""; 
	my $ec = system("fastq-dump --split-files $gzip_flag --outdir '$OUTDIR/sra2fastq_temp' $sra_file");
	die "ERROR: Failed to dump FASTQ from SRA file.\n" if $ec > 0;
}

sub checkRequiredExec {
	die "ERROR: 'curl' not found.\n" unless `which curl`;
	die "ERROR: 'fastq-dump' not found. Please install SRA-Toolkit properly and add executables to the path.\n" unless `which fastq-dump`;
}

sub usage {
print <<__END__;

[DESCRIPTION]
    A script retrieve SRA files from NCBI SRA database using `curl` and convert SRA to
FASTQ format with SRA-toolkit. Input SRA accessions support studies (SRP*/ERP*/DRP*), 
experiments (SRX*/ERX*/DRX*), samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or 
submissions (SRA*/ERA*/DRA*).

[USAGE]
    $0 [OPTIONS] <SRA Accession> (<SRA Accession 2> <SRA Accession 3>...)
    
    Or
    
    $0 [OPTIONS] --local </path/to/SRA_File.sra> (<SRA File2>,<SRA File3>...)

[OPTIONS]
    --local|l              Convert local SRA files to FASTQ files
    --outdir|d             Output directory
    --concat|c             Concatenate FASTQ files to a single (SE) or two (PE) files
    --concat-prefix|p      The filename PREFIX of concatenated FASTQ files.
                           Default is "output".
    --platform-restrict|r  Only convert a specific platform of SRAs (case insensitive). The
                           shorten platform name is allowed, e.g. "illu" for "Illumina".
    --clean                clean up temp directory
    --gzip                 compressed output fastq file
    --help/h/?             display this help          

__END__
exit();
}
