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
my $proxy = $ENV{HTTP_PROXY} || $ENV{http_proxy};
$proxy = "--proxy $proxy" if ($proxy);
my $gzip;

my $res=GetOptions(
    'local|l'                  => \$local_mode,
    'outdir|d=s'               => \$OUTDIR,
    'concat|c'                 => \$CONCAT,
    'concat-prefix|p=s'        => \$CONCAT_PREFIX,
    'platform-restrict|r=s'    => \$PLAT_R,
    'clean|n'                => \$CLEAN,
    'help|?'        => sub{usage()}
) || &usage();

if ( !scalar @ARGV ) { &usage(); }

## init temp directory
`mkdir -p $OUTDIR/sra2fastq_temp`;
`rm -rf $OUTDIR/sra2fastq_temp/*`;

## Main
my @srr;

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
		
		my $path = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$sra";
		
		print STDERR "The input SRA# is $sra_type.\n";
		print STDERR "Trying to Retrieve SRR numbers belong to $sra...\n";

		my $list = `/usr/bin/curl $proxy \"$path\" 2>/dev/null`;

		my @sra_srr_list = $list =~ /\n(SRR\d+),/mg;

		print STDERR scalar(@sra_srr_list)," number found.\n";
	
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

print STDERR "\nDownloading...\n";
my @sra_files = getSRAFastq(@srr);
foreach my $sra_file (@sra_files){
	die "ERROR: SRA file '$sra_file' doesn't exist.\n" unless -e $sra_file;
}

`mv $OUTDIR/sra2fastq_temp/*.gz $OUTDIR`;

if( $CLEAN ){
	print STDERR "\nCleaning up...";
	`rm -rf $OUTDIR/sra2fastq_temp`;
	print STDERR "Done.\n";
}

## Subroutines ########################################################################

sub getSRAFastq {
	my @sra_acc = @_;
	my @sra_files;
	my $sralist = join ",", @sra_acc;
	my $url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=$sralist&format=fastq";
	my $agent = "Mozilla/5.0";

	system("/usr/bin/curl $proxy -A \"$agent\" -o $OUTDIR/sra2fastq_temp/$sralist.fastq.gz -L \"$url\"");
	die "ERROR: Failed to download SRA file from ftp://ftp-trace.ncbi.nih.gov. Please check http://www.ncbi.nlm.nih.gov/sra/?term=$sralist\n" unless -s "$OUTDIR/sra2fastq_temp/$sralist.fastq.gz";

	#Paired-end flag
	my $peflag = 0;
	my $sraHead = `gzip -dc $OUTDIR/sra2fastq_temp/$sralist.fastq.gz | head -5`;
	$peflag = 1 if $sraHead =~ /\.1\.1 .*\n.*\n.*\n.*\n.*\.1\.2 /m; 

	if( $peflag ){
		print STDERR "Paired-end reads found. Deinterleaving...";
		my $di_flag = system("gzip -dc $OUTDIR/sra2fastq_temp/$sralist.fastq.gz | deinterleave_fastq.sh $OUTDIR/sra2fastq_temp/$sralist.1.fastq.gz $OUTDIR/sra2fastq_temp/$sralist.2.fastq.gz compress");
		die "ERROR: Failed to deinterleave FASTQ file.\n" if $di_flag > 0;
		print STDERR "Done.\n";
		`rm -f $OUTDIR/sra2fastq_temp/$sralist.fastq.gz`;

		return ("$OUTDIR/sra2fastq_temp/$sralist.1.fastq.gz","$OUTDIR/sra2fastq_temp/$sralist.2.fastq.gz");
	}
	else{
		return ("$OUTDIR/sra2fastq_temp/$sralist.fastq.gz");
	}
}

sub checkRequiredExec {
	die "ERROR: 'gzip' not found.\n" unless `which gzip`;
	die "ERROR: 'curl' not found.\n" unless `which curl`;
	die "ERROR: 'deinterleave_fastq.sh' not found.\n" unless `which curl`;
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
    --outdir|d             Output directory
    --clean                clean up temp directory
    --help/h/?             display this help          

__END__
exit();
}
