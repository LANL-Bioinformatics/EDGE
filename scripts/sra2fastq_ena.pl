#!/usr/bin/perl
# sra2fastq.pl
# ver 0.1
# 2014/12/19
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# ver 0.3 (2016/10/28)
# - Switch sequecne source to ERA
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
my $PLAT_R;
my $CLEAN;
my $FLSZ_R;
my $proxy = $ENV{HTTP_PROXY} || $ENV{http_proxy};
$proxy = "--proxy $proxy " if ($proxy);
my $gzip;

my $res=GetOptions(
    'outdir|d=s'               => \$OUTDIR,
    'platform-restrict|r=s'    => \$PLAT_R,
    'filesize-restrict|r=s'    => \$FLSZ_R,
    'clean|n'                => \$CLEAN,
    'help|?'        => sub{usage()}
) || &usage();

if ( !scalar @ARGV ) { &usage(); }

## init temp directory
`rm -rf $OUTDIR/sra2fastq_temp/`;
`mkdir -p $OUTDIR/sra2fastq_temp/merged/`;

## Main
my $enaInfo;
my $sra_type = "ByRun";

foreach my $acc ( @ARGV ){
	die "ERROR: $acc is not a valid SRA/ERA/DRA number." if $acc !~ /^(SR|ER|DR)/;
	
	$sra_type = "ByExp"    if $acc =~ /^(SRX|ERX|DRX)/;
	$sra_type = "BySample" if $acc =~ /^(SRS|ERS|DRS)/;
	$sra_type = "ByStudy"  if $acc =~ /^(SRP|ERP|DRP)/;
	
	print STDERR "$acc -- $sra_type:\n";

	my $path = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$acc&result=read_run&fields=run_accession,instrument_platform,fastq_ftp,fastq_md5,fastq_bytes";
	print STDERR "Retrieving run acc# from $path...\n";
	
	my $result = `/usr/bin/curl $proxy "$path" 2>/dev/null`;

	print STDERR "\n$result\n";

	my @lines = split '\n', $result;

	foreach my $line (@lines){
		next if $line =~ /^run_accession/;
		chomp;

		my @f = split '\t', $line;
		$enaInfo->{$acc}->{$f[0]}->{instrument_platform} = $f[1];

		my @url  = split ';', $f[2];
		my @md5  = split ';', $f[3];
		my @size = split ';', $f[4];

		for( my $i=0; $i<=$#url; $i++ ){
			$enaInfo->{$acc}->{$f[0]}->{fastq_ftp}->{$i}->{url}  = $url[$i];
			$enaInfo->{$acc}->{$f[0]}->{fastq_ftp}->{$i}->{md5}  = $md5[$i];
			$enaInfo->{$acc}->{$f[0]}->{fastq_ftp}->{$i}->{size} = $size[$i];
		}
	}

	if( $result !~ /^run_accession/ ){
		die "ERROR: ENA errors. Please check if $acc is a valid SRA/ERA/DRA number or your internet connection.\n";
	}
	if( ! defined $enaInfo->{$acc} ){
		die "ERROR: No run found for $acc.\n";
	}
}

print STDERR "Downloading...\n";
getSRAFastq($enaInfo);

`mv $OUTDIR/sra2fastq_temp/merged/*.gz $OUTDIR`;

if( $CLEAN ){
	print STDERR "\nCleaning up...";
	`rm -rf $OUTDIR/sra2fastq_temp`;
	print STDERR "Done.\n";
}

## Subroutines ########################################################################

sub getSRAFastq {
	my $enaInfo = shift;
	my $agent = "Mozilla/5.0";

	foreach my $acc (keys %$enaInfo){
		foreach my $run_acc (keys %{$enaInfo->{$acc}}){
			print STDERR "$acc => run_acc# $run_acc...\n";

			if( $PLAT_R && $enaInfo->{$acc}->{$run_acc}->{instrument_platform} !~ /$PLAT_R/i ){
				die "ERROR: $enaInfo->{$acc}->{$run_acc}->{instrument_platform} platfrom detected. Only $PLAT_R is allowed.\n";
			}
			
			foreach my $i ( keys %{$enaInfo->{$acc}->{$run_acc}->{fastq_ftp}} ){
				my $url  = $enaInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{url};
				my $md5  = $enaInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{md5};
				my $size = $enaInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{size};
				
				my ($filename) = $url =~ /\/([^\/]+)$/;
								
				print "Downloading $url...\n";
				system("/usr/bin/curl $proxy -A \"$agent\" -o $OUTDIR/sra2fastq_temp/$filename \"http://$url\"");
				#system("wget --no-proxy -O $OUTDIR/sra2fastq_temp/$filename \"ftp://$url\"");
				
				
				#check downloaded file
				my $filesize = -s "$OUTDIR/sra2fastq_temp/$filename";

				if( !$filesize ){
					die "ERROR: Failed to download $filename from $url.\n";
				}
				if( $filesize != $size ){
					die "ERROR: $OUTDIR/sra2fastq_temp/$filename incompleted/corrupted -- file sizes mismatch.\n";
				}

				#check md5
				my $md5sum = `md5sum $OUTDIR/sra2fastq_temp/$filename`;
				if(  $md5sum !~ /^$md5/ ){
					die "ERROR: $OUTDIR/sra2fastq_temp/$filename corrupted -- md5 checksum mismatch.\n";
				}

				#merging file
				if( $filename =~ /_1/ ){
					`cat $OUTDIR/sra2fastq_temp/$filename >> $OUTDIR/sra2fastq_temp/merged/${acc}.1.fastq.gz`;
				}
				elsif( $filename =~ /_2/ ){
					`cat $OUTDIR/sra2fastq_temp/$filename >> $OUTDIR/sra2fastq_temp/merged/${acc}.2.fastq.gz`;
				}
				else{
					`cat $OUTDIR/sra2fastq_temp/$filename >> $OUTDIR/sra2fastq_temp/merged/${acc}.fastq.gz`;
				}
				`rm $OUTDIR/sra2fastq_temp/$filename`;
			}
		}

		my $acc_filesize = -s "$OUTDIR/sra2fastq_temp/merged/${acc}.1.fastq.gz" + -s "$OUTDIR/sra2fastq_temp/merged/${acc}.2.fastq.gz" + -s "$OUTDIR/sra2fastq_temp/merged/${acc}.fastq.gz";
		if( $FLSZ_R && $FLSZ_R < $acc_filesize/1024/1024 ){
			die "ERROR: downloaded file size is too large.\n";
		}
		else{
			print STDERR "Succesfully downloaded ${acc}.\n"
		}
	}
}

sub checkRequiredExec {
	die "ERROR: 'gzip' not found.\n" unless `which gzip`;
	die "ERROR: 'curl' not found.\n" unless `which curl`;
}

sub usage {
print <<__END__;

[DESCRIPTION]
    A script retrieve FASTQ files from EBI ENA database using `curl`. Input 
accessions support studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), 
samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions 
(SRA*/ERA*/DRA*).

[USAGE]
    $0 [OPTIONS] <Accession> (<Accession 2> <Accession 3>...)

[OPTIONS]
    --outdir|d             Output directory
    --clean                clean up temp directory
    --platform-restrict    Only allow a specific platform
    --filesize-restrict    (in MB) Only allow downloaded file smaller than a specific size.
    --help/h/?             display this help          

__END__
exit();
}
