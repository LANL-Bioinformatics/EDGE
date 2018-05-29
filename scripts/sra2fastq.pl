#!/usr/bin/perl
# sra2fastq.pl
# ver 0.5
# 2014/12/19
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# ver 0.5 (2016/11/02)
# - support downloading sequences from NCBI-SRA, EBI-ENA and DDBJ.
# - auto switch download sites.
# ver 0.3 (2016/10/28)
# - Switch sequecne source to ERA
# ver 0.2 (2015/01/06)
# - Input SRA accessions support studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions (SRA*/ERA*/DRA*)
# - Provide "--platform-restrict" option to limit the platform of SRAs
# - Provide "--concat" option to concatenate multiple FASTQ files into a singal (single-end) or two (paired-end) files
# - Remove dependency of File::Which
# ver 0.1
# - Initial release

use strict;
use Getopt::Long;
use FindBin qw($Bin);

$|=1;
$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";

my $local_mode = 0;
my $OUTDIR = ".";
my $PLAT_R;
my $CLEAN;
my $FLSZ_R;
my $RUN_R;
my $Download_tool='curl';
my $user_proxy;
my $no_proxy;
my $http_proxy = $ENV{HTTP_PROXY} || $ENV{http_proxy};
my $ftp_proxy = $ENV{FTP_PROXY} || $ENV{ftp_proxy};
$http_proxy = "--proxy \'$http_proxy\' " if ($http_proxy);
$ftp_proxy = "--proxy \'$ftp_proxy\' " if ($ftp_proxy);
my $gzip;

checkRequiredExec();

my $res=GetOptions(
    'outdir|d=s'             => \$OUTDIR,
    'platform-restrict|pr=s' => \$PLAT_R,
    'filesize-restrict|fr=s' => \$FLSZ_R,
    'runs-restrict|r|rr=s'   => \$RUN_R,
    'download-interface=s'   => \$Download_tool,
    'proxy'                  => \$user_proxy,
    'no_proxy'               => \$no_proxy,
    'clean|n'                => \$CLEAN,
    'help|?'        => sub{usage()}
) || &usage();

if ( !scalar @ARGV ) { &usage(); }

$http_proxy ="" if $no_proxy;
$ftp_proxy ="" if $no_proxy;
$http_proxy = "--proxy \'$user_proxy\' " if ($user_proxy);
$ftp_proxy = "--proxy \'$user_proxy\' " if ($user_proxy);

my $curl = ($Download_tool =~ /wget/)?  "wget -v -U \"Mozilla/5.0\" ": "curl -k -A \"Mozilla/5.0\" -L";

## init temp directory
`rm -rf $OUTDIR/sra2fastq_temp/`;
`mkdir -p $OUTDIR/sra2fastq_temp/merged/`;

## Main
my $readInfo;
my $total_runs=0;
my $total_size=0;
my ($dl_status, $dl_runs, $dl_size);

foreach my $acc ( @ARGV ){
	die "ERROR: $acc is not a valid SRA/ERA/DRA number." if $acc !~ /^(SR|ER|DR)/;

	my $sra_type = "ByRun";
	$sra_type = "ByExp"    if $acc =~ /^(SRX|ERX|DRX)/;
	$sra_type = "BySample" if $acc =~ /^(SRS|ERS|DRS)/;
	$sra_type = "ByStudy"  if $acc =~ /^(SRP|ERP|DRP)/;
	
	print STDERR "Processing $acc ($sra_type)...\n";

	#gathering reads information from NCBI-SRA / EBI-ENA
	$readInfo = getReadInfo( $acc, $readInfo); #retrieve runs infor for each accession

	if( !defined $readInfo->{$acc} ){
		die "ERROR: No sequence found. Please check if $acc is a valid SRA/ERA/DRA number or your internet connection.\n";
	}

	foreach my $run_acc (keys %{$readInfo->{$acc}})
	{
		#check if runs exceed the limit
		$total_runs++;
		if( $RUN_R && $RUN_R < $total_runs ){
			die "ERROR: Run(s) exceed the limit ($RUN_R MB).\n";
		}

		#check platform
		my $platform = $readInfo->{$acc}->{$run_acc}->{platform};
		if( $PLAT_R && $platform !~ /$PLAT_R/i ){
			warn "WARN: $readInfo->{$acc}->{$run_acc}->{platform} platfrom detected. Only $PLAT_R is allowed.\n";
			next;
		}

		#download FASTQ
		$dl_status = getSraFastqToolkits( $readInfo->{$acc}->{$run_acc}, $run_acc ); 
		$dl_status = getSraFastq( $readInfo->{$acc}->{$run_acc}, $run_acc ) if $dl_status eq 'failed';
		$dl_status = getDdbjFastq( $readInfo->{$acc}->{$run_acc}, $run_acc ) if $dl_status eq 'failed';
		$dl_status = getEnaFastq( $readInfo->{$acc}->{$run_acc}, $run_acc ) if $dl_status eq 'failed';
		die "ERROR: Please check your internet connection.\n" if $dl_status eq 'failed';

		#merging fastqs in multiple runs into a single file
		if( -e "$OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz" ){
			`cat $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz >> $OUTDIR/sra2fastq_temp/merged/${acc}.1.fastq.gz`;
			$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz";
			`rm -rf ${run_acc}_1.fastq.gz`;
		}
		if( -e "$OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz" ){
			`cat $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz >> $OUTDIR/sra2fastq_temp/merged/${acc}.2.fastq.gz`;
			$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz";
			`rm -rf ${run_acc}_2.fastq.gz`;
		}

		if( -e "$OUTDIR/sra2fastq_temp/${run_acc}.fastq.gz" ){
			`cat $OUTDIR/sra2fastq_temp/${run_acc}.fastq.gz >> $OUTDIR/sra2fastq_temp/merged/${acc}.fastq.gz`;
			$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}.fastq.gz";
			`rm -rf ${run_acc}.fastq.gz`;
		}
		

		if( $FLSZ_R && $FLSZ_R < $total_size/1024/1024 ){
			die "ERROR: downloaded file size exceed limitation ($FLSZ_R MB).\n";
		}
		else{
			print STDERR "Succesfully downloaded $run_acc.\n";
		}
	}
	
	print STDERR "Finished downloading acc# $acc.\n";
}

`mv $OUTDIR/sra2fastq_temp/merged/*.gz $OUTDIR`;

if( $CLEAN ){
	print STDERR "\nCleaning up...";
	`rm -rf $OUTDIR/sra2fastq_temp`;
	print STDERR "Done.\n";
}

## Subroutines ########################################################################

sub getSraFastq {
	my ($info, $run_acc) = @_;

	print STDERR "Retrieving FASTQ for $run_acc from NCBI SRA (online converting)...\n";

	my $platform = $info->{platform};
	my $library  = $info->{library};

	my $url  = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=$run_acc&format=fastq";

	print STDERR "Downloading $url...\n";
	my $cmd = ($Download_tool =~ /wget/)?"$curl -O $OUTDIR/sra2fastq_temp/$run_acc.fastq.gz \"$url\"":"$curl $http_proxy -o $OUTDIR/sra2fastq_temp/$run_acc.fastq.gz \"$url\"";
	my $ec = system("$cmd");

	if( $ec > 0 ){
		print STDERR "Failed to download SRA file from $url.\n";
		return "failed";
	}
	#Deinterleaving if paired-end reads
	if( $platform =~ /illu/i && $library =~ /pair/i ){
		print STDERR "Paired-end reads found. Deinterleaving...";
		my $di_flag = system("gzip -dc $OUTDIR/sra2fastq_temp/$run_acc.fastq.gz | deinterleave_fastq.sh $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz compress");
		return "failed" if $di_flag > 0;
		print STDERR "Done.\n";

		`rm -f $OUTDIR/sra2fastq_temp/$run_acc.fastq.gz`;
	}
	
	return "success";
}

sub getSraFastqToolkits {
	my ($info, $run_acc) = @_;

	print STDERR "Retrieving FASTQ for $run_acc with NCBI SRA Toolkit...\n";

	my $platform = $info->{platform};
	my $url      = $info->{url};
	my $filename = $run_acc;

	print STDERR "Downloading $url...\n";
	my $cmd = ($Download_tool =~ /wget/)? "$curl -O $OUTDIR/sra2fastq_temp/$filename \"$url\"":
					"$curl $http_proxy -o $OUTDIR/sra2fastq_temp/$filename \"$url\"";
	my $ec = system("$cmd");
	
	if( $ec > 0 ){
		print STDERR "Failed to download SRA file from $url.\n";
		return "failed";
	}
	print STDERR "Done.\n";
	
	#check downloaded file
	my $filesize = -s "$OUTDIR/sra2fastq_temp/$filename";

	if( !$filesize ){
		print STDERR "Failed to download SRA file from $url.\n";
		return "failed";
	}

	#dump fastq from SRA file
	my $options = "-gzip ";
	$options .= "--split-files " if $platform =~ /illu/i;
	$options .= "--split-files -B " if $platform =~ /solid/;
	print STDERR "Running fastq-dump with options $options...\n";
	$ec = system("fastq-dump $options --outdir '$OUTDIR/sra2fastq_temp' $OUTDIR/sra2fastq_temp/$filename 2>/dev/null");
	
	if( $ec > 0 ){
		print STDERR "Failed to run fastq-dump from $OUTDIR/sra2fastq_temp/$filename.\n";
		return "failed";
	}
	print STDERR "Done.\n";

	#clean up temp files
	`rm $OUTDIR/sra2fastq_temp/$filename`;
	return "success";
}

sub getDdbjFastq {
	my ($info, $run_acc) = @_;

	print STDERR "Retrieving FASTQ for $run_acc from DDBJ...\n";

	my $platform = $info->{platform};
	my $library  = $info->{library};
	my $exp_acc  = $info->{exp_acc};
	my $sub_acc  = $info->{sub_acc};
	my $platform = $info->{platform};
	my ($sra_acc_first6) = $sub_acc =~ /^(\w{3}\d{3})/;
	my $ec;
	my $cmd;
	my $MinSize=10000;
	
	my $url = "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/$sra_acc_first6/$sub_acc/$exp_acc";
	
	if( $platform =~ /illu/i && $library =~ /pair/i ){
		print STDERR "Downloading $url/${run_acc}_1.fastq.bz2...\n";
		$cmd = ($Download_tool =~ /wget/)? "$curl -O $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.bz2 \"$url/${run_acc}_1.fastq.bz2\"":"$curl $ftp_proxy -o $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.bz2 \"$url/${run_acc}_1.fastq.bz2\"";
		$ec = system("$cmd");
		print STDERR "finished.\n";
		print STDERR "Downloading $url/${run_acc}_2.fastq.bz2...\n";
		$cmd = ($Download_tool =~ /wget/)? "$curl -O $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.bz2 \"$url/${run_acc}_2.fastq.bz2\"": "$curl $ftp_proxy -o $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.bz2 \"$url/${run_acc}_2.fastq.bz2\"";
		$ec = system("$cmd");
		print STDERR "finished.\n";
	}
					
	print STDERR "Downloading $url/$run_acc.fastq.bz2...\n";
	$ec = system("$curl $ftp_proxy -o $OUTDIR/sra2fastq_temp/$run_acc.fastq.bz2 \"$url/$run_acc.fastq.bz2\"");
	print STDERR "finished.\n";

	my $total_size = 0;
	if( -s "$OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.bz2" > $MinSize){
		print STDERR "Convering bz2 to gz...\n";
		$ec = system("bunzip2 -c < $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.bz2 | gzip -c > $OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz");
		if( $ec > 0 ){
			print STDERR "failed to convert bz2 to gz.\n";
			return "failed";
		}
		$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.gz";
		print STDERR "Done.\n";
	}else{
		unlink "$OUTDIR/sra2fastq_temp/${run_acc}_1.fastq.bz2";
	}
	if( -s "$OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.bz2" > $MinSize){
		print STDERR "Convering bz2 to gz...\n";
		$ec = system("bunzip2 -c < $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.bz2 | gzip -c > $OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz");
		if( $ec > 0 ){
			print STDERR "failed to convert bz2 to gz.\n";
			return "failed";
		}
		$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.gz";
		print STDERR "Done.\n";
	}else{
		unlink "$OUTDIR/sra2fastq_temp/${run_acc}_2.fastq.bz2";
	}
	if( -s "$OUTDIR/sra2fastq_temp/${run_acc}.fastq.bz2" > $MinSize){
		print STDERR "Convering bz2 to gz...\n";
		$ec = system("bunzip2 -c < $OUTDIR/sra2fastq_temp/${run_acc}.fastq.bz2 | gzip -c > $OUTDIR/sra2fastq_temp/${run_acc}.fastq.gz");
		if( $ec > 0 ){
			print STDERR "failed to convert bz2 to gz.\n";
			return "failed";
		}
		$total_size += -s "$OUTDIR/sra2fastq_temp/${run_acc}.fastq.gz";
		print STDERR "Done.\n";
	}else{
		unlink "$OUTDIR/sra2fastq_temp/${run_acc}.fastq.bz2";
	}
	print $total_size,"\n";

	if( $total_size < 50 ){
		`rm $OUTDIR/sra2fastq_temp/${run_acc}*gz`;
		print STDERR "failed to download FASTQ and convert FASTQ files from DDBJ.\n";
		return "failed";
	}

	return "success";
}

sub getEnaFastq {
	my ($info, $run_acc) = @_;

	print STDERR "Retrieving FASTQ for $run_acc from EBI-ENA...\n";
	
	foreach my $i ( keys %{$info->{ena_fastq_ftp}} ){
		my $url  = $info->{ena_fastq_ftp}->{$i}->{url};
		my $md5  = $info->{ena_fastq_ftp}->{$i}->{md5};
		my $size = $info->{ena_fastq_ftp}->{$i}->{size};
		
		my ($filename) = $url =~ /\/([^\/]+)$/;
						
		print STDERR "Downloading $url...\n";
		my $cmd = ($Download_tool =~ /wget/)? "$curl -O $OUTDIR/sra2fastq_temp/$filename \"ftp://$url\"":"$curl $ftp_proxy -o $OUTDIR/sra2fastq_temp/$filename \"ftp://$url\"";
		system("$cmd");
		
		#check downloaded file
		my $filesize = -s "$OUTDIR/sra2fastq_temp/$filename";

		if( !$filesize ){
			print STDERR "Failed to download $filename from $url.\n";
			return "failed";
		}
		if( $filesize != $size ){
			print STDERR  "$OUTDIR/sra2fastq_temp/$filename incompleted/corrupted -- file sizes mismatch.\n";
			return "failed";
		}

		#check md5
		my $md5sum = `md5sum $OUTDIR/sra2fastq_temp/$filename`;
		if(  $md5sum !~ /^$md5/ ){
			print STDERR "$OUTDIR/sra2fastq_temp/$filename corrupted -- md5 checksum mismatch.\n";
			return "failed";
		}
		
		print STDERR "Done.\n";
	}
	
	return "success";
}


sub getReadInfo {
	my ($acc, $readInfo) = @_;
	my $url;
	my $web_result;
	my @lines;

	print STDERR "Retrieving run(s) information from NCBI-SRA...\n";

	#get info from NCBI-SRA
	$url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$acc";
	print STDERR "Retrieving run acc# from NCBI-SRA $url...\n";
	my $cmd = ($Download_tool =~ /wget/)? "$curl -O - \"$url\" 2>/dev/null": "$curl $http_proxy \"$url\" 2>/dev/null";
	$web_result = `$cmd`;
	
	my @lines = split '\n', $web_result;
	print STDERR "$#lines run(s) found from NCBI-SRA.\n";

	foreach my $line (@lines){
		next if $line =~ /^Run/;
		chomp;

		my @f = split ',', $line;

		my $sub_acc  = $f[42]; #submission
		my $exp_acc  = $f[10]; #experiment_accession
		my $run_acc  = $f[0];  #run
		my $size_MB  = $f[7];  
		my $platform = $f[18]; #platform
		my $library  = $f[15]; #LibraryLayout
		my $url      = $f[9];  #download_path

		print STDERR "Run $run_acc has size 0 MB\n" if (!$size_MB);
	
		$readInfo->{$acc}->{$run_acc}->{exp_acc}  = $exp_acc;
		$readInfo->{$acc}->{$run_acc}->{sub_acc}  = $sub_acc;
		$readInfo->{$acc}->{$run_acc}->{platform} = $platform;
		$readInfo->{$acc}->{$run_acc}->{library}  = $library;
		$readInfo->{$acc}->{$run_acc}->{url}      = $url;

	}

	
	
	##### get info from EBI-ENA when NCBI-SRA fails ####
	
	print STDERR "Retrieving run(s) information from EBI-ENA...\n";
	
	$url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$acc&result=read_run";
	print STDERR "Retrieving run acc# from EBI-ENA $url...\n";
	$cmd = ($Download_tool =~ /wget/)? "$curl -O - \"$url\" 2>/dev/null": "$curl $http_proxy \"$url\" 2>/dev/null";
	$web_result = `$cmd`;
	die "ERROR: Failed to retrieve sequence information for $acc.\n" if $web_result !~ /^study_accession/;

	my @lines = split '\n', $web_result;
	print STDERR "$#lines run(s) found from EBI-ENA.\n";

	foreach my $line (@lines){
		next if $line =~ /^study_accession/;
		chomp;

		my @f = split '\t', $line;

		my $sub_acc  = $f[6]; #submission_accession
		my $exp_acc  = $f[4]; #experiment_accession
		my $run_acc  = $f[5]; #run_accession
		my $platform = $f[9]; #instrument_platform
		my $library  = $f[13]; #library_layout

		my @url  = split ';', $f[29]; #fastq_ftp
		my @md5  = split ';', $f[28]; #fastq_md5
		my @size = split ';', $f[27]; #fastq_bytes

		$readInfo->{$acc}->{$run_acc}->{platform} = $platform;
		$readInfo->{$acc}->{$run_acc}->{exp_acc}  = $exp_acc;
		$readInfo->{$acc}->{$run_acc}->{sub_acc}  = $sub_acc;
		$readInfo->{$acc}->{$run_acc}->{library}  = $library;
	
		for( my $i=0; $i<=$#url; $i++ ){
			$readInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{url}  = $url[$i];
			$readInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{md5}  = $md5[$i];
			$readInfo->{$acc}->{$run_acc}->{fastq_ftp}->{$i}->{size} = $size[$i];
		}
	}

	return $readInfo;
}
sub checkRequiredExec {
	die "ERROR: 'gzip' not found.\n" unless `which gzip`;
	if ($Download_tool=~ /curl/){
		die "ERROR: 'curl' not found.\n" unless `which curl`;
	}
	if ($Download_tool=~ /wget/){
		die "ERROR: 'wget' not found.\n" unless `which wget`;
	}
}

sub usage {
print <<__END__;

[DESCRIPTION]
    A script retrieves sequence project in FASTQ files from 
NCBI-SRA/EBI-ENA/DDBJ database using `curl` or `wget`. Input accession number
supports studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), 
samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions 
(SRA*/ERA*/DRA*).

[USAGE]
    $0 [OPTIONS] <Accession#> (<Accession# 2> <Accession# 3>...)

[OPTIONS]
    --outdir|d             Output directory
    --clean                clean up temp directory
    --platform-restrict    Only allow a specific platform
    --filesize-restrict    (in MB) Only allow to download less than a specific
	                       total size of files.
    --run-restrict         Only allow download less than a specific number
	                       of runs.
    --download-interface   curl or wget [default: curl]
    --help/h/?             display this help

__END__
exit();
}
