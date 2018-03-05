#!/usr/bin/perl
# contig_classifier.pl
# ver 0.1
# 2013/09/25
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# 2015 Feb
# - added accumulated coverage
# - added LCA

use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

$ENV{PATH} = "$Bin:$Bin/../../bin/:$Bin/../:$ENV{PATH}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'pacbio',
    'threads|t=s',
    'prefix|p=s',
    'db|d=s',
    'debug',
    'help|?') || &usage();

my ($INPUT, $PREFIX, $THREADS, $DB)
                      = ($opt{input}, $opt{prefix}, $opt{threads}, $opt{db});

if ( $opt{help} || !-e $INPUT ) { &usage(); }

my ($FILENAME) = $INPUT =~ /([^\/]+)$/;
my ($BASENAME,$EXT) = $FILENAME =~ /(.+)\.(\w+)$/;

#init options
$THREADS ||= 2;
$PREFIX  ||= "output";
$DB      ||= "/opt/apps/edge/database/bwa_index/NCBI-Bacteria-Virus.fna";

my $time = time;
my ($period,$mem);

my $abs_path = abs_path($INPUT);
executeCommand("mkdir -p temp$$");
executeCommand("ln -s $abs_path temp$$/$FILENAME");

$period = &timeInterval($time);
$mem    = &memUsage();
print "[$period] Running BWA\n";
if( $opt{pacbio} ){
    executeCommand("bwa mem -B5 -Q2 -E1 -a -M -t$THREADS $DB temp$$/$FILENAME > temp$$/$FILENAME.sam 2>>$PREFIX.log");
}
else{
    executeCommand("bwa mem -a -M -t $THREADS $DB temp$$/$FILENAME > temp$$/$FILENAME.sam 2>>$PREFIX.log");
}

$period = &timeInterval($time);
print "[$period] Filtering unmapped contigs\n";
executeCommand("samtools view -F4  temp$$/$FILENAME.sam > temp$$/$FILENAME.mapped.sam 2>>$PREFIX.log");

$period = &timeInterval($time);
print "[$period] Splitting SAM file\n";
executeCommand("split_sam_by_lines.pl --line 20000 --input temp$$/$FILENAME.mapped.sam");

$period = &timeInterval($time);
print "[$period] Classifying contigs\n";
executeCommand("cd temp$$; parallel --results ${PREFIX}_para_log -j $THREADS 'seq_coverage.pl --input {}' ::: *.part* > ../$PREFIX.ctg_class.csv 2>>../$PREFIX.log");
#executeCommand("parallel --results temp$$/${PREFIX}_para_log -j $THREADS 'seq_coverage.pl --input {}' ::: temp$$/*.part* > $PREFIX.ctg_class.csv 2>>$PREFIX.log");

$period = &timeInterval($time);
print "[$period] Reporting unclassified contigs\n";
executeCommand("samtools view -f4 -S temp$$/$FILENAME.sam 2>>$PREFIX.log | awk -F\\\\t '{print \$1\"\\tunclassified\\tunclassified\\t\\t\\t\",length(\$10),\"\\t\\t\\t\\t\\t0\\t\\t0\\t\\t\\t\",length(\$10)}' >> $PREFIX.ctg_class.csv");
executeCommand("samtools view -f4 -S temp$$/$FILENAME.sam 2>>$PREFIX.log | awk '{print \">\"\$1\"\\n\"\$10}' > $PREFIX.unclassified.fasta");

$period = &timeInterval($time);
print "[$period] Merging classification\n";
executeCommand("merge_coverage.pl $PREFIX.ctg_class.csv $PREFIX > $PREFIX.assembly_class.csv 2>>$PREFIX.log");

$period = &timeInterval($time);
print "[$period] Reporting classification (BEST hit)\n";
executeCommand("class_top_hit_summary.pl < $PREFIX.assembly_class.csv > $PREFIX.assembly_class.top.csv 2>>$PREFIX.log");
executeCommand("class_top_hit_summary.pl < $PREFIX.ctg_class.csv > $PREFIX.ctg_class.top.csv 2>>$PREFIX.log");

$period = &timeInterval($time);
print "[$period] Reporting classification (LCA)\n";
executeCommand("report_LCA.pl < $PREFIX.ctg_class.csv > $PREFIX.ctg_class.LCA.csv 2>>$PREFIX.log");
executeCommand("(head -n 1 $PREFIX.ctg_class.LCA.csv && tail -n +2 $PREFIX.ctg_class.LCA.csv | sort -t '\t' -k6nr) > $PREFIX.ctg_class.LCA.csv.sort");
executeCommand("mv $PREFIX.ctg_class.LCA.csv.sort $PREFIX.ctg_class.LCA.csv");

my ($tol_contigs_count, $tol_contigs_bases, $classified_contigs_count, $classified_contigs_bases, $unclassified_contigs_count, $unclassified_contigs_bases)=&countResult("$PREFIX.ctg_class.top.csv");
$period = &timeInterval($time);
printf ("[$period] Total Contigs: %d (%d bp); Classified Contigs: %d (%d bp); Unclassified Contigs: %d (%d bp);\n",$tol_contigs_count,$tol_contigs_bases,$classified_contigs_count,$classified_contigs_bases,$unclassified_contigs_count,$unclassified_contigs_bases);

unless( $opt{debug} ){
    $period = &timeInterval($time);
    print "[$period] Cleaning temporary files\n";
    executeCommand("rm -rf temp$$/");
}


$period = &timeInterval($time);
print "[$period] Finished. Please find the result for contig in $PREFIX.ctg_class.csv and $PREFIX.assembly_class.csv\n";

sub countResult{
	my $file=shift;
	my ($unclassified_contigs_count,$unclassified_contigs_bases,$classified_contigs_count,$classified_contigs_bases) = (0,0,0,0);
	my $tol_contigs_bases=0;
	open(my $fh, $file) or die "$!\n";
	while(<$fh>){
		chomp;
		next if /^#/;
		my @fields = split /\t/,$_;
		if ($fields[1] eq "superkingdom"){
			$classified_contigs_bases+=$fields[-1];
			$classified_contigs_count++;
			$tol_contigs_bases += $fields[5];
		}
		if ($fields[1] eq "unclassified" && $fields[5] == $fields[-1]){
			$unclassified_contigs_count++;
			$unclassified_contigs_bases+=$fields[-1];
			$tol_contigs_bases += $fields[5];
		}
		
	}
	close $fh;
	my $tol_contigs_count = $classified_contigs_count + $unclassified_contigs_count;
	return ($tol_contigs_count, $tol_contigs_bases, $classified_contigs_count, $classified_contigs_bases, $unclassified_contigs_count, $unclassified_contigs_bases);
}


sub timeInterval{
    my $now = shift;
    $now = time - $now;
    return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}

sub memUsage {
    open STAT , "</proc/$$/statm" or die "Unable to open stat file";
    my @stat = split /\s+/ , <STAT>;
    close STAT;

    return sprintf "%.2f", ($stat[1]/1024/1024);
}

sub executeCommand 
{
    my $command = shift;
	print "[DEBUG] $command\n" if $opt{debug};
    system($command) == 0
         || die "COMMAND FAILED: $command\n";
}


sub usage {
print <<__END__;

$0 [OPTIONS] --input <FILENAME>

    --input   | -i <STRING>  the sequence of contigs in FASTA format

[OPTIONS]

    --threads | -t <NUM>     the number of threads (default: 2)
    --prefix  | -p <STRING>  the prefix of the output file (default: output)
    --db      | -d <STRING>  the path of BWA index
                             (default: /opt/apps/edge/database/bwa_index/NCBI-Bacteria-Virus.fna)
    --pacbio                 use this option to treat pacbio reads as contigs 
    --debug                  keep all temp files
    --help/h/?               display this help                   

__END__
exit();
}
