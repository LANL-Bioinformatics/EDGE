#!/usr/bin/perl -w

######################################################
# Written by Sanaa Ahmed
# Nov. 30, 2012

# Given a directory containing paired read files, runs bowtie
######################################################

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../lib";
use Parallel::ForkManager;

my ($indir,$reference,$prefix,$thread,$list,$aligner);
my @command;
my $outdir=`pwd`;
   $outdir =~ s/\n//;

GetOptions(
   'q|querydir=s'    => \$indir,
   'r|reference=s'   => \$reference,
   'd|outdir=s'      => \$outdir,
   't|thread=i'      => \$thread,
   'l|list=s'        => \$list,
   'a|aligner=s'     => \$aligner,
   'h|help'          => sub{usage()}
);

&usage() unless ($indir && $aligner);
my $max_thread = ($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
if ($thread<1 || $thread>$max_thread){die("-thread value must be between than 1 and $max_thread.\n");}

if (! -e $outdir){mkdir "$outdir";}
if ($outdir=~ /.+\/$/){my $temp= chop($outdir);}
chdir $outdir;

if ($indir=~ /.+\/$/){my $temp= chop($indir);}
my $querydir= $indir;
my $bindir=getBinDirectory();
 
read_directory($querydir);

my $pm=new Parallel::ForkManager(1);

# data structure retrieval and handling from multi-threading 
$pm->run_on_finish ( # called BEFORE the first call to start()
   sub{
      my ($pid,$ident)=@_;
   #   print "Process (Ident: $ident) (pid: $pid) core dumped.\n";
#      print "Exit code\t$exit_code\nExit signal\t$exit_signal\nCore dump\t$core_dump\n";
      opendir (TEMP, "$outdir");
      my $lines =0;
      while (my $gaps= readdir(TEMP)){
         if ($gaps=~ /.+\.gaps$/){
            my $gap_file= $outdir.'/'.$gaps;
            open(FILE, "$gap_file");
            while (<FILE>){$lines++;}
            if ($lines==1){`rm $gap_file`; print "Removed $gap_file\n"; $lines=0;}
            else {`mv $gap_file $outdir/gaps`;}
         } 
         if ($gaps=~ /.+\.vcf$/){
            my $gap_file= $outdir.'/'.$gaps;
            `mv $gap_file $outdir/snps`;
         }
      }
   }
);

foreach my $command (@command){
   $pm->start and next;
   print $command,"\n";
   if (system($command)){die "Error running $command.\n";}
   $pm->finish;
}
$pm->wait_all_children;
print "Read Mapping complete\n";

sub read_directory
{
my $dir=shift;
my ($ref_name,$path,$suffix,$file,$query);
my $ref=0;
#my @fastq;
#my %queries;
if ($reference){($ref_name,$path,$suffix)=fileparse("$reference",qr/\.[^.]*/);}

open (IN,$list);
while (<IN>){
   chomp;
 #  $queries{$_}++;
   my ($query_prefix,$query_suffix) = $_ =~ /(.+)_read\.(\w+)$/;
   create_bowtie_commands($query_prefix,$query_suffix,$ref_name);
}
close IN;

#opendir(PARENT,$dir);
#while (my $files= readdir(PARENT)){  
#   next if ($files=~ /^..?$/);
#   if ($files=~ /(.+)_R.\.fastq$/ || $files=~ /(.+)_SE\.fastq$/ )
#   {
#      my $temp= $1.'_read';
#      if (exists $queries{$temp}){
#         $query=$dir.'/'.$files;
#      print "$query\n";
#         my ($qname,$qpath,$qsuffix)=fileparse("$query",qr/\.[^.]*/);
#         $prefix= "$ref_name";
#      print "$qname\t$prefix\n";
#         create_bowtie_commands($query,$prefix);
#      }
 #  } 
#}
#closedir(PARENT);
}

sub create_bowtie_commands
{
    my $query_prefix=shift;
    my $query_suffix=shift;
    my $ref_prefix=shift;
    my $prefix = "${ref_prefix}_$query_prefix";
    my ($read1, $read2, $se_read);
    if ( -e "$indir/${query_prefix}_R1.$query_suffix")
    {
        $read1 = "$indir/${query_prefix}_R1.$query_suffix" ;
        $read2 = "$indir/${query_prefix}_R2.$query_suffix" ;
    }
    $se_read = "$indir/${query_prefix}_SE.$query_suffix" if  (-e "$indir/${query_prefix}_SE.$query_suffix" );
    my $bowtie_command= "$bindir/runReadsToGenome.pl -ref $reference -pre $prefix -d $outdir -no_plot -aligner $aligner -bowtie_options \"-p $thread -a  \" " ;
    $bowtie_command .= " -p $read1,$read2 " if ($read1);
    $bowtie_command .= " -u $se_read " if ($se_read);
#   print "READ1:  $read1\nREAD2:  $read2\n$bowtie_command\n";
    push (@command,$bowtie_command);

}

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

sub usage
{
print <<Usage;
$0 -r <reference_fasta> -q <query_dir> -d <output_directory> -t <# of threads> -a <aligner bwa|bowtie>

Usage
exit;
}
