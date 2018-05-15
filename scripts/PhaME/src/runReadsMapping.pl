#!/usr/bin/perl -w

######################################################
# Written by Sanaa Ahmed
# Nov. 30, 2012

# Given a directory containing paired read files, runs bowtie
######################################################

use strict;
use FindBin qw($RealBin);
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

# set up environments
$ENV{PATH}="$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

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
my $max_thread=($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
if ($thread<1 || $thread>$max_thread){die("-thread value must be between than 1 and $max_thread.\n");}

if (! -e $outdir){mkdir "$outdir";}
if ($outdir=~ /.+\/$/){my $temp= chop($outdir);}
chdir $outdir;

if ($indir=~ /.+\/$/){my $temp= chop($indir);}
my $querydir= $indir;
 
read_directory($querydir);

my $pm=new Parallel::ForkManager($thread);

# data structure retrieval and handling from multi-threading 
$pm->run_on_finish ( # called BEFORE the first call to start()
   sub{
      my ($pid,$ident)=@_;
 #     print "Process (Ident: $ident) (pid: $pid) core dumped.\n";
#      print "Exit code\t$exit_code\nExit signal\t$exit_signal\nCore dump\t$core_dump\n";
      opendir (TEMP, "$outdir");
      my $lines =0;
      while (my $file= readdir(TEMP)){
         if ($file=~ /.+\.gaps$/){
            my $gap_file= $outdir.'/'.$file;
            open(FILE, "$gap_file");
            while (<FILE>){$lines++;}
            if ($lines==1){`rm -f $gap_file`; print "Removed $gap_file\n"; $lines=0;}
            else {`mv $gap_file $outdir/gaps`;}
         }
#         if ($file=~/.+\.vcf$/){
#            print "Moving vcf files to snps/ directory\n";
#            my $vcf_file=$outdir.'/'.$file;
#            `mv $vcf_file $outdir/snps`;
#         }
      }
   }
);

foreach my $command (@command){
   $pm->start and next;
   if (system($command)){die "Error running $command.\n";}
   $pm->finish;
}
$pm->wait_all_children;

print "Read Mapping complete\n";

sub read_directory
{
my $dir=shift;
my ($name,$path,$suffix,$file,$query);
my $ref=0;
my @fastq;
my %queries;
if ($reference){($name,$path,$suffix)=fileparse("$reference",qr/\.[^.]*/);}

open (IN,$list);
while (<IN>){
   chomp;
   $queries{$_}++;
}
close IN;

my $temp=0;
my %check;
opendir(PARENT,$dir);
while (my $files= readdir(PARENT)){  
   next if ($files=~ /^..?$/);
   if ($files!~ /$name/ && ($files=~ /(.+)[_.]R1.*\.fa?s?t?q$/ ||  $files=~ /(.+)[_.]1.*\.fa?s?t?q$/)){
      $temp= $1.'_pread';
#      print "$temp\n";
      if (exists $queries{$temp} && !exists $check{$temp}){
         $check{$temp}++;
         $query=$dir.'/'.$files;
         $prefix= "$name";
         create_bowtie_commands($query,$prefix,$temp);
      }
      $temp= $1.'_read';
      if (exists $queries{$temp} && !exists $check{$temp}){
         $check{$temp}++;
         $query=$dir.'/'.$files;
         $prefix= "$name";
         create_bowtie_commands($query,$prefix,$temp);
      }
      $temp= $1.'_sread';
      if (exists $queries{$temp} && !exists $check{$temp}){
         $check{$temp}++;
         $query=$dir.'/'.$files;
         $prefix= "$name";
         create_bowtie_commands($query,$prefix,$temp);
      }

   }
 #  if ($files!~ /$name/ && $files=~ /(.+)[_.]1.*\.fa?s?t?q$/){
  #    $temp=$1.'_pread';
  #    if (exists $queries{$temp} && !exists $check{$temp}){
   ##      $check{$temp}++;
    #     $query=$dir.'/'.$files;
     #    $prefix=$name;
#         print "$temp\n";
      #   create_bowtie_commands($query,$prefix,$temp);
     # }
  # }
   if ($files!~ /$name/ && $files=~ /(.+)\.fa?s?t?q$/){
      $temp= $1.'_sread';
      if (exists $queries{$temp} && !exists $check{$temp}){
         $check{$temp}++;
         $query=$dir.'/'.$files;
         $prefix= "$name";
         create_bowtie_commands($query,$prefix,$temp);
      }
   } 
}
closedir(PARENT);
}

sub create_bowtie_commands
{
my $read=shift;
my $prefix=shift;
my $temp=shift;
my $read1;
my $read2;
my $readu;

my ($name,$path,$suffix)=fileparse("$read",qr/\.[^.]*/);
if ($temp=~/pread/i){
#   print "$read\n";
   if ($name=~ /(.+)([_.]R)(\d)(.*)$/){
      $prefix.= "\_$1";
      $read1=$path.$1.$2.$3.$4.$suffix;
      $read2=$path.$1.$2.'2'.$4.$suffix;
   }
   if ($name=~ /(.+)([_.])(1)(.*)$/){
      $prefix.= "\_$1";
      $read1=$path.$1.$2.$3.$4.$suffix;
      $read2=$path.$1.$2.'2'.$4.$suffix;
   }

   my $bowtie_command= "runReadsToGenome.pl -p $read1,$read2 -ref $reference -pre $prefix -d $outdir -aligner $aligner" ;
#   print "READ1:  $read1\nREAD2:  $read2\n$bowtie_command\n\n\n";
   push (@command,$bowtie_command);
}
elsif ($temp=~/sread/i){
   $prefix.="\_$name";
   my $bowtie_command= "runReadsToGenome.pl -u $read -ref $reference -pre $prefix -d $outdir -aligner $aligner" ;
    $bowtie_command= "runReadsToGenome.pl -long $read -ref $reference -pre $prefix -d $outdir -aligner $aligner" if ($aligner =~ /minimap/);
#   print "READ1:  $read1\n$bowtie_command\n\n\n";
   push (@command,$bowtie_command);
}elsif ($temp=~/_read/i){
#   print "$read\n";
   if ($name=~ /(.+)([_.]R)(\d)(.*)$/){
      $prefix.= "\_$1";
      $read1=$path.$1.$2.$3.$4.$suffix;
      $read2=$path.$1.$2.'2'.$4.$suffix;
      $readu=$path.$1.$suffix;
   }
   elsif ($name=~ /(.+)([_.])(1)(.*)$/){
      $prefix.= "\_$1";
      $read1=$path.$1.$2.$3.$4.$suffix;
      $read2=$path.$1.$2.'2'.$4.$suffix;
      $readu=$path.$1.$suffix;
   }

   my $bowtie_command= "runReadsToGenome.pl -p $read1,$read2 -u $readu -ref $reference -pre $prefix -d $outdir -aligner $aligner" ;
#   print "READ1:  $read1\nREAD2:  $read2\n$bowtie_command\n\n\n";
   push (@command,$bowtie_command);
}
}

sub usage
{
print <<Usage;
$0 -r <reference_fasta> -q <query_dir> -d <output_directory> -t <# of threads> -a <aligner bwa|bowtie>

Usage
exit;
}
