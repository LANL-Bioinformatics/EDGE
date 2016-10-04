#!/usr/bin/env perl

use strict;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../../lib";
use SNPphy;

$ENV{PATH} = "$RealBin:$RealBin/../../bin/:$ENV{PATH}";
$ENV{PERL5LIB} = "$RealBin/../../lib:$ENV{PERL5LIB}";

$|=1;
my $refdir;
my $workdir;
my $outdir;
my $outfile;
my $rsignal=1;
my $reference;
my $name;
my $suffix;
my $path;
my $gsignal=0;
my $genbank;
my $time=2;
my $data;
my $tree;
my $modeltest;
my $pselection;
my $genefile;
my $clean;
my $threads;
my $maxthreads;
my $cutoff;
my $bindir;
my %fasta_list;
my %contig_list;
my %read_list;
my $incompleteRun=0; # can be deleted after correction of "run coninuation".
my $nucmer=0;
my $contig_nucmer=0;
my $read_mapping=0;
my $buildSNP=0;
my $buildtree=0;
my $evolanal=0;
my $mappingGaps;
my $ptree;
my $ps=0;

#my $control= "SNPphy.ctl";
my $control = $ARGV[0];
$bindir=getBinDirectory();

open(CTL, "$control")||die "Please provide a control file in the working directory";
while (<CTL>){
   if (/refdir\s*=\s*(\S+)\s*#{0,1}.*$/){
      $refdir=$1;
      if ($refdir=~/.+\/$/){my $temp= chop($refdir);}
   }
   if (/workdir\s*=\s*(\S+)\s*#{0,1}.*$/){
      $workdir=$1;
      if ($workdir=~/.+\/$/){my $temp= chop($workdir);}
      $outdir=$workdir.'/results';
      if (! -e $outdir){mkdir "$outdir";}
   }
   if (/outfile\s*=\s*(\S+)\s*#{0,1}.*$/){$outfile=$1;}
   if (/reference\s*=\s*(0|1)\s*#{0,1}.*$/){$rsignal=$1;}
   if ($rsignal==1 && /reffile\s*=\s*(\S+)\s*#{0,1}.*$/){
      $reference="$refdir/$1";
      ($name,$path,$suffix) = fileparse("$reference", qr/\.[^.]*/);
   }
   elsif($rsignal==0){  # random select a ref file
      opendir(DIR, $refdir);
      my @reffiles = grep { /.[fna|fa|fasta|fsa]$/ && -f "$refdir/$_" } readdir(DIR);
      $reference= "$refdir/". $reffiles[int(rand(scalar(@reffiles)))];
      ($name,$path,$suffix) = fileparse("$reference", qr/\.[^.]*/);
      close DIR;
   }
   if (/cdsSNPS\s*=\s*(0|1)\s*#{0,1}.*$/){$gsignal=$1;}
#     0=No; 1=Yes    
   if ($gsignal==1){
      if (/gbfile\s*=\s*(\S+)\s*#{0,1}.*$/){
         $genbank=$1;
         $genbank=$refdir.'/'.$genbank;
      }
   }
   if (/FirstTime\s*=\s*(1|2)\s*#{0,1}.*$/){
      $time=$1;
#      print "$time\n";
      if($time==2){if (/alignfile\s*=\s*(\S+)\s*#{0,1}.*$/){$outfile=$1;}}
   }
   if (/data\s*=\s*(\d)\s*#{0,1}.*$/){$data=$1;}
#     0=Finished; 1=contig; 2=reads; 3=F+C; 4=F+R; 5=C+R 6=F+C+R
   if (/tree\s*=\s*(0|1|2|3)\s*#{0,1}.*$/){$tree=$1;}
#     0=no tree; 1=fasttree; 2=raxml; 3=both
   if (/modelTest\s*=\s*(0|1)\s*#{0,1}.*$/){$modeltest=$1;}
   if (/PosSelect\s*=\s*(0|1|2|3)\s*#{0,1}.*$/){$pselection=$1;} 
   if (/genefile\s*=\s*(\S+)\s*#{0,1}.*$/){
      if ($1!~/^$refdir/){$genefile=$refdir.'/'.$1;}
      else {$genefile=$1;}
   }
   if (/clean\s*=\s*(0|1)\s*#{0,1}.*$/){$clean=$1;}
#     0=no clean; 1=clean
   if (/threads\s*=\s*(\d+)\s*#{0,1}.*$/){
      $threads=$1;
      my $maxthreads = ($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
      if ($threads <1 || $threads>$maxthreads){die ("-thread value must be between 1 and $maxthreads.\n");}
   }
   if (/cutoff\s*=\s*(\d+)\s*#.*$/){$cutoff=$1;}
}
close CTL;
my $error=$outdir.'/error.log';
print "Checking directories and files... ";
my $check=SNPphy::check($workdir,$refdir,$time,$data);
if ($check==0){
   if ($data==0){$nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==1){$contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==2){$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==3){$nucmer=1;$contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==4){$nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==5){$contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==6){$nucmer=1;$contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==7){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
}
if($time==1){
   if ($check==1){
      print "Continuing run...\n";
      if ($data<=2){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==3){$contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==4||$data==5){$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){$contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==2){
      print "Continuing run...\n";
      if ($data<=2){$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==3||$data==4||$data==5){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==3){
      if ($data<=2){print "Run is complete... Nothing to do\n";}
      if ($data==3||$data==4||$data==5){print "Continuing run...\n";$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){print "Continuing run...\n";$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==4){
      if ($data<=5){print "Run is complete... Nothing to do\n";}
      if ($data==6){print "Continuing run...\n";$buildtree=1;if ($pselection>0){$ps=1;}}
   }
}
if($time==2 && $check==1){exit;}
if($time==1 && $check==0){print "Complete\n";}

if ($check==0){print "Preparing Files... ";}
#if ($time==1 && $check==0 && !-e $workdir.'/files'){`mkdir $workdir/files`;}
if ($check==0 && !-e $workdir.'/files'){`mkdir $workdir/files`;}

opendir(DIR,$refdir);
while (my $files=readdir(DIR)){
   next if ($files=~/^..?$/);
   my ($rname,$rpath,$rsuffix)=fileparse("$files",qr/\.[^.]*/);
   if ($nucmer==1){
      if ($files=~/.+\.fn|s?a?s?t?a$/){
         #next if ($files=~ /contigs?/||$files=~ /fa?s?t?q/);
         next if ($files=~ /fa?s?t?q/);
         my $fasta=$refdir.'/'.$files;
         $fasta_list{SNPphy::prepareComplete($refdir,$workdir,$fasta,$rname)}++;
      }
   }
}
close DIR;

opendir(DIR,$workdir);
while (my $files=readdir(DIR)){
   next if ($files=~/^..?$/);
   my ($qname,$qpath,$qsuffix)=fileparse("$files",qr/\.[^.]*/);
   if ($contig_nucmer==1){
      if ($files=~ /.+\.contigs?$/){
         next if ($files=~ /\.fn|s?a?s?t?a$/ || $files=~ /\.fa?s?t?q$/);
         my $contig=$workdir.'/'.$files;
         $contig_list{SNPphy::prepareContig($workdir,$contig,$qname)}++;
      }
   }
   if ($read_mapping==1){
      if ($files=~ /.+\.(fa?s?t?q)$/){
         my $suffix = $1;
         next if ($files=~ /.+\.fn|s?a?s?t?a$/ || $files=~ /\.contigs?/);
         #my $fastq=$refdir.'/'.$files;
         if ($qname=~/(.+)_R\d$/ || $qname=~/(.+)_SE$/){
            my $read=$1."_read.$suffix";
            $read_list{$read}++;
            print $read,"\n";
         }
      }
   }
}
print "Complete\n";

my $snpdir=$outdir.'/snps';
my $gapdir=$outdir.'/gaps';
my $statdir=$outdir.'/stats';

if ($time==1){open(ALL,">$workdir/working_list.txt")||die "$!";}
#if ($time==1 && $check==1){open(ALL,">$workdir/working_list.txt")||die "$!";}
elsif($time==2){open (ALL, ">>$workdir/working_list.txt")||die "$!";}

if ($nucmer==1){
   open (FAS, ">$workdir/fasta_list.txt")||die "$!";
   foreach my $names(keys %fasta_list){print FAS "$names\n";}
   foreach my $names(keys %fasta_list){print ALL "$names\n";}
#  foreach my $names(keys %fasta_list){print "$names\n";}
   print "Running NUCmer on complete genomes\n";
   ($time==2)?
   SNPphy::completeNUCmer($workdir,$bindir,$reference,"$workdir/working_list.txt",$threads,$error):
   SNPphy::completeNUCmer($workdir,$bindir,$reference,"$workdir/fasta_list.txt",$threads,$error);
}
if ($contig_nucmer==1){
   open (CON, ">$workdir/contigs_list.txt")||die "$!";
   foreach my $names(keys %contig_list){print CON "$names\n";}
   foreach my $names(keys %contig_list){print ALL "$names\n";}
   print "Running NUCmer on contigs\n";
   SNPphy::contigNUCmer($workdir,$bindir,"$workdir/contigs_list.txt",$threads,$reference,$error);
}
if ($read_mapping==1){
   open (READ, ">$workdir/reads_list.txt")||die "$!";
   foreach my $names(keys %read_list){print ALL "$names\n";}
   foreach my $names(keys %read_list){print READ "$names\n";}
   print "Preparing to map reads\n";
   
   if ($time==2){
      print "Identifying all gaps\n";
      $mappingGaps=SNPphy::identifyGaps($outdir,"$workdir/fasta_list.txt",$name,"map");
      SNPphy::removeGaps($bindir,$reference,$mappingGaps);
      SNPphy::readsMapping($workdir,$bindir,"$workdir/reads_list.txt",$threads,$name,$error);
   }
   elsif ($time==1){
      my $tempdir=$outdir.'/temp/';
      `cp $reference $tempdir`;
      SNPphy::readsMapping($workdir,$bindir,"$workdir/reads_list.txt",$threads,$name,$error);
   }
}

if ($gsignal==1 && $buildSNP==1){
   print "GFF file provided, SNPs will be differentiated as coding vs noncoding\n";
   my ($genname,$genpath,$gensuffix)=fileparse("$genbank",qr/\.[^.]*/);
   SNPphy::genbank($outdir,$genbank,$genname);
}
if ($buildSNP==1){
   SNPphy::identifyGaps($outdir,"$workdir/working_list.txt",$name,"snp");
   SNPphy::buildSNPDB($outdir,$bindir,$reference,"$workdir/working_list.txt",$outfile,$gsignal,$error);
}
if ($buildtree==1){
#   if ($tree==2||$tree==3){SNPphy::convertFasta($outdir,$outfile);}
#   if ($gsignal==1&&($tree==2||$tree==3)){SNPphy::convertFasta($outdir,$outfile);}
   if (($tree==2||$tree==3)&&($modeltest==1)){SNPphy::modeltest($refdir,$outfile,$threads,$error);}
   SNPphy::buildTree($outdir,$threads,$outfile,$tree,$name,$error);
}

my $tinfo=$refdir.'/RAxML_info.'.$name;
my $tlog=$refdir.'/RAxML_log.'.$name;
my $tparsimony=$refdir.'/RAxML_parsimonyTree.'.$name;
my $tresult=$refdir.'/RAxML_result.'.$name;
my $tbest=$refdir.'/RAxML_bestTree.'.$name;
if (-e $tinfo){`mv $tinfo $outdir/`;}
if (-e $tlog){`mv $tlog $outdir/`;}
if (-e $tparsimony){`mv $tparsimony $outdir/`;}
if (-e $tresult){`mv $tresult $outdir/`;}
if (-e $tbest){`mv $tbest $outdir/`;}
if ($tree==2||$tree==3){$tbest=$outdir.'/RAxML_bestTree.'.$name;}
if ($tree==1){$tbest=$outdir.'/snp_alignment_all_'.$name.'.fasttree';}

if ($clean){SNPphy::clean($outdir);}
if ($ps==1){
   my $stats_file=$outdir.'/snp_stats_all_'.$name.'.txt';
   my $genedir= $outdir.'/PSgenes';
   my $pamldir= $outdir.'/paml';
   if (!-d $pamldir){`mkdir $pamldir`;}
   `cp $tbest $pamldir`;
   if ($tree==2||$tree==3){$ptree= $pamldir.'/RAxML_bestTree.'.$name;}
   if ($tree==1){$ptree=$pamldir.'/snp_alignment_all_'.$name.'.fasttree';}
   elsif ($tree==0){print "Need to build a tree before any evolutionary analyses can be performed.\n"; exit;}
   my $gapfile=$outdir.'/all_gaps.txt';

   SNPphy::extractGenes($outdir,$stats_file,$reference,$bindir,"$workdir/working_list.txt",$threads,$gapfile,$genefile,$error);
   SNPphy::translateGenes($outdir,$bindir,$threads,"translate",$error);
   SNPphy::alignGenes($outdir,$bindir,$threads,"mafft",$error);
   SNPphy::revTransGenes($outdir,$bindir,$threads,"pal2nal",$error);
   SNPphy::positiveSelection($outdir,$bindir,$ptree,0,"Sites","1,2",$threads,$error);
   SNPphy::positiveSelection($outdir,$bindir,$ptree,2,"BrSites",2,$threads,$error);
}

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

