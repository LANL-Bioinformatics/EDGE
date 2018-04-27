#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin $RealBin);
use lib "$Bin";
use lib "$RealBin/../lib/";
use lib "$RealBin/../ext/lib/perl5";
use File::Basename;
use PhaME;

$|=1;

# set up environments
$ENV{PATH}="$RealBin:$RealBin/../ext/bin:$ENV{PATH}";
$ENV{PERL5LIB} = ($ENV{PERL5LIB})?"$ENV{PERL5LIB}:$RealBin/../lib:$RealBin/../ext/lib/perl5":"$RealBin/../lib:$RealBin/../ext/lib/perl5";

=head

DESCRIPTION
   runPhaME.pl is a wrapper to run the PhaME analysis pipeline 

USAGE
   Control file "phame.ctl" needs to be copied and editted in the working director. 
   This wrapper runs the PhaME pipeline based on the settings of the control file. 

REQUIREMENTS
   See PhaME manual for a complete list of system and software requirements

KNOWN BUGS
   Working script 

COPYRIGHT/LICENCE

AVAILABILITY
   
AUTHOR
   Sanaa Ahmed
   2015/01/15
   B-11
   Los Alamos National Laboratory
   sahmed@lanl.gov
=cut

my $refdir;
my $workdir;
my $outdir;
my $project;
my $rsignal=0;
my $reference;
my $aligner="bowtie";
my $name;
my $path;
my $suffix;
my $gsignal=0;
my $time=1;
my $data=0;
my $reads=2;
my $tree=2;
my $modeltest=0;
my $bsignal=0;
my $bootstrap=100;
my $pselection=0;
my $code=0;
my $type;
my $clean=1;
my $threads=1;
my $cutoff=0;
my $annotation;
my $genefile;
my $runtime=time;
my $nucmer=0;
my $contig_nucmer=0;
my $read_mapping=0;
my $buildSNP=0;
my $buildtree=0;
my $ps=0;
my $bs=0;
my %fasta_list;
my %contig_list;
my %read_list;
my $mappingGaps;
my $ptree;

my $control= $ARGV[0] || "phame.ctl";
my $bindir=getBinDirectory();

## Read in control file 
open(CTL, "$control")||die "Please provide a control file";
while (<CTL>){
   if (/refdir\s*=\s*(\S+)\s*#{0,1}.*$/){
      $refdir=$1;
      if ($refdir=~/.+\/$/){my $temp= chop($refdir);}
   }
   if (/workdir\s*=\s*(\S+)\s*#{0,1}.*$/){
      $workdir=$1;
      if ($workdir=~/.+\/$/){my $temp= chop($workdir);}
      $outdir=$workdir.'/results';
      if (! -e $outdir){`mkdir -p $outdir`;}
   }
   if (/project\s*=\s*(\S+)\s*#{0,1}.*$/){$project=$1;}
   if (/aligner\s*=\s*(\S+)\s*#{0,1}.*$/){$aligner=$1;}
   if (/reference\s*=\s*(0|1)\s*#{0,1}.*$/){$rsignal=$1;}
   if ($rsignal==1 && /reffile\s*=\s*(\S+)\s*#{0,1}.*$/){$reference="$refdir/$1";}
   elsif($rsignal==0){
      opendir(DIR, $refdir);
      my @reffiles = grep { /.[fna|fa|fasta|fsa]$/ && -f "$refdir/$_" } readdir(DIR);
      $reference= "$refdir/". $reffiles[int(rand(scalar(@reffiles)))];
      closedir DIR;
   }
   if (/cdsSNPS\s*=\s*(0|1)\s*#{0,1}.*$/){$gsignal=$1;}
#     0=No; 1=Yes    
   if (/FirstTime\s*=\s*(1|2)\s*#{0,1}.*$/){$time=$1;}
   if (/data\s*=\s*(\d)\s*#{0,1}.*$/){$data=$1;}
#     0=Finished; 1=contig; 2=reads; 3=F+C; 4=F+R; 5=C+R 6=F+C+R
   if (/\s+reads\s*=\s*(1|2|3)\s*#{0,1}.*$/){$reads=$1;}
   if (/tree\s*=\s*(0|1|2|3)\s*#{0,1}.*$/){$tree=$1;}
#     0=no tree; 1=fasttree; 2=raxml; 3=both
   if (/modelTest\s*=\s*(0|1)\s*#{0,1}.*$/){$modeltest=$1;}
   if (/bootstrap\s*=\s*(0|1)\s*#{0,1}.*$/){$bsignal=$1;}
   if (/N\s*=\s*(\d+)\s*#{0,1}.*$/){$bootstrap=$1;}
   if (/PosSelect\s*=\s*(0|1|2|3)\s*#{0,1}.*$/){$pselection=$1;}

   if (/code\s*=\s*(0|1|2)\s*#{0,1}.*$/){$code=$1;}
   if ($code==1){$type="virus";}
   elsif ($code==2){$type="eukaryote";}
   else{$type="bacteria";}
#     0=Bacteria; 1=Virus; 2=Eukaryote  

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

($name,$path,$suffix)=fileparse("$reference",qr/\.[^.]*/);
if ($gsignal==1){$annotation="$refdir/$name.gff";}
if ($pselection>0){$genefile="$refdir/$name.gff";}

my $error="$outdir/$project.error";
my $logfile="$outdir/$project\_PhaME.log";

&print_timeInterval($runtime,"\tLoading information\n");
print "\tRefd:\t$refdir\n";
print "\tWorkd:\t$workdir\n";
print "\tOutd:\t$outdir\n";
if (!-d $refdir || !-d $workdir || !-d $outdir){
   print "One or more of the above directories does not exist.\nPlease provide correct directory paths\n";
}

print "\tReference:\t$reference\n";
#if (!-e "$reference"){print "File $reference does not exist.\nPlease provide correct reference\n";}
if ($gsignal==1){
   print "\tAnnotation:\t$annotation\n";
   if ( !-e $annotation){print "File $annotation does not exist.\nPlease provide correct annotation file\n";}
}
if ($pselection>0){
   print "\tGenes:\t$genefile\n";
   if (-e $reference && !-e $genefile){print "File $genefile does not exist.\nPlease provide correct functional protein file\n";}
}

print "\tCode:\t$type\n";
print "\tLog file:\t$logfile\n";
print "\tError file:\t$error\n";

&print_timeInterval($runtime,"\tChecking directories and files... ");
my $check=PhaME::check($workdir,$refdir,$time,$data,$name,$logfile,$project);
if ($data==7){$check=0;}
#print "\nCHECK:\t$check\n";
my $snpdir=$outdir.'/snps';
my $gapdir=$outdir.'/gaps';
my $statdir=$outdir.'/stats';
my $summary=$outdir.'/'."$project\_summaryStatistics.txt";
`mkdir -p $snpdir $gapdir $statdir`;

if ($time==1 && $check==0 && $data!=7){
   open(ALL,">$workdir/working_list.txt")||die "$!";
   open(STAT,">$summary")||die "$!";
}

if ($time==1 && $check>0){open(ALL,">>$workdir/working_list.txt")||die "$!";}
elsif($time==2 || $data==7){
   open (ALL, ">>$workdir/working_list.txt")||die "$!";
   open(STAT,">>$summary")||die "$!";
}

if ($check==0){
   if ($data==0){$nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==1){
      my ($list,$genome_size)=PhaME::prepareComplete($workdir,$reference,$name) if ($time ne "2");
      $fasta_list{$list}=$genome_size if ($list);
      $contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}
   }
   if ($data==2){
      my ($list,$genome_size)=PhaME::prepareComplete($workdir,$reference,$name) if ($time ne "2");
      $fasta_list{$list}=$genome_size if ($list);
      $read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}
   }
   if ($data==3){$nucmer=1;$contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==4){$nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==5){
      my ($list,$genome_size)=PhaME::prepareComplete($workdir,$reference,$name) if ($time ne "2");
      $fasta_list{$list}=$genome_size if ($list);
      $contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}
   }
   if ($data==6){$nucmer=1;$contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   if ($data==7){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
}

if ($time==1){
   if ($check==1){
      print "\n\tContinuing run...\n";
      if ($data<=2){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==3){$contig_nucmer=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==4 || $data==5){$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){$contig_nucmer=1;$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==2){
      print "\n\tContinuing run...\n";
      if ($data<=2){$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data>2 && $data<6){$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){$read_mapping=1;$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==3){
      if ($data<=2){
         if ($pselection>0){print "\n\tContinuing run...\n";$ps=1;}
         else{&print_timeInterval($runtime,"\nRun is complete... Nothing to do here\n");}
      }
      if ($data>2 && $data<6){print "\n\tContinuing run...\n";$buildtree=1;if ($pselection>0){$ps=1;}}
      if ($data==6){print "\n\tContinuing run...\n";$buildSNP=1;$buildtree=1;if ($pselection>0){$ps=1;}}
   }
   if ($check==4){
      if ($data>2 && $data<6){
         if ($bsignal==1){print "\n\tContinuing run...\n";$bs=1;}
         if ($pselection>0){print "\n\tContinuing run...\n";$ps=1;}
         else{&print_timeInterval($runtime,"\nRun is complete... Nothing to do here\n");}
      }
      if ($data==6){print "\n\tContinuing run...\n";$buildtree=1;if ($pselection>0){$ps=1;}}
   }
}
if ($time==2 && $check==1){exit;}
if ($check==0){print "Complete.\n";}

&print_timeInterval($runtime,"\tPreparing files... ");
if ($check==0){
   `mkdir -p $workdir/files`;
   
   if ($time==1){
      opendir(DIR,$refdir);
      while(my $files=readdir(DIR)){
         next if ($files=~/^\..?/);
         my ($rname,$rpath,$rsuffix)=fileparse("$files",qr/\.[^.]*/);
         if ($nucmer==1){
            if ($files=~/.+\.fn|s?a?s?t?a$/){
               my $fasta="$refdir/$files";
               my ($list,$genome_size)=PhaME::prepareComplete($workdir,$fasta,$rname);
               $fasta_list{$list}=$genome_size;
            }
         }
      }
      closedir DIR;
   }

   opendir(DIR,$workdir);
   while(my $files=readdir(DIR)){
      next if ($files=~/^\..?/);
      my ($qname,$qpath,$qsuffix)=fileparse("$files",qr/\.[^.]*/);
      if ($nucmer==1){
         if ($files=~/.+\.fn|s?a?s?t?a$/ && $files!~/contigs?/ && $files!~/fa?s?t?q/){
            my $fasta="$workdir/$files";
            if (!exists $fasta_list{$qname}){
               my ($list,$genome_size)=PhaME::prepareComplete($workdir,$fasta,$qname);
               $fasta_list{$list}=$genome_size;
            }
         }
      }
      if ($contig_nucmer==1){
         if ($files=~ /.+\.contigs?/ && $files!~ /fa?s?t?q/){
            my $contig=$workdir.'/'.$files;
            my ($list,$genome_size)=PhaME::prepareContig($workdir,$contig,$qname);
            $contig_list{$list}=$genome_size;
         }
      }
      if ($read_mapping==1){
         if ($files=~ /.+\.f{1}a?s?t?q$/ && $files!~ /.+\.f{1}n|s?a?s?t?a$/ && $files!~ /\.contigs?/){
            my $fastq=$refdir.'/'.$files;
            print "$qname\n";
            my $read_list_name;
            if ($qname=~/(.+)[_.]R?[12]$/){
               if ($reads==2){
                  $read_list_name=$1.'_pread';
#                   print "$read\n";
               }
               if ($reads==1 && !exists $read_list{"$1_pread"}){
                  $read_list_name=$1.'_sread';
               }
               if ($reads==3){
                  delete $read_list{"$1_sread"};
		  $read_list_name=$1.'_read';
                  $read_list_name=$1.'_pread' if ( ! -e "$workdir/$1.fastq" && ! -e "$workdir/$1.fq");
	       }
                  $read_list{$read_list_name}++;
            }
            else{
               next if ($reads==2 || exists $read_list{"${qname}_read"});
               $read_list_name=$qname.'_sread';
               $read_list{$read_list_name}++;
            } 
         }
      }
   } # workdir
   closedir DIR;
   print "Complete\n";
} # check=0

$name =~ s/\W/_/g;
$reference = "$workdir/files/$name.fna";

if ($nucmer==1){
   open (FAS, ">$workdir/fasta_list.txt")||die "$!";
   open (CON, ">$workdir/contigs_list.txt")||die "$!";

   foreach my $names(sort keys %fasta_list){
      print FAS "$names\n";
      print ALL "$names\n";
      print STAT "$names\tTotal_length\t",$fasta_list{$names},"\n" if ($time ne "2");
#      print "$names\t",$fasta_list{$names},"\n"; 
   }

   &print_timeInterval($runtime,"\tRunning NUCmer on complete genomes\n");
   PhaME::completeNUCmer($workdir,$bindir,"$workdir/fasta_list.txt",$type,$threads,$error,$logfile);
#   &print_timeInterval($runtime,"\tNUCmer on genomes complete\n");
}

if ($contig_nucmer==1){
   open (CON, ">$workdir/contigs_list.txt")||die "$!";
   foreach my $names(sort keys %contig_list){
      print CON "$names\n";
      print ALL "$names\n";
      print STAT "$names\tTotal_length\t",$contig_list{$names},"\n";
   }

   &print_timeInterval($runtime,"Running NUCmer on contigs\n");
   PhaME::contigNUCmer($workdir,$bindir,"$workdir/contigs_list.txt",$code,$threads,$reference,$time,$error,$logfile);
#   PhaME::contigNUCmer($workdir,$bindir,"$workdir/contigs_list.txt",$threads,$reference,"2",$error);
}
close STAT;

if ($read_mapping==1){
   open (READ, ">$workdir/reads_list.txt")||die "$!";
   foreach my $names(keys %read_list){print ALL "$names\n";}
   foreach my $names(keys %read_list){print READ "$names\n";}
   &print_timeInterval($runtime,"Preparing to map reads\n");

   if ($time==2){
      &print_timeInterval($runtime,"Identifying all gaps\n");
      $mappingGaps=PhaME::identifyGaps($outdir,"$workdir/fasta_list.txt",$name,"map",$project);
      PhaME::removeGaps($bindir,$reference,$mappingGaps);
      &print_timeInterval($runtime,"Mapping reads to reference\n");
      my $end=PhaME::readsMapping($workdir,$bindir,"$workdir/reads_list.txt",$threads,$name,$error,$aligner,$logfile);
      &print_timeInterval($runtime,"$end\n");
   }
   elsif ($time==1){
      my $tempdir=$outdir.'/temp/';
      `mkdir -p $tempdir; cp $reference $tempdir`;
      &print_timeInterval($runtime,"Mapping reads to reference\n");
      my $end=PhaME::readsMapping($workdir,$bindir,"$workdir/reads_list.txt",$threads,$name,$error,$aligner,$logfile);
      &print_timeInterval($runtime,"$end\n");
   }
}

if ($buildSNP==1){
   if ($gsignal==1){
      &print_timeInterval($runtime,"Preparing to identify SNPs\n");
      print "\tGFF annotation file provided, SNPs will be differentiated as coding vs noncoding\n";
      my ($genname,$genpath,$gensuffix)=fileparse("$annotation",qr/\.[^.]*/);
      PhaME::codingRegions($outdir,$annotation,$genname);
   }
   &print_timeInterval($runtime,"Identifying SNPs\n");
   PhaME::identifyGaps($outdir,"$workdir/working_list.txt",$name,"snp",$project);
   my $end=PhaME::buildSNPDB($outdir,$bindir,$reference,"$workdir/working_list.txt",$project,$gsignal,$error,$logfile);
   &print_timeInterval($runtime,"$end\n");
}

if ($buildtree==1|| $bs==1){
   if (($tree==2||$tree==3)&&($modeltest==1)){
       my $jmodeltest_jar = "$RealBin/../ext/opt/jmodeltest-2.1.10/jModelTest.jar";
       PhaME::modeltest($jmodeltest_jar,$outdir,$project,$threads,$error,$logfile);
   }
   my $end=PhaME::buildTree($bindir,$outdir,$threads,$tree,"$project\_all",$error,$logfile);
   if ($gsignal==1){PhaME::buildTree($bindir,$outdir,$threads,$tree,"$project\_cds",$error,$logfile);}
   &print_timeInterval($runtime,"$end\n");
   if ($bsignal==1){
      my $end=PhaME::bootstrap($bindir,$outdir,$threads,$tree,"$project\_all",$bootstrap,$error,$logfile);
      &print_timeInterval($runtime,"$end\n");
   }
}

my $tinfo=$refdir."/RAxML_info.$project\_all";
my $tlog=$refdir."/RAxML_log.$project\_all";
my $tparsimony=$refdir."/RAxML_parsimonyTree.$project\_all";
my $tresult=$refdir."/RAxML_result.$project\_all";
my $tbest=$refdir."/RAxML_bestTree.$project\_all";
if (-e $tinfo){`mv $tinfo $outdir/`;}
if (-e $tlog){`mv $tlog $outdir/`;}
if (-e $tparsimony){`mv $tparsimony $outdir/`;}
if (-e $tresult){`mv $tresult $outdir/`;}
if (-e $tbest){`mv $tbest $outdir/`;}
if ($tree==2||$tree==3){$tbest=$outdir."/RAxML_bestTree.$project\_cds";}
if ($tree==1){$tbest=$outdir."/$project\.fasttree";}

my $pamldir;
if ($ps==1){
   my $end=0;
   my $stats_file=$outdir."/$project\_stats.txt";
   my $genedir= $outdir.'/PSgenes';
   if ($pselection==1 || $pselection==3){
      $pamldir= $outdir.'/paml';
      if (!-d $pamldir){`mkdir -p $pamldir`;}
      `cp $tbest $pamldir`;
      if ($tree==2||$tree==3){$ptree= $pamldir."/RAxML_bestTree.$project\_cds";}
      if ($tree==1){$ptree=$pamldir."$outdir./$project\.fasttree";}
      elsif ($tree==0){print "Need to build a tree before any evolutionary analyses can be performed.\n"; exit;}
   }
   my $gapfile=$outdir."/$project\_gaps.txt";

   $end=PhaME::extractGenes($outdir,$stats_file,$reference,$bindir,"$workdir/working_list.txt",$threads,$gapfile,$genefile,$error,$logfile);
   &print_timeInterval($runtime,"$end\n");

   $end=PhaME::translateGenes($outdir,$bindir,$threads,"translate",$error,$logfile);
   &print_timeInterval($runtime,"$end\n");

   $end=PhaME::alignGenes($outdir,$bindir,$threads,"mafft",$error,$logfile);
   &print_timeInterval($runtime,"$end\n");

   $end=PhaME::revTransGenes($outdir,$bindir,$threads,"pal2nal",$error,$logfile);
   &print_timeInterval($runtime,"$end\n");

   my $core=$genedir."/".$project."_core_genome.cdn";
   $end=PhaME::core($outdir,$bindir,$core,$error,$logfile);
   &print_timeInterval($runtime,"$end\n");

   if ($pselection==1 || $pselection==3){
      PhaME::paml($outdir,$bindir,$ptree,0,"Sites","1,2",$core,$threads,$error,$logfile);
      PhaME::paml($outdir,$bindir,$ptree,2,"BrSites",2,$core,$threads,$error,$logfile);
   }
   if ($pselection==2 || $pselection==3){
      my $rootedtree;
      if ($tree==2||$tree==3){$rootedtree= "$outdir/RAxML_rootedTree.$project\_cds_r";}
      if ($tree==1){$rootedtree="$outdir/$project\_cds_rooted.fasttree";}
      print "$rootedtree\n";
      PhaME::hyphy($outdir,$bindir,$tbest,$rootedtree,$core,$threads,"bsrel",$error,$logfile);
   }
}

if ($clean){PhaME::clean($outdir);}

&print_timeInterval($runtime,"End of run!\n");

close STDOUT;
###########################################################
#
# SUBROUTINES
#
###########################################################

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

sub print_timeInterval
{
my $now=shift;
my $msg=shift;

$now=time - $now;

my $string=sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
print "[$string]  $msg";
}

