#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;
use FindBin qw($RealBin);

# set up environments
$ENV{PATH}="$RealBin:$RealBin/../ext/bin:$ENV{PATH}";

my $dir;
my $thread;
my $tree;
my $analysis;
my $seqfile;
my @geneList;
my %models;
my $tree_output;
my $core;

GetOptions(
   'i|dir=s'       => \$dir,
   't|thread=i'    => \$thread,
   'r|tree=s'      => \$tree,
   'a|analysis=s'  => \$analysis,
   'o|output=s'    => \$tree_output,
   'c|core=s'      => \$core,
);


if ($dir=~ /.+\/$/){my $temp= chop($dir);}
my ($name,$path,$suffix)=fileparse("$tree",qr/\.[^.]*/);
if ($suffix=~ /.+\/$/){my $temp= chop($suffix);}

=headmy $Rscript = "$dir/Rscript$$";
open (Rscript, ">$Rscript");
print Rscript "
library(ape)
library(phangorn)

tree <- read.tree(\"$tree\")
rtree <- midpoint(tree)
write.tree(rtree, file = \"$tree_output\")
";
close Rscript;
system ("R --vanilla --slave --silent < $Rscript 2>/dev/null");
=cut

my $genedir=$dir."/PSgenes";
my $hyphydir=$dir."/hyphy";
if (!-d $hyphydir){`mkdir $hyphydir`;}
my $modelsFile=$hyphydir.'/models.txt';
my $hyphy="$RealBin/../ext/lib/hyphy/TemplateBatchFiles";
chdir $hyphy;

read_directory($genedir);

if ($analysis=~/model/i || $analysis=~/all/i){chooseModel();}
if (-e $modelsFile){
   open (IN,"$modelsFile")||die "Can't run SLAC analysis without models.\n";
   while (<IN>){
      chomp;
      my ($gene,$model)=split /\s+/,$_;
      $models{$gene}=$model;
   }
   close IN;
}
elsif ($analysis!~/BSrel/i && !-e $modelsFile){print "Can't run subsequent analyses (SLAC, REL, FEL) without models.\n"; exit;}
if ($analysis=~/^SLAC/i || $analysis=~/^all/i){runSLAC();}
if ($analysis=~/^FEL/i || $analysis=~/^all/i){runFEL();}
if ($analysis=~/^REL/i || $analysis=~/^all/i){runREL();}
if ($analysis=~/^BSrel/i | $analysis=~/^all/i){runBranchSiteREL();}

sub read_directory
{
my $dir=shift;

opendir(DIR,$dir);
foreach my $file(sort{$a cmp $b}readdir(DIR)){
   if ($file=~ /\.cdn$/){
      $seqfile=$genedir.'/'.$file;
      push(@geneList,$seqfile);
   }
}
closedir DIR;
}

sub chooseModel
{
my $models=$hyphydir.'/models.txt';
if (-e $models){`rm $models`;}

my $pm=new Parallel::ForkManager($thread);
$pm->run_on_finish (
   sub{
      my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$file)=@_;
      if (defined($file)){
         my $result=${$file};

         my ($name,$path,$suffix)=fileparse("$result",qr/\.[^.]*/);
         $name=~ s/_model//;

         open (OUT, ">>$models");
         open (IN, "$result")||die "$!";
         while (<IN>){
            chomp;
            if (/AIC based winner: \((\d+)\).+/){
               print OUT "$name\t$1\n";
            }
         }
      }
   }
);

for (my $i=0;$i<=$#geneList; $i++){
   $pm->start($i) and next;
   my ($name,$path,$suffix)=fileparse("$geneList[$i]",qr/\.[^.]*/);
   my $batch=$hyphydir.'/'.$name.'_model.bf';
   my $output=$hyphydir.'/'.$name.'_model.txt';

   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"NucModelCompare.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"Local\";
inputRedirect[\"02\"]=\"$geneList[$i]\";
inputRedirect[\"03\"]=\"$tree_output\";
inputRedirect[\"04\"]=\"0.05\";
inputRedirect[\"05\"]=\"No\";
inputRedirect[\"06\"]=\"$output\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";} 

   $pm->finish(0,\$output);
}
$pm->wait_all_children;
}

sub runSLAC
{
my $pm=new Parallel::ForkManager($thread);

for (my $i=0;$i<=$#geneList; $i++){
   $pm->start($i) and next;
   my ($name,$path,$suffix)=fileparse("$geneList[$i]",qr/\.[^.]*/);
   my $model=$models{$name};

   my $batch=$hyphydir.'/'.$name.'_slac.bf';
   my $nwkOutput=$hyphydir.'/'.$name.'_slac.nwk';
   my $output=$hyphydir.'/'.$name.'_slac.tab';
   
   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"QuickSelectionDetection.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"Universal\";
inputRedirect[\"02\"]=\"New Analysis\";
inputRedirect[\"03\"]=\"$geneList[$i]\";
inputRedirect[\"04\"]=\"Custom\";
inputRedirect[\"05\"]=\"$model\";
inputRedirect[\"06\"]=\"$tree_output\";
inputRedirect[\"07\"]=\"$nwkOutput\";
inputRedirect[\"08\"]=\"Estimate dN/dS only\";
inputRedirect[\"09\"]=\"Single Ancestor Counting\";
inputRedirect[\"10\"]=\"Full tree\";
inputRedirect[\"11\"]=\"Averaged\";
inputRedirect[\"12\"]=\"Approximate\";
inputRedirect[\"13\"]=\"0.05\";
inputRedirect[\"14\"]=\"Export to File\";
inputRedirect[\"15\"]=\"$output\";
inputRedirect[\"16\"]=\"Count\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";}

   $pm->finish();
}
$pm->wait_all_children;
}

sub runFEL
{
my $pm=new Parallel::ForkManager($thread);

for (my $i=0;$i<=$#geneList; $i++){
   $pm->start($i) and next;
   my ($name,$path,$suffix)=fileparse("$geneList[$i]",qr/\.[^.]*/);
   my $model=$models{$name};

   my $batch=$hyphydir.'/'.$name.'_fel.bf';
   my $nwkOutput=$hyphydir.'/'.$name.'_fel.nwk';
   my $output=$hyphydir.'/'.$name.'_fel.lrt';

   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"QuickSelectionDetectionMF.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"Universal\";
inputRedirect[\"02\"]=\"New Analysis\";
inputRedirect[\"03\"]=\"Custom\";
inputRedirect[\"04\"]=\"$model\";
inputRedirect[\"05\"]=\"1\";
inputRedirect[\"06\"]=\"$geneList[$i]\";
inputRedirect[\"07\"]=\"$tree_output\";
inputRedirect[\"08\"]=\"$nwkOutput\";
inputRedirect[\"09\"]=\"Estimate dN/dS only\";
inputRedirect[\"10\"]=\"FEL\";
inputRedirect[\"11\"]=\"0.05\";
inputRedirect[\"12\"]=\"All\";
inputRedirect[\"13\"]=\"$output\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";}

   $pm->finish();
}
$pm->wait_all_children;
}

sub runREL
{
my $pm=new Parallel::ForkManager($thread);

for (my $i=0;$i<=$#geneList; $i++){
   $pm->start($i) and next;
   my ($name,$path,$suffix)=fileparse("$geneList[$i]",qr/\.[^.]*/);
   my $model=$models{$name};

   my $batch=$hyphydir.'/'.$name.'_rel.bf';
   my $output=$hyphydir.'/'.$name.'_rel.txt';

   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"dNdSRateAnalysis.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"$geneList[$i]\";
inputRedirect[\"02\"]=\"Universal\";
inputRedirect[\"03\"]=\"$tree_output\";
inputRedirect[\"04\"]=\"Codon Model\";
inputRedirect[\"05\"]=\"MG94xCustom\";
inputRedirect[\"06\"]=\"$model\";
inputRedirect[\"07\"]=\"Run All\";
inputRedirect[\"08\"]=\"Independent Discrete\";
inputRedirect[\"09\"]=\"Default\";
inputRedirect[\"10\"]=\"3\";
inputRedirect[\"11\"]=\"3\";
inputRedirect[\"12\"]=\"$output\";
inputRedirect[\"13\"]=\"None\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";}

   $pm->finish();
}
$pm->wait_all_children;
}

sub runBranchSiteREL
{
my $pm=new Parallel::ForkManager($thread);

=head
if (defined $core){
   $pm->start($core) and next;
   my ($name,$path,$suffix)=fileparse("$core",qr/\.[^.]*/);

   my $batch=$hyphydir.'/'.$name.'_BSrel.bf';
   my $output=$hyphydir.'/'.$name.'_BSrel';

   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"BranchSiteREL.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"Universal\";
inputRedirect[\"02\"]=\"Yes\";
inputRedirect[\"03\"]=\"No\";
inputRedirect[\"04\"]=\"$core\";
inputRedirect[\"05\"]=\"$tree_output\";
inputRedirect[\"06\"]=\"All\";
inputRedirect[\"07\"]=\"\";
inputRedirect[\"08\"]=\"$output\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";}

   $pm->finish();
}
=cut

for (my $i=0;$i<=$#geneList; $i++){
   $pm->start($i) and next;
   my ($name,$path,$suffix)=fileparse("$geneList[$i]",qr/\.[^.]*/);

   my $batch=$hyphydir.'/'.$name.'_BSrel.bf';
   my $output=$hyphydir.'/'.$name.'_BSrel';

   open (OUT, "> $batch")||die "$!";
   print OUT "
fileToExecute = HYPHY_BASE_DIRECTORY + DIRECTORY_SEPARATOR + \"BranchSiteREL.bf\"; 

inputRedirect={};
inputRedirect[\"01\"]=\"Universal\";
inputRedirect[\"02\"]=\"Yes\";
inputRedirect[\"03\"]=\"No\";
inputRedirect[\"04\"]=\"$geneList[$i]\";
inputRedirect[\"05\"]=\"$tree_output\";
inputRedirect[\"06\"]=\"All\";
inputRedirect[\"07\"]=\"\";
inputRedirect[\"08\"]=\"$output\";

ExecuteAFile (fileToExecute,inputRedirect);
";
   my $command="HYPHYMP $batch";
   if (system($command)){die "Error running $command.\n";}

   $pm->finish();
}
$pm->wait_all_children;
}

