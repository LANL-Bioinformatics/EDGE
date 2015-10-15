#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

my $dir;
my $program="mafft";
my $muscle_options= "-diags";
my $mafft_options= "";
my @files;
my $thread;

GetOptions(
   'd=s'      => \$dir,
   't=i'      => \$thread,
   'm=s'      => \$program,
   'h|help'   => sub {usage()},
);

if (!$dir) {die &usage;}
if ($dir=~ /.+\/$/){my $tmp= chop($dir);}

sub usage
{
print STDERR"
Usage: $0 -d (directory) -m (script)
   required options
     -d                 File directory
     -m                 program to run
                          default: mafft
     -t                 Number of threads
   additional options
     -muscle <string>   muscle options 
                          default: \"-diags\"
                          type \"muscle\" to see additional muscle options 
     -mafft <string>    mafft options
                          type \"mafft\" to see additional mafft options 
   list of progams 
     -mafft:      Multiple alignment program
     -muscle:     Multiple alignment program
     -translate:  Translate DNA to amino acid
     -oneline:     Write multi-line fasta sequences to one line
     -pal2nal:    Convert amino acid alignment to codon alignment, needs DNA fasta file
";
exit;
} 
my $bindir=getBinDirectory();
my $pm= new Parallel::ForkManager($thread);
$pm->run_on_finish(sub{my ($pid,$ident)=@_;});

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

opendir(DIR, $dir);
while (my $file= readdir(DIR)){
   next if ($file=~ /^..?$/);
   $file= $dir.'/'.$file;
   if ($program=~ /pal2nal/ ||$program=~ /oneline/){
      if ($file=~ /mft$/){push (@files,$file);}
   }
   if ($program=~ /translate/){
      if ($file=~ /fna$/){push (@files,$file);}
   }
   if ($program=~ /mafft/ ||$program=~ /muscle/){
      if ($file=~ /faa$/){push (@files,$file);}
   }
}

for (my $i=0; $i<=$#files; $i++){
   $pm->start and next;   
   my ($ref_file_name,$ref_file_path,$ref_file_suffix)=fileparse("$files[$i]", qr/\.[^.]*/);
   if ($program=~/oneline/){`$bindir/fasta_oneline.pl $files[$i] > $ref_file_path/$ref_file_name.faa`}
   if ($program=~ /mafft/i){`mafft $mafft_options $files[$i] > $ref_file_path$ref_file_name.mft`;}
   if ($program=~/pal2nal/){`pal2nal.pl $files[$i] $ref_file_path/$ref_file_name.fna -output paml > $ref_file_path/$ref_file_name.cdn`;}
#   if ($program=~/pal2nal/){`$bindir/pal2nal.pl $files[$i] $ref_file_path/$ref_file_name.fna -output fasta > $ref_file_path/$ref_file_name.cdn`;}
   if ($program=~ /translate/){`$bindir/translate.pl $files[$i] > $ref_file_path/$ref_file_name.faa`;}
   elsif ($program=~ /muscle/i){`muscle3.8.31_i86linux64 $muscle_options -in $files[$i] -out $ref_file_path$ref_file_name.mft`;}
$pm->finish;
}
$pm->wait_all_children;

