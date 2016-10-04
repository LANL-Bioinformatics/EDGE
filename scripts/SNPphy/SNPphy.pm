#!/usr/bin/env perl

use strict;
use warnings;
use FileHandle;

package SNPphy;

# Check to see if this really is a first time run and if the previous run finished properly 
sub check
{
my $dir=shift;
my $refdir=shift;
my $time=shift;
my $data=shift;
my $wdir=$dir.'/results';
my $log=$wdir.'/currentRun.log';
my $list=$dir.'/working_list.txt';
my $progress=0;
my $alignment=0;
my %refcheck;
=head
Another check to make sure that the first run finished properly
if (-d $wdir){ 
   opendir (DIR,$wdir) or die "$!";
   if(grep ! /^\.\.?$/,readdir DIR){print 'stuff in here';}
}
=cut
if (-e $log && !-z $log){
   open (LOG,"$log")||die "$!";
   while (<LOG>){ 
      chomp;
      if (/NUCmer on all reference genomes complete/){$progress+=1;}
      if (/NUCmer on all contig\/draft genomes complete/){$progress+=1;}
      if (/Read Mapping complete/){$progress+=1;}
      if (/SNP alignment complete/){$progress+=1;$alignment=1;}
      if (/Tree phylogeny complete/ && $alignment==1){$progress+=1;}
   }
}
else {return("0");}

if ($time==1){
   if ($data<=2){
      if ($progress==3){return($progress);}
      if ($progress<3){
         print "\nWarning: Previous run not complete.\n";
         return($progress);
      }
   }
   if ($data>=3 && $data<=5){
      if ($progress==4){return($progress);}
      if ($progress<4){
         print "\nWarning: Previous run not complete.\n";      
         return($progress);
      }
   }
   if ($data==6){
      if ($progress==5){return($progress);}
      if ($progress<5){
         print "\nWarning: Previous run not complete.\n";
         return($progress);
      }
   }
}
#if ($time==2 && $progress==0){
#   print "\nWarning: Previous run not complete.\nSNP and gap coordinate files are needed to update SNP alignment.\n";
#   return ("1");
#}
if ($time==2){
   open (LIST,$list);
   while (<LIST>){chomp;$refcheck{$_}++;}
   close LIST;
   opendir(DIR,"$refdir");
   while (my $file= readdir(DIR)){
      next if ($file=~ /^..?$/);
      
      if ($file=~/(.+)\.fn|s?a?s?t?a$/){
         $file=$refdir.'/'.$file;
         if(!exists $refcheck{$1}){`cp $file $dir/files/`;}
      }
      #print "$file\n";
   }

return 0;
}# check for second run and snp_alignment available. Can be done after running buildSNPdb.pl
}

# Concatenate multiple reference chromosomes from one genome into one mega chromosome 
sub prepareComplete
{
my $rdir=shift;
my $wdir=shift;
my $file=shift;
my $name=shift;

open (OUT,">$wdir/files/$name.fna")||die "$!";
print OUT ">$name\n";
open (IN,$file)||die "$!";
while (<IN>){
   chomp;
   if (!/^>/){print OUT $_;}
}
print OUT "\n";
close OUT;
close IN;
return $name;
}

# Change contigs names
sub prepareContig
{
my $dir=shift;
my $file=shift;
my $name=shift;
my ($header,@seq);
my $sequence;
my $count=1;
my $contig=$name.'_contig';
my $outfile=$dir.'/files/'.$name.'_contigs.fna';

open (OUT, ">$outfile")||die "$!";
my $fh=FileHandle->new($file)|| die "$!";
if ($fh->open("<$file")){
   $/=">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless($_){next;};
      ($header,@seq)=split /\n/,$_;
      $sequence= join "",@seq;
      $header= $name.'_'.$count;
#      print ">$header\n$sequence\n";
      print OUT ">$header\n$sequence\n";
      $count++;  
   }
   $/="\n";
   $fh->close;
   close OUT;
}
return $contig;
}

# Identifies gap coords in reference genomes
# Gaps identified in NUCmer
sub identifyGaps
{
my $dir=shift;
my $list=shift;
my $name=shift;
my $type=shift;
my $gapdir=$dir.'/gaps';
my $repeatdir=$dir.'/stats';
my %query;
my $line=0;
my $gapfile;
my $gap_start;
my $gap_end;

if ($type=~/map/){$gapfile=$dir.'/mapping_gaps.txt';}
elsif ($type=~/snp/){$gapfile=$dir.'/all_gaps.txt';}
open (GAP,">>$gapfile")||die "$!";

open (LIST,"$list")||die "$!";
while (<LIST>){
   chomp;
   $_ =~ s/(.+_read)\.\w+$/$1/;
   $query{$_}++;
}
close LIST;

opendir(DIR,"$gapdir");
while (my $gaps= readdir(DIR))
{
#   if ($gaps=~ /^$name\_(.+)\_norepeats\.gaps/ || $gaps=~ /^$name\_(.+_contig)s\.gaps/ && $gaps!~ /^$name\_norepeats/){
   if ($gaps=~ /^$name\_norepeats\_(.+)\_norepeats\.gaps/ || $gaps=~ /^$name\_(.+_contig)s\.gaps/){
      if (exists $query{$1}){
         $line=0;
         my $gapfile= "$gapdir/$gaps";
         open (IN,$gapfile)|| die "$!";
         #print "Nucmer Gaps\n";
         while (<IN>){$line++; print GAP "$_";}
         close IN;
         if ($line ==0){print "Empty File: $gapfile\n"; $line=0;}
      }
   }
   if ($type=~/snp/)
   {
#      if ($gaps=~ /^$name\_norepeats\_(.+)\_?$name?\_?[\d+\_\d+]?\.gaps$/){
      if ($gaps=~ /^$name\_(.+)\_$name(\_\d+\_\d+)?\.gaps$/)
     #  if ($gaps=~ /^$name\_(.+)\.gaps$/)
       {
         my $query=$1;
         my $tmp= $query."_read";
         if (exists $query{$tmp})
         {
            $line=0;
            my $gapfile= "$gapdir/$gaps";
            open (IN,$gapfile)|| die "$!";
#         print "Read Mapping Gaps\n";
            while (<IN>)
            {
               chomp;
               $line++;
               next if (/Start\s+End\s+Length.+/);
               my ($start,$end,$length,$ref)= split "\t",$_;
               if ($ref=~ /$name\_(\d+)\_(\d+)$/){
                  $gap_start=$start+$1-1;
                  $gap_end=$1+$end-1;
               }
               else{
                  $gap_start=$start;
                  $gap_end=$end;
               }
               print GAP "$name\t$gap_start\t$gap_end\t$length\tREADS_$query\n";
            }
            close IN;
            if ($line == 1){`rm $gapfile`; $line=0;}
         }
      }
   }
#   last OUTER;
}
closedir DIR;

opendir(REPEAT,"$repeatdir");
while (my $repeats= readdir(REPEAT)){
   if ($repeats=~ /$name\_repeat_coords\.txt/){
#      print "$repeats\n";
#   if ($gaps=~ /repeats_coords\.txt/){
      my $repeatfile= "$repeatdir/$repeats";
      open (IN,$repeatfile)|| die "$!";
#      print "Repeats\n";
      while (<IN>){
         chomp;
         my ($ref,$temp,$start,$end,$length)= split "\t",$_;
         print GAP "$ref\t$start\t$end\t$length\tREPEATS\n";
      }
      close IN;
   }
#   last OUTER;
}
closedir REPEAT;
return $gapfile; 
}

# Identify CDS coords 
sub genbank
{
my $dir=shift;
my $genbank=shift;
my $name=shift;
my $start;
my $end;
my $gap_start=1;
my $gap_end;
my $source_start=1;
my $source_end=0;
my %CDS;
my $outfile=$dir."/noncoding.txt";
my $coding=$dir."/CDScoords.txt";
my $line;
my $temp;

open (OUT,">$outfile");
open (CDS,">$coding");

my $first=1;
my $permutation=0;

open (IN, "$genbank")|| die "$!";
while(<IN>){
   chomp;
    if (/##sequence-region/){
      $permutation=$source_end;
      ($line,$temp,$source_start,$source_end)=split " ",$_;
   }
   if (!/^#/){
      my ($name,$source,$method,$start,$stop,$score,$strand,$phase,$field)=split "\t",$_;
      my @group=split ";",$field;

      if ($method=~/CDS/){
         $start=$start+$permutation;
         $stop=$stop+$permutation;

         print CDS "$name\t$start\t$stop\t";
         $CDS{$start}=$stop;
         foreach (@group){if (/product=(.+)/){print CDS $1,"\n";}}
      }
   }

}
my $prev=0;
foreach my $begin (sort{$a<=>$b} keys %CDS){
   my $end=$CDS{$begin};
   if ($first){
      if ($begin==1){$gap_start=$end+1;}
      else{
         $gap_end=$begin-1;
         if ($gap_start<$gap_end){print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";}
      }
      $first =0;
   }
   else {
      $gap_end=$begin-1;
      if ($gap_start<$gap_end){print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";}
   }
   if ($begin==$prev){
      $gap_start=$prev-1;
      $gap_end=$prev-1;
      print OUT "$name\t$gap_start\t$gap_end\tnoncoding\n";
   }

   $gap_start=$end+1;
   $prev=$end+2;
}
}

sub clean
{
my $dir=shift;
`rm -f $dir/*.pileup $dir/*bam $dir/*bcf $dir/*mgaps $dir/*ntref $dir/*sam $dir/*delta $dir/*gapfilter $dir/*snpfilter`;
}

# Run NUCmer on reference genomes 
sub completeNUCmer
{
my $indir=shift;
my $bindir=shift;
my $reference=shift;
my $list=shift;
my $thread=shift;
my $error=shift;
my $outdir=$indir.'/results';
my $log=$outdir.'/currentRun.log';

print "\n";
my $nucmer="time $bindir/runNUCmer.pl -r $reference -q $indir -d $outdir -t $thread -l $list 2>$error > $log\n\n";
print $nucmer;
if (system ($nucmer)){die "Error running $nucmer.\n";}
}

sub contigNUCmer
{
my $indir=shift;
my $bindir=shift;
my $list=shift;
my $thread=shift;
my $reference=shift;
#my $type=shift;
my $error=shift;
my $outdir=$indir.'/results';
my $log=$outdir.'/currentRun.log';

print "\n";
#my $con_nucmer="time $bindir/runContigNUCmer.pl -r $reference -q $indir -d $outdir -l $list -t $thread -y $type 2>>$error >> $log\n\n";
my $con_nucmer="time $bindir/runContigNUCmer.pl -r $reference -q $indir -d $outdir -l $list -t $thread 2>>$error >> $log\n\n";

print $con_nucmer;
if (system ($con_nucmer)){die "Error running $con_nucmer.\n";}
}

# Removes gaps using a gap file
sub removeGaps
{
my $bindir=shift;
my $reference=shift;
my $readgaps=shift;

print "\n";
my $remove="time $bindir/removeGaps.pl $reference $readgaps\n\n";
print $remove;
if (system ($remove)){die "Error running $remove.\n";}
}

# Runs bowtie on paired-end reads
sub readsMapping
{
my $indir=shift;
my $bindir=shift;
my $list=shift;
my $thread=shift;
my $name=shift;
my $error=shift;
my $outdir=$indir."/results";
my $log=$outdir.'/currentRun.log';
my $reference= $outdir.'/temp/'.$name.'.fna';

if(!-e $reference){$reference=$indir.'/files/'.$name.'.fna';} 
print "\n";
my $map="time $bindir/runReadsMapping.pl -r $reference -q $indir -d $outdir -t $thread -l $list -a bowtie 2>>$error >> $log\n\n";
print $map;
if (system ($map)){die "Error running $map.\n";}
}

sub buildSNPDB
{
my $outdir=shift;
my $bindir=shift;
my $reference=shift;
my $list=shift;
my $outfile=shift;
my $signal=shift;
my $error=shift;
my $log=$outdir.'/currentRun.log';

print "\n";
my $SNPdb="time $bindir/buildSNPDB.pl -i $outdir -r $reference -l $list -o $outfile -g $signal 2>>$error >> $log\n\n";
print $SNPdb;
if (system ($SNPdb)){die "Error running $SNPdb.\n";}
}

sub buildTree
{
my $outdir=shift;
my $thread=shift;
my $file=shift;
my $tree=shift;
my $name=shift;
my $error=shift;
my $log=$outdir.'/currentRun.log';

my @type=("cds","all");
foreach my $type (@type)
{
	next if ( ! -e "$outdir/${file}_${type}.aln.fasta");
    if ($tree==1||$tree==3){
       print "\n";
       #my $fasttree="time FastTreeMP -nt -gtr < $outdir/snp_alignment_all_$file > $outdir/snp_alignment_all_$file.fasttree 2>>$error >> $log\n\n";
       my $fasttree="export OMP_NUM_THREADS=$thread; time FastTreeMP -nt -gtr < $outdir/${file}_${type}.aln.fasta > $outdir/FastTree.${type} 2>>$error\n\n";
       print $fasttree;
       if (system ($fasttree)){die "Error running $fasttree.\n";}
    }
    if ($tree==2||$tree==3){
       print "\n";
       my $raxml="time raxmlHPC-PTHREADS -p 10 -T $thread -m GTRGAMMAI -s $outdir/${file}_${type}.aln.fasta -w $outdir -n ${type} 2>>$error >> $log\n\n";
       print $raxml;
       if (system ($raxml)){die "Error running $raxml.\n";}
    }
}
open (OUT, ">>$log");
print OUT "Tree phylogeny complete.\n";
close OUT;
}

sub extractGenes
{
my $dir=shift;
my $stat=shift;
my $file=shift;
my $bindir=shift;
my $list=shift;
my $thread=shift;
my $gapfile=shift;
my $genefile=shift;
my $error=shift;
my $log=$dir.'/currentRun.log';
my $genedir=$dir.'/PSgenes';

print "\n";
my $extract="time $bindir/extractGenes.pl -d $dir -t $thread -l $list -s $stat -f $file -p $gapfile -g $genefile 2>>$error >> $log\n\n";
print $extract;
if (system ($extract)){die "Error running $extract.\n";}

opendir (DIR,"$genedir") ||die "$!";
OUTER:while (my $files= readdir(DIR)){
   next if ($files=~/^..?$/);
   if ($files=~ /.+\.fna$/){
      my $file= $genedir.'/'.$files;
      open (IN,"$file")|| die"$!";
      while (<IN>){
         if (!/^>/){
            if (!/^ATG/){
               print "$file\n"; 
               `rm $file`;
               next OUTER;
            }
         }
      }
      close IN;
   }
}
close DIR;
}

sub translateGenes
{
my $dir=shift;
my $bindir=shift;
my $thread=shift;
my $program=shift;
my $error=shift;
my $log=$dir.'/currentRun.log';
my $genedir=$dir.'/PSgenes';

print "\n";
my $translate="time $bindir/parallel_run.pl -d $genedir -t $thread -m $program 2>>$error >> $log\n\n";
print $translate;
if (system ($translate)){die "Error running $translate.\n";}
}

sub alignGenes
{
my $dir=shift;
my $bindir=shift;
my $thread=shift;
my $program=shift;
my $error=shift;
my $log=$dir.'/currentRun.log';
my $genedir=$dir.'/PSgenes';

my $align="time $bindir/parallel_run.pl -d $genedir -t $thread -m $program 2>>$error >> $log\n\n";
print $align;
if (system ($align)){die "Error running $align.\n";}
}

sub revTransGenes
{
my $dir=shift;
my $bindir=shift;
my $thread=shift;
my $program=shift;
my $error=shift;
my $log=$dir.'/currentRun.log';
my $genedir=$dir.'/PSgenes';

my $revTrans="time $bindir/parallel_run.pl -d $genedir -t $thread -m $program 2>>$error >> $log\n\n";
print $revTrans;
if (system ($revTrans)){die "Error running $revTrans.\n";}
}

sub positiveSelection
{
my $dir=shift;
my $bindir=shift;
my $tree=shift;
my $model=shift;
my $suffix=shift;
my $NSsites=shift;
my $thread=shift;
my $error=shift;
my $log=$dir.'/currentRun.log';
my $branch=1;
my $pamldir=$dir.'/paml';

if ($model==0){
   print "\n";
   my $ps="time $bindir/runPAML.pl -i $dir -t $thread -r $tree -m $model -n $NSsites -s $suffix 2>>$error >> $log\n\n";
   print $ps;
   if (system ($ps)){die "Error running $ps.\n";} 
   
   `mv $pamldir/*/*$suffix $pamldir`;
   print "\n";
   my $parse="time $bindir/parseSitePAML.pl $pamldir $NSsites 2>>$error >> $log\n\n";
   print $parse;
   if (system($parse)){die "Error running $parse. \n";}
}

if ($model==2){
   print "\n";
   my $edit="time $bindir/ParseTree.pl $tree 2>>$error >> $log\n\n";
   print $edit;
   if (system ($edit)){die "Error running $edit.\n";}

   print "\n";
   my $ps="time $bindir/runPAML.pl -i $dir -t $thread -r $tree -m $model -n $NSsites -s $suffix 2>>$error >> $log\n\n";
   print $ps;
   if (system ($ps)){die "Error running $ps.\n";}

   `mv $pamldir/*/*$suffix $pamldir`;
   print "\n";
   my $parse="time $bindir/parseSitePAML.pl $pamldir 0,1,2,7,8,$NSsites 2>>$error >> $log\n\n";
   print $parse;
   if (system($parse)){die "Error running $parse. \n";}

}
}

1;
