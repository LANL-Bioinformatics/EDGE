#!/usr/bin/perl -w

use strict;
use FileHandle;
use Getopt::Long;
use Parallel::ForkManager;

my $dir;
my $stat;
my $file;
my $list;
my $thread;
my %sequence;
my $seq;
my $header;
my @strings;
my $count=0;
my ($refid,$queryid,$rpos,$qpos,$snp,$start,$end);
my %coords;
my $fragment;
my @header_list;
my ($first,$second);
my $genefile;
my $first_line=1;
my $outfile;
my @indexArray;
my %reference;
my $name;
my %alternative;
my $gapfile;
my %gaps;
my $prev=0;

GetOptions(
   'd|dir=s'      => \$dir,
   's|stat=s'     => \$stat,
   'f|file=s'     => \$file,
   'l|list=s'     => \$list,
   't|thread=i'   => \$thread,
   'g|gene=s'     => \$genefile,
   'p|gap=s'      => \$gapfile,
);

my $maxthreads = ($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
if ($threads <1 || $threads>$maxthreads){die ("-thread value must be between 1 and $maxthreads.\n");}

my $genedir=$dir.'/PSgenes';
if (!-d $genedir){`mkdir $genedir`;}

my $pm=new Parallel::ForkManager($thread);
$pm->run_on_finish (sub{my ($pid,$ident)=@_;});

open (IN,$gapfile)|| die "$!";
while (<IN>){
   chomp;
   my ($gapstart,$gapend)=split " ",$_;
   my $length=$gapend-$gapstart;
#   print "$temp1\t$gapstart\t$gapend\t$length\t$temp2\n";
   if ($length>=100){$gaps{$gapstart}=$gapend;}
}

my $fh=FileHandle->new($file)|| die "$!";
if ($fh->open("< $file")){
   $/=">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless ($_){next;};
      ($header,@strings)=split /\n/,$_;
      $seq=join "",@strings;
      $reference{$header}=$seq;
      $name=$header;
#      print "$name\n$seq\n";
   }
   $/="\n";
   $fh->close;
}

my $fh1=FileHandle->new($genefile)|| die "$!";
if ($fh1->open("< $genefile")){
   $/=">";
   while (<$fh1>){
      $_=~ s/\>//g;
      unless ($_){next;};
      ($header,@strings)=split /\n/,$_;
      $seq=join "",@strings;
      $sequence{$header}=$seq;
#      print "$header\n$seq\n";
   }
   $/="\n";
   $fh1->close;
}

open (LIST,"$list")||die "$!";
while (<LIST>){
   chomp;
   push(@header_list,"$name:$_");
#   print "$name:$_\n";
}
close LIST;

open (IN,"$stat")||die "$!";
while (<IN>){
   chomp;
   if ($first_line){$first_line=0;}
   else{
      if (/^(\S+)\s+(\S+)\s+coding SNP\s+(\d+)\s+(\d+)\s+.+\s+(\S)\s+(\d+)\s+(\d+)$/){
         ($refid,$queryid,$rpos,$qpos,$snp,$start,$end)=($1,$2,$3,$4,$5,$6,$7);
#         print "$refid, $queryid, $rpos, $qpos, $snp, $start, $end\n";
         if ($start==0){$start=1;}
         if (abs $start-$end> 300){
#            print "$start\t$end\n";
            $coords{$start}=$end;
            $alternative{$rpos}{"$refid:$queryid"}=$snp;
         }
      }
   }
}
close IN;

my $temp=keys %coords;
my $size = length $temp;
my $index_file=$dir.'/gene_index.txt';
my $gapgenes=$dir.'/gene_gaps.txt';
my $compliment;

open (GAP,">$gapgenes")||die "$!";
OUTER:foreach my $start(sort{$a<=>$b}keys %coords ){
   my $end=$coords{$start};

   foreach my $gap (keys %gaps){
      if ($gap>=$start && $gaps{$gap}<=$end){
         if ($gaps{$gap}-$gap %3 != 0){
            print GAP "Gene with coords $start\_$end has a gap from $gap to $gaps{$gap} in the middle.\n";
            next OUTER;
         }
      }
      if ($gap<$start && $gaps{$gap}>$start && $gaps{$gap}<$end){
         if ($gaps{$gap}-$start %3 != 0){
            print GAP "Gene with coords $start\_$end has a gap from $start to $gaps{$gap} in the beginning.\n";
            next OUTER;
         }
      }
      if ($gap>$start && $gap<$end && $gaps{$gap}>$end){
         if ($end-$gap %3 != 0){
            print GAP "Gene with coords $start\_$end has a gap from $gap to $end at the end.\n";
            next OUTER;
         }
      }
   }

   $count++;
   $outfile=$genedir.'/Gene';
   $outfile.=sprintf "%0${size}d", $count;
   my $index=sprintf "%0${size}d", $count;
   $outfile.='.fna';
#   print "$index\n";
   push (@indexArray,"$index\t$start\t$end");
   $pm->start and next;
   open (OUT,">$outfile")||die "$!";
#   print "$start\t$end\n";

   foreach my $comparison(@header_list){
      if ($comparison=~/(.+):(.+)/){
         ($first,$second)=($1,$2);
         print OUT ">$second\n";
#          print "$comparison\t$first\n";
      }

      for my $entry(keys %sequence){
         my $gene=$sequence{$entry};
#         print "$entry\n";
         if ($entry=~ /.+\_$start\_$end\_(-?1)$/){
            my $strand=$1;
            if ($strand<0){
               $compliment=reverse($gene);
               $compliment =~ tr/ACGTacgt/TGCAtgca/;
               $gene=$compliment;
            }
#            print "$first\t$start\t$end\n";
            foreach my $position(sort {$a<=>$b}keys %alternative){
               if (defined $alternative{$position}{$comparison}){
                  if ($position>=$start && $position<=$end){
                     my $snp=$position-$start; 
#                     print "$start\t$end\t$position\t$snp\t",$alternative{$position}{$comparison},"\n";
                     substr($gene,$snp,1,$alternative{$position}{$comparison});
                  }
               }
            }
            if ($strand<0){
               $compliment=reverse($gene);
               $compliment =~ tr/ACGTacgt/TGCAtgca/;
               $gene=$compliment;
            }
            print OUT "$gene\n";
#            print "$gene\n";
         }
      }
   }
   close OUT;
   $pm->finish;
}
$pm->wait_all_children;

open (IND,">$index_file")||die "$!";
foreach my $index(@indexArray){print IND $index,"\n";}
close IND;

