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
my %positions;
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
my $string=0;
my $frame;

GetOptions(
   'd|dir=s'      => \$dir,
   's|stat=s'     => \$stat,
   'f|file=s'     => \$file,
   'l|list=s'     => \$list,
   't|thread=i'   => \$thread,
   'g|gene=s'     => \$genefile,
   'p|gap=s'      => \$gapfile,
);

my $genedir=$dir.'/PSgenes';
if (!-d $genedir){`mkdir $genedir`;}
#print "$outfile\n";

my $pm=new Parallel::ForkManager($thread);
$pm->run_on_finish (sub{my ($pid,$ident)=@_;});

open (IN,$gapfile)|| die "$!";
while (<IN>){
   chomp;
   my ($temp1,$gapstart,$gapend,$length,$temp2)=split "\t",$_;
#   print "$temp1\t$gapstart\t$gapend\t$temp2\t$temp3\n";
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
#      print "$seq\n";
      $reference{$header}=$seq;
      $name=$header;
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
#      print "$seq\n";
      $sequence{$header}=$seq;
   }
   $/="\n";
   $fh1->close;
}

open (LIST,"$list")||die "$!";
while (<LIST>){
   chomp;
   push(@header_list,"$name:$_");
}
close LIST;

open (IN,"$stat")||die "$!";
while (<IN>){
   chomp;
   if ($first_line){$first_line=0;}
   else{
      if (/^(\S+)\s+(\S+)\s+coding SNP\s+(\d+)\s+(\d+)\s+.+\s+(\S)\s+(\d+)\s+(\d+)$/){
         ($refid,$queryid,$rpos,$qpos,$snp,$start,$end)=($1,$2,$3,$4,$5,$6,$7);
#         print "$refid, $queryid, $pos, $start, $end\n";
         if ($start==0){$start=1;}
         if (abs $start-$end> 300){
            $coords{$start}=$end;
            if ($start != $prev){$positions{$start}{"$refid:$queryid"}=$qpos;$prev=$start;}
            $alternative{$rpos}{"$refid:$queryid"}=$snp;
         }
      }
   }
}
close IN;

my $temp=keys %coords;
my $size = length $temp;
my $index_file=$dir.'/gene_index.txt';
my $nohit;
my $gapgenes=$dir.'/gene_gaps.txt';

open (GAP,">$gapgenes")||die "$!";
OUTER:foreach my $start(sort{$a<=>$b}keys %coords ){
   my $end=$coords{$start};
   foreach my $gap (keys %gaps){
      if ($gap>=$start && $gaps{$gap}<=$end){
         print GAP "Gene $start\t$end has a gap from $gap to $gaps{$gap} in the middle.\n";
         next OUTER;
      }
      if ($gap<$start && $gaps{$gap}>$start && $gaps{$gap}<$end){
         if ($gaps{$gap}-$start >10){
            print GAP "Gene $start\t$end has a gap from $start to $gaps{$gap} in the beginning.\n";
            next OUTER;
         }
      }
      if ($gap>$start && $gap<$end && $gaps{$gap}>$end){
         if ($end-$gap >10){
            print GAP "Gene $start\t$end has a gap from $gap to $end at the end.\n";
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
      if ($comparison=~/(.+):(.+)/){($first,$second)=($1,$2);}
#      print "$comparison\n";
      $nohit=1;
      for my $gene(keys %sequence){if ($gene=~ /$first\_$start\_$end\_(-?1)$/){$frame=$1;}}
      $string=$reference{$first};

      if (defined $positions{$start}{$comparison}){
         my $position=$positions{$start}{$comparison};
         for my $gene(keys %sequence){
            if ($gene=~ /(\S+_\S+)_(\d+)_(\d+)_(-1|1)$/){
               my ($genestr,$geneend)=($2,$3);
               if ($second eq $1){
                  $nohit=0;
                  if ($position>=$genestr && $position<=$geneend){$string=$sequence{$gene};}  
               }
            }
         }
         if ($nohit){
            $string=$reference{$first};
            foreach my $position(sort {$a<=>$b}keys %alternative){
               if (defined $alternative{$position}{$comparison}){
                  if ($position>=$start && $position<=$end){
                     substr($string,$position-1,1,$alternative{$position}{$comparison});
                  }
               }
            }
            $fragment=substr($string,$start-1,$end-$start+1);
#            $string=$fragment;
            if ($frame==-1){
               $string=reverse $fragment;
               $string=~ tr/ATGCatgc/TACGTACG/;
#               print "$frame\n$fragment\n$string\n";
            }
            else {$string=$fragment;}
         }
      }
      else{$string=$sequence{"$first\_$start\_$end\_$frame"};}
      my $genestring=substr($string,0,-3);
      print OUT ">$second\n$genestring\n";
   }
   close OUT;
   $pm->finish;
} 
$pm->wait_all_children;

open (IND,">$index_file")||die "$!";
foreach my $index(@indexArray){print IND $index,"\n";}
close IND;

