#!/usr/bin/perl -w

use strict;
use FindBin qw($RealBin);
use lib "$RealBin";
use lib "$RealBin/../lib/";
use lib "$RealBin/../ext/lib/perl5";
use FileHandle;
use Getopt::Long;
use Parallel::ForkManager;

my $dir;
my $stat;
my $file;
my $list;
my $threads=1;
my $seq;
my $header;
my @strings;
my $count=0;
my ($refid,$queryid,$rpos,$qpos,$snp,$start,$end);
my %coords;
my %genes;
my $fragment;
my @header_list;
my ($first,$second);
my $gff_file;
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
   't|thread=i'   => \$threads,
   'g|gene=s'     => \$gff_file,
   'p|gap=s'      => \$gapfile,
);

my $maxthreads = ($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
if ($threads <1 || $threads>$maxthreads){die ("-thread value must be between 1 and $maxthreads.\n");}

my $genedir=$dir.'/PSgenes';
if (!-d $genedir){`mkdir $genedir`;}

my $pm=new Parallel::ForkManager($threads);
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

my $fh1=FileHandle->new($gff_file)|| die "$!";
if ($fh1->open("< $gff_file")){
	while (<$fh1>){
		chomp;
		next if (/^#/);
		my @array = split /\t/,$_;
		if ( scalar(@array) ==  9 ){
			#$id,$source,$type,$start,$end,$score,$strand,$phase,$Attributes
			my $chromosome_id = $array[0];
			my $type=$array[2];
			my $start=$array[3];
			my $end = $array[4];
			my $strand = ($array[6] eq '+')? 1 : -1;
			my $Attributes = $array[8];
			if ($type eq "CDS"){
				my %annotations=map { split /=/;} split /;/,$Attributes;
				my $gene_id =  $annotations{"ID"} ||  $annotations{"Name"};
				$gene_id =~ s/\W/_/g;
				$genes{"$start:$end"}->{id} = $gene_id;		
				$genes{"$start:$end"}->{strand} = $strand;		
			}	
		}
		
		
	}
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
my $stat_header = <IN>;
while (<IN>){
	chomp;
	my ($refid,$queryid,$type,$rpos,$qpos,$rbase,$qbase,$start,$end) = split /\t/,$_;
	if ($type eq "coding SNP"){
		if ($start==0){$start=1;}
		if (abs $start-$end> 300){
  	         #print "$start\t$end\n";
			$coords{"$start:$end"}=$start;
			$alternative{$rpos}{"$refid:$queryid"}=$qbase;
		}
	}
}
close IN;

my $temp=keys %coords;
my $size = length $temp;
my $index_file=$genedir.'/gene_index.txt';
my $gapgenes=$genedir.'/gene_gaps.txt';
my $compliment;

open (GAP,">$gapgenes")||die "$!";
OUTER:foreach my $coord (sort{ $coords{$a} <=> $coords{$b} }keys %coords ){
   my ($start,$end)= $coord =~ /(\d+):(\d+)/;

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
   my $gene_id=$genes{"$start:$end"}->{id};
   $outfile = "$genedir/${gene_id}_${start}_${end}.fna";
   #$outfile=$genedir.'/Gene';
   #$outfile.=sprintf "%0${size}d", $count;
   my $index=sprintf "%0${size}d", $count;
   #$outfile.='.fna';
#   print "$index\n";
   push (@indexArray,"$index\t$gene_id\t$start\t$end") if ($genes{"$start:$end"});
   $pm->start and next;
   open (my $fh,">$outfile")||die "$!";
#   print "$start\t$end\n";

   foreach my $comparison(@header_list){
      if ($comparison=~/(.+):(.+)/){
         ($first,$second)=($1,$2);
         print $fh ">$second\n";
#          print "$comparison\t$first\n";
      }

	
      if ($genes{"$start:$end"}){
            my $gene_len = $end - $start + 1;
            my $gene=substr($reference{$name},$start-1,$gene_len);
#         print "$entry\n";
            my $strand=$genes{"$start:$end"}->{strand};
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
            print $fh "$gene\n";
#            print "$gene\n";
      }
   }
   close $fh;
   $pm->finish;
}
$pm->wait_all_children;

open (IND,">$index_file")||die "$!";
foreach my $index(@indexArray){print IND $index,"\n";}
close IND;

