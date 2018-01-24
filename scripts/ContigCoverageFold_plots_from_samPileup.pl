#! /usr/bin/perl

use strict;
use File::Basename;
$| = 1;

my $file1=$ARGV[0]; # reference or contig file

my $file2=$ARGV[1];  # samtools pileup -f ref.fa aln.bam > rawraw.pileup  

my $outDir=$ARGV[2];

my $cov_cut_off= $ARGV[3];
                
my ($file_name, $path, $suffix)=fileparse("$file2", qr/\.[^.]*/);

unless ($file1 && $file2 && $outDir){              
   print "Usage:\nperl $0 <contigs.fasta> <samtools_pileup_file> <output directory> <coverage_cut_off>\n\n" ;
#   print "      Will output contig coverage >=80%(default) fasta file.\n";
   exit;
}
my $len_cut_off=1; # contig length cutoff
#$cov_cut_off=0.0000001 if (!$cov_cut_off or $cov_cut_off == 0);  # contig coverage cutoff
my ($id, $seq, %seq_hash, $GC_num, $GC_content,$seq_len,$total_base);
open (IN,"$file1") or die "Can't open $file1:$!";
while (<IN>)
{
    chomp;
    if(/^>(\S+)/)
    {
       if ($seq){
            $seq_len= length($seq);
            $GC_num = $seq=~ tr/GCgc/GCgc/; 
            $GC_content = $GC_num/$seq_len;
            $seq_hash{$id}->{len}= $seq_len;
            $seq_hash{$id}->{GC}=$GC_content;
            $seq_hash{$id}->{seq}=$seq if ($cov_cut_off>0);
            #for my $pos (1..$seq_len){
            #  $seq_hash{$id}->{$pos}=0;
            #}
            #print "$id\t$seq_len\n"
       }
       $id =$1;
       $seq ="";
    }
    else
    {
       $seq .= $_;
    }
}
      if ($seq){
            $seq_len= length($seq);
            $GC_num = $seq=~ tr/GCgc/GCgc/; 
            $GC_content = $GC_num/$seq_len;
            $seq_hash{$id}->{len}= $seq_len;
            $seq_hash{$id}->{GC}=$GC_content;
            $seq_hash{$id}->{seq}=$seq if ($cov_cut_off>0);
            #$total_base += $seq_len if ($seq_len>$len_cut_off);
            #for my $pos (1..$seq_len){
            #  $seq_hash{$id}->{$pos}=0;
            #}
      }


close IN;

my (@array, $total_fold,$total_covered_base,$pileup_id,$len);
my @fold_array;
my $pileup_old_id="";
my ($within_contig_sd,$covered_base,$each_id_avg_fold,$fold)=(0,0,0,0);
open (IN2,"$file2") or die "Can't open $file2:$!";
while(<IN2>)
{
   chomp; 
   @array=split /\t/,$_;
   $pileup_id = $array[0];
   next if ($array[2] eq "N" );
   if ($pileup_old_id ne $pileup_id) 
   { 
        $len=$seq_hash{$pileup_old_id}->{len};
      if ($pileup_old_id){
        ($within_contig_sd,$covered_base,$each_id_avg_fold,$fold)= &within_contig_standard_deviation($len,@fold_array);
        $seq_hash{$pileup_old_id}->{within_contig_sd}=$within_contig_sd;
        $seq_hash{$pileup_old_id}->{each_id_avg_fold}=$each_id_avg_fold;
        $seq_hash{$pileup_old_id}->{total_fold}=$fold;
        $seq_hash{$pileup_old_id}->{covered_base}=$covered_base;
      }
        undef @fold_array;
   }
 #  $seq_hash{$id}->{fold} += $array[3];
   push @fold_array, $array[3];
   #$seq_hash{$array[0]}->{$array[1]}=$array[3];
 #  $seq_hash{$id}->{covered_base}++;
   $pileup_old_id=$pileup_id;
  # $total_fold += $array[3] if ($seq_hash{$array[0]}->{len} > $cut_off);
  # $total_covered_base++ if ($seq_hash{$array[0]}->{len} > $len_cut_off);
}
    $len=$seq_hash{$pileup_old_id}->{len};
    ($within_contig_sd,$covered_base,$each_id_avg_fold,$fold)= &within_contig_standard_deviation($len,@fold_array);
    $seq_hash{$pileup_old_id}->{within_contig_sd}=$within_contig_sd;
    $seq_hash{$pileup_old_id}->{each_id_avg_fold}=$each_id_avg_fold;
    $seq_hash{$pileup_old_id}->{total_fold}=$fold;
    $seq_hash{$pileup_old_id}->{covered_base}=$covered_base;
close IN2;

open (OUT,">$outDir/${file_name}_coverage.table");
#open (OUT2,">${file_name}_base_coverage.txt");
#open (OUT3,">${file_name}_contig_coord.txt");
print OUT "ID\tLength\tGC%\tAvg_fold\tFold_std\tBase_Coverage%\n";
open (FASTA, ">$outDir/Final_contigs.fasta") if ($cov_cut_off > 0);
my $end;
foreach my $id (sort {$seq_hash{$b}->{len}<=>$seq_hash{$a}->{len}}  keys %seq_hash)
{
     my $len=$seq_hash{$id}->{len};
     if ($seq_hash{$id}->{total_fold})
     {
        $each_id_avg_fold=$seq_hash{$id}->{each_id_avg_fold};
        $within_contig_sd=$seq_hash{$id}->{within_contig_sd};
        $fold=$seq_hash{$id}->{total_fold};
        $covered_base= $seq_hash{$id}->{covered_base};
     }
     else
     {
        ($within_contig_sd,$covered_base,$each_id_avg_fold,$fold)=(0,0,0,0);
     }
     next if ($len < $len_cut_off);
     next if (!$id);
 #    print $id,"\t",$len,"\n";
     my $cov= $covered_base/$len*100;
   
     #next if ($cov < $cov_cut_off);
    
     if ($cov != 0)
     {
         
     }
     #print $id,"\t",$within_contig_sd,"\n";
    # my $fold=$seq_hash{$id}->{fold};
     $total_fold += $fold;
     $total_covered_base += $covered_base;
     $total_base += $len;
     if ($cov_cut_off>0 and $cov >= $cov_cut_off){ 
       my $seq = $seq_hash{$id}->{seq};
       $seq =~ s/(.{80})/$1\n/g; 
       chomp $seq;
       print FASTA ">$id\n$seq\n";
     }
     
     printf OUT ("%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
           $id,
           $len,
           $seq_hash{$id}->{GC} * 100,
           $each_id_avg_fold,
           $within_contig_sd,
           $cov) if ($cov >= $cov_cut_off);
           
     if ($seq_hash{$id}->{len} > $len_cut_off){      
       #print OUT3 1+$end,"\t", $end+$seq_hash{$id}->{len},"\n";
       #for my $pos (1..$seq_hash{$id}->{len}){
       #   print OUT2 $seq_hash{$id}->{$pos},"\n";
       #}
       #$end += $seq_hash{$id}->{len};
     }
     
}
close FASTA  if ($cov_cut_off> 0);
close OUT;
#close OUT2;
#close OUT3;
my $avg_cov_fold= $total_fold/$total_base;
my $coverage=$total_covered_base/$total_base * 100;
printf "Avg_coverage_fold:\t%.4f\n", $avg_cov_fold;
printf "Coverage:\t%.4f%%\n", $coverage;


open (Rscript, ">$path/Rscript$$");

  print Rscript "
#jpeg(filename=\"${file_name}_base_coverage.jpg\",width=1024,height=640,quality=100)
# init device
#par(mfrow=c(1,1))
# setup plotting area
#par(mar = c(5, 5, 5, 5), xpd=TRUE, cex.main=1.2)
#a<-read.table(file=\"${file_name}_base_coverage.txt\")
#data.contigs <- read.table(file=\"${file_name}_contig_coord.txt\")
#plot(1:length(a\$V1),a[,1],type=\"l\",col=\"blue\",cex=2,xlab=\"Contig\",ylab=\"Coverage fold (x)\",main=\"\",log=\"y\")
#pa<-par(\'usr\');
#for(i in 1:dim(data.contigs)[1]){
#  adj<-30
#  direction <- 3
#  arrows(data.contigs[i,1],round(pa[3])-adj,data.contigs[i,2],round(pa[3])-adj,code=direction, col=\"black\",lwd=2,length=0.05)
#}

#leg.txt<-paste(\"Coverage: \",format($coverage,digit=4))
#legend(\"topright\",leg.txt)
#tmp<-dev.off()

#jpeg(filename=\"${file_name}_avg_coverage_histogram.jpg\",quality=100)
#a<-read.table(file=\"${file_name}_avg_coverage.table\",header=TRUE)
#b<-round (6*sd(a\$Avg_fold))
#c<-a\$Avg_fold[a\$Avg_fold<b]
#d<-length(a\$Avg_fold[a\$Avg_fold>=b])
#h<-hist(c,breaks=b,plot=FALSE)
#new_a<-c(h\$count,d)
#plot(new_a,type=\'h\',lwd=3, col=\'black\',main=\"\",xlab=\"Coverage(fold)\",ylab=\'Frequency\',log=\"y\",xaxt=\"n\")
#x<-seq(0,round (5*sd(a\$Avg_fold)),round (5*sd(a\$Avg_fold))/5)
#x.text <- format(x,digit=1)
#axis(1,labels=x.text,at=x+1,tick=TRUE)
#axis(1,labels=paste(\">\",b),at=b+1,tick=TRUE)
#leg.txt<-paste(\"Average fold: \",format($avg_cov_fold,digit=4));
#legend(\"topright\",leg.txt);
#tmp<-dev.off();

pdf(file=\"$path/${file_name}_plots.pdf\",width=10,height=8);
par(mar=c(5,6,4,2))
a<-read.table(file=\"$outDir/${file_name}_coverage.table\",header=TRUE)

# Contigs average fold coverage vs. Contigs Length
#ylog<-\"\"
#if (max(a\$Avg_fold)>500) { ylog<-\"y\" }
#Avg_cov_depth_sd<-sd(a\$Avg_fold)
#png(filename=\"$outDir/${file_name}_avg_fold_vs_len.png\",height=800,width=800,type=\"cairo\")
plot(a\$Length/1000,a\$Avg_fold,pch=8,xlab=\"Contig Length (kbp)\",ylab=\"Average Coverage Fold (x)\", main=\"Contig Average Fold Coverage vs. Contig Length\",log=\"y\")
#DepthCovMean<-sprintf(\"Mean %.2f (x)\",mean(a\$Avg_fold))
DepthCovMean<-sprintf(\"Mean %.2f (x)\",$avg_cov_fold)
DepthCovSD<-sprintf(\"Std %.2f\",sd(a\$Avg_fold))
leg.txt<-c(DepthCovMean,DepthCovSD)
legend(\"topright\",leg.txt,inset=0.02);
#tmp<-dev.off();

# Contigs Linear Coverage vs. Contigs Length
#png(filename=\"$outDir/${file_name}_contig_coverage_vs_len.png\",height=800,width=800,type=\"cairo\")
plot(a\$Length,a\$Base_Coverage,pch=8,xlab=\"Contig Length (bp)\",ylab=\"Contig Coverage (%)\",main=\"Contig Coverage vs. Contig Length\")
OverallLinearCov <- sprintf(\"Overall Cov %.2f %%\" ,sum(a\$Length*a\$Base_Coverage/100)/sum(a\$Length)*100)
legend(\"bottomright\",OverallLinearCov)
#tmp<-dev.off();

# Contigs average fold coverage vs. Contig %GC
#png(filename=\"$outDir/${file_name}_gc_vs_avg_fold.png\",height=1024,width=1024,type=\"cairo\");
#scale<-a\$Length/sd(a\$Length)
#scale[scale<0.2]=0.2
scale<-log10(a\$Length) - 1.5
plot(a\$Avg_fold,a\$GC,ylab=\"GC (%)\",xlab=\"Average Coverage fold (x)\", main=\"Contig Average Fold Coverage vs. Contig %GC\",pch=20,col=\"blue\",cex=scale,log=\"x\");
grid(col=\"grey\")

tmp<-dev.off();

quit();
";

  if (system ("R --vanilla --slave --silent < $path/Rscript$$ 2>/dev/null")) {warn "$!\n"};
  unlink "$path/Rscript$$";
  unlink "Rplots.pdf" if ( -e "Rplots.pdf");


sub within_contig_standard_deviation {
my($len,@numbers) = @_;
#Prevent division by 0 error in case you get junk data
return undef unless(scalar(@numbers));

my $covered_base = scalar(@numbers);
if (scalar(@numbers) != $len) 
{
   for (1..($len-$covered_base)){push @numbers, 0;}
}

# Step 1, find the mean of the numbers
my $total1 = 0;
foreach my $num (@numbers) {
$total1 += $num;
}
my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
my $total2 = 0;
foreach my $num (@numbers) {
$total2 += ($mean1-$num)**2;
}
my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
my $std_dev = sqrt($mean2);
return ($std_dev,$covered_base,$mean1,$total1);
}


