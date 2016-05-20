#! /usr/bin/perl -w
# chienchi@lanl.gov
# generate some stats for Newbler or Velvet (>0.7.6) contigs. 
# Assume contigs fasta file is in the assembly output folder 
#   which contains other files, like Log (velvet), stats.txt(velvet), 
#   and 454ContigGraph.txt(Newbler) 
# 20100423
# no flag -t can generate some statistics with contig file only.

use strict;
use File::Basename;
use Getopt::Long;
my $plot;
my $assembly_tool="";
my $min_threshold=100;

GetOptions( "plot"     => \$plot,
            "tool=s"   => \$assembly_tool,
            "min=i"    => \$min_threshold,
            "help|?"   => sub {Usage()});
sub Usage
{
    print STDERR "perl $0 [options] <contigs.fasta/q> ...\n";
    print STDERR "     Can input several files\n";
    print STDERR "     -tool       The assembly tool: [Newbler] or [Velvet] for extra stats info\n";
    print STDERR "     -min        <int>  Lower threshold for contig length [default: 100]\n";
    print STDERR "     -plot       Produce statistical plots. (need R installed)\n";
    print STDERR "     -help       Print usage\n";
    exit;
}
if (scalar(@ARGV) < 1){&Usage;}
#if ($assembly_tool !~ /Newbler/i && $assembly_tool !~ /Velvet/i){&Usage;}
my %statH;
foreach my $file_index (0..$#ARGV)
{
    my $file=$ARGV[$file_index];
    my $isFasta=&isFasta($file);
    my $isFastq=&isFastq($file);
    if ( -z $file or ! -e $file) {print "The input file not exists.\n"; exit;}
    my ($file_name, $path, $suffix)=fileparse($file, qr/\.[^.]*/);
    my ($len,$total)=(0,0);
    my @x;
    my $seq_num;
    my $seq;
    my $reads_num;
    my $kmer_cov;
    my $GC_num;
    my $GC_content;
    my $velvet_id;
    my $Newbler_id;
    my $id_to_reads;
    my $id_to_cov;
    my $used_percent;
    my $singleton;
    my $exp_cov;
    my $junk;
    my $fastq;
    my $for_next_seq_check_id;
    my $seq_name;
    my ($over100k_bases,$over50k_bases,$over25k_bases,$over10k_bases,$over5k_bases,$over3k_bases,$over2k_bases,$over1k_bases)=(0,0,0,0,0,0,0,0,0);
    my ($over100k_reads,$over50k_reads,$over25k_reads,$over10k_reads,$over5k_reads,$over3k_reads,$over2k_reads,$over1k_reads)=(0,0,0,0,0,0,0,0,0);
    
    ($id_to_reads,$used_percent,$singleton,$exp_cov)=&read_velvet_stats_and_log($path) if ($assembly_tool =~ /Velvet/i);
    if ($assembly_tool =~ /Newbler/i)
    {
        my $assembled_reads=`grep -c "Assembled" $path/454ReadStatus.txt`;
        my $total_reads=`wc -l $path/454ReadStatus.txt | awk '{print \$1 - 1}'`;
        $singleton=`grep -c "Singleton" $path/454ReadStatus.txt`;
        chomp $singleton;
        $used_percent= $assembled_reads/$total_reads*100;
        if ($plot) {$id_to_cov=&read_Newbler_ContigGraph($path);}
    }
    if ($plot){
        open (PLOT,">${file}_len_gc_cov.txt");
        print PLOT "Len\tGC\tNumReads\tCov\n";
    } 
    $/= ">" if $isFasta;
    open (my $fh, $file) or die "$! $file\n";
    while(<$fh>){
            chomp;
            if ($isFastq)   # fastq format
            {  
                 $fastq=1;
                 $seq_name=$_;
                 $seq=<$fh>;
                 $seq =~ s/\n//g;
                 while ($seq !~ /\+/)
                 {
                    $seq .= <$fh>;
                    $seq =~ s/\n//g;
                 }
                 my $q_id_pos=index($seq,"+");
                 my $q_id = substr($seq,$q_id_pos);
                 $seq = substr($seq, 0, $q_id_pos);
                 $len = length $seq;
                 my $qual_seq=<$fh>;
                 $qual_seq =~ s/\n//g;
                 my $qual_seq_len = length $qual_seq;
                 while ( $qual_seq_len < $len )
                 {
                     last if ( $qual_seq_len == $len);
                     $qual_seq .= <$fh>;
                     $qual_seq =~ s/\n//g;
                     $qual_seq_len = length $qual_seq;
                 }
                 #print $seq,"\n";
                 next if ($len < $min_threshold);
                 $seq_num++;
	         $total+=$len;
                 $over100k_bases+=$len if ($len>100000);
                 $over50k_bases+=$len if ($len>50000);
                 $over25k_bases+=$len if ($len>25000);
                 $over10k_bases+=$len if ($len>10000);
                 $over5k_bases+=$len if ($len>5000);
                 $over3k_bases+=$len if ($len>3000);
                 $over2k_bases+=$len if ($len>2000);
                 $over1k_bases+=$len if ($len>1000);

                 if ($assembly_tool =~ /velvet|newbler/i){
                      $over100k_reads+=$reads_num if ($len>100000);
                      $over50k_reads+=$reads_num if ($len>50000);
                      $over25k_reads+=$reads_num if ($len>25000);
                      $over10k_reads+=$reads_num if ($len>10000);
                      $over5k_reads+=$reads_num if ($len>5000);
                      $over3k_reads+=$reads_num if ($len>3000);
                      $over2k_reads+=$reads_num if ($len>2000);
                      $over1k_reads+=$reads_num if ($len>1000);
                  }

    	         push @x,$len;
    	         $GC_num = $seq =~ tr/GCgc/GCgc/;
                 $GC_content = $GC_num/$len;
                 printf PLOT ("%d\t%.4f\t%s\t%s\n",$len,$GC_content,"NA","NA") if ($plot && $assembly_tool !~ /Velvet|Newbler/i);
                 last if (eof);
            }
            if ($isFasta)
            {
                    $_ =~ s/\>//g;
                    my ($id, @seq) = split /\n/, $_;
                    next if (!$id);
                    $seq = join "", @seq;
                    $len = length $seq;
                    if ($assembly_tool =~ /Velvet/i){
                        ($velvet_id,$kmer_cov) = $id =~ /NODE_(\d+)_length_\d+_cov_(.*)/;
                        $reads_num = $id_to_reads->{$velvet_id}; 
                    }
                    if ($assembly_tool =~ /Newbler/i){
                        ($reads_num) = $id =~ /numreads=(\d+)/;
                        ($Newbler_id) = $id =~ /^>(\S+)/;          
                    } 
       
    	            if($len>=$min_threshold)
                    {
                        $seq_num++;
	                $total+=$len;
                        $over100k_bases+=$len if ($len>100000);
                        $over50k_bases+=$len if ($len>50000);
                        $over25k_bases+=$len if ($len>25000);
                        $over10k_bases+=$len if ($len>10000);
                        $over5k_bases+=$len if ($len>5000);
                        $over3k_bases+=$len if ($len>3000);
                        $over2k_bases+=$len if ($len>2000);
                        $over1k_bases+=$len if ($len>1000);
                    
                        if ($assembly_tool =~ /velvet|newbler/i){
                             $over100k_reads+=$reads_num if ($len>100000);
                             $over50k_reads+=$reads_num if ($len>50000);
                             $over25k_reads+=$reads_num if ($len>25000);
                             $over10k_reads+=$reads_num if ($len>10000);
                             $over5k_reads+=$reads_num if ($len>5000);
                             $over3k_reads+=$reads_num if ($len>3000);
                             $over2k_reads+=$reads_num if ($len>2000);
                             $over1k_reads+=$reads_num if ($len>1000);
                         }

    			push @x,$len;
    			$GC_num = $seq =~ tr/GCgc/GCgc/;
                        $GC_content = $GC_num/$len;
                        printf PLOT ("%d\t%.4f\t%s\t%s\n",$len,$GC_content,"NA","NA") if ($plot && $assembly_tool !~ /Velvet|Newbler/i);
                        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$kmer_cov) if ($plot && $assembly_tool =~ /Velvet/i);
                        printf PLOT ("%d\t%.4f\t%d\t%.2f\n",$len,$GC_content,$reads_num,$id_to_cov->{$Newbler_id}) if ($plot && $assembly_tool =~ /Newbler/i);
                    }
    	      }
          
    }# end while
    # last fasta seq
    $/="\n" if ($isFasta);
    close PLOT if ($plot);
    close $fh;
    @x=sort{$b<=>$a} @x; 
    my $N50=0;
    my $N90=0;
    my ($count,$half)=(0,0);
    my ($top10, $top20, $top40, $top100);
    for (my $j=0;$j<@x;$j++){
            $top10+= $x[$j] if ($j<=9);
            $top20+= $x[$j] if ($j<=19);
            $top40+= $x[$j] if ($j<=39);
            $top100+= $x[$j] if ($j<=99);
    	$count+=$x[$j];
    	if (($count>=$total/2)&&($half==0)){
            $N50=$x[$j];
            $half=$x[$j]
    	}elsif ($count>=$total*0.9){
            $N90=$x[$j];
    	}
    }
    my $median=0;
    if ($seq_num % 2) {
     $median = $x[int($seq_num/2)];
    } else {
     $median = ($x[$seq_num/2] + $x[$seq_num/2 - 1]) / 2;
    }
    
    my $mean=0;
    $mean = int ($total/$seq_num);
     
    $statH{$file_index}{expCov}=$exp_cov;
    $statH{$file_index}{usedPec}=$used_percent;
    $statH{$file_index}{singleton}=$singleton;
    $statH{$file_index}{seq_num}=$seq_num;
    $statH{$file_index}{N50}=$N50;
    $statH{$file_index}{N90}=$N90;
    $statH{$file_index}{MAX}=$x[0];
    $statH{$file_index}{MIN}=$x[-1];
    $statH{$file_index}{mean}=$mean;
    $statH{$file_index}{median}=$median;
    $statH{$file_index}{total}=$total;
    $statH{$file_index}{top10}=$top10;
    $statH{$file_index}{top20}=$top20;
    $statH{$file_index}{top40}=$top40;
    $statH{$file_index}{top100}=$top100;
    $statH{$file_index}{over100kb}=$over100k_bases;
    $statH{$file_index}{over50kb}=$over50k_bases;
    $statH{$file_index}{over25kb}=$over25k_bases;
    $statH{$file_index}{over10kb}=$over10k_bases;
    $statH{$file_index}{over5kb}=$over5k_bases;
    $statH{$file_index}{over3kb}=$over3k_bases;
    $statH{$file_index}{over2kb}=$over2k_bases;
    $statH{$file_index}{over1kb}=$over1k_bases;
    $statH{$file_index}{over100kR}=$over100k_reads;
    $statH{$file_index}{over50kR}=$over50k_reads;
    $statH{$file_index}{over25kR}=$over25k_reads;
    $statH{$file_index}{over10kR}=$over10k_reads;
    $statH{$file_index}{over5kR}=$over5k_reads;
    $statH{$file_index}{over3kR}=$over3k_reads;
    $statH{$file_index}{over2kR}=$over2k_reads;
    $statH{$file_index}{over1kR}=$over1k_reads;
    
} # end foreach @ARGV
    
&plot(\%statH,@ARGV) if ($plot);

# print header
    print "File\t";
    print "Expected_coverage\t" if ($assembly_tool =~ /Velvet/i);
    if ($assembly_tool =~ /velvet|newbler/i){
      print "Assembled_reads\t";
      print "Singleton:\t";
    }
    print "Contigs_number\t";
    print "N50\t";
    print "N90\t";
    print "Max\t";
    print "Min\t";
    print "Mean\t";
    print "Median\t";
    print "Total_bases\t";
    print "Top10_bases\t";
    print "Top20_bases\t";
    print "Top40_bases\t";
    print "Top100_bases\t";
    print ">100kb_bases\t";
    print ">50kb_bases\t";
    print ">25kb_bases\t";
    print ">10kb_bases\t";
    print ">5kb_bases\t";
    print ">3kb_bases\t";
    print ">2kb_bases\t";
    print ">1kb_bases";
    if ($assembly_tool =~ /velvet|newbler/i){
      print "\t>100kb_reads\t";
      print ">50kb_reads\t";
      print ">25kb_reads\t";
      print ">10kb_reads\t";
      print ">5kb_reads\t";
      print ">3kb_reads\t";
      print ">2kb_reads\t";
      print ">1kb_reads\t";
    }
    print "\n";


## print out
foreach my $file_index (sort keys %statH)
{
    my $file = $ARGV[$file_index];
    print $file,"\t";
    print $statH{$file_index}{expCov},"\t" if ($assembly_tool =~ /Velvet/i);
    print $statH{$file_index}{usedPec},"\t",$statH{$file_index}{singleton},"\t" if ($assembly_tool =~ /velvet|newbler/i);
    print $statH{$file_index}{seq_num},"\t",
    $statH{$file_index}{N50},"\t",
    $statH{$file_index}{N90},"\t",
    $statH{$file_index}{MAX},"\t",
    $statH{$file_index}{MIN},"\t",
    $statH{$file_index}{mean},"\t",
    $statH{$file_index}{median},"\t",
    $statH{$file_index}{total},"\t",
    $statH{$file_index}{top10},"\t",
    $statH{$file_index}{top20},"\t",
    $statH{$file_index}{top40},"\t",
    $statH{$file_index}{top100},"\t",
    $statH{$file_index}{over100kb},"\t",
    $statH{$file_index}{over50kb},"\t",
    $statH{$file_index}{over25kb},"\t",
    $statH{$file_index}{over10kb},"\t",
    $statH{$file_index}{over5kb},"\t",
    $statH{$file_index}{over3kb},"\t",
    $statH{$file_index}{over2kb},"\t",
    $statH{$file_index}{over1kb};
    if ($assembly_tool =~ /velvet|newbler/i){
    print
    $statH{$file_index}{over100kR},"\t",
    $statH{$file_index}{over50kR},"\t",
    $statH{$file_index}{over25kR},"\t",
    $statH{$file_index}{over10kR},"\t",
    $statH{$file_index}{over5kR},"\t",
    $statH{$file_index}{over3kR},"\t",
    $statH{$file_index}{over2kR},"\t",
    $statH{$file_index}{over1kR};
    }
    print "\n";
    unlink "${file}_len_gc_cov.txt";
}



sub plot {
  my $statRef=shift;
  my @files=@_;
  my $numOfFile=scalar(@files);
  my %stat=%{$statRef};
  my ($maxTotalLen,$maxTotalLenIndex)=0;
  my ($maxSeqNum,$maxSeqNumIndex)=0;
  open (Rscript, ">Rscript$$");
  print Rscript 
"pdf(file=\"contigs_stats.pdf\",width=10,height=8);
par(mar=c(5,6,4,2))
LenCumSumList<-list()
length(LenCumSumList)<-$numOfFile
NxVectorList<-list()
length(NxVectorList)<-$numOfFile
myColor<-rainbow($numOfFile)
";

if ($numOfFile >1)
{
  print Rscript "
plot(1:$numOfFile*30,xaxt=\'n\',yaxt=\'n\',type=\'n\',ylab=\'\',xlab=\'\',main=\"Files\")
fileListtxt<-NULL
";
  foreach my $file_index (0..$#files)
  {
      my $file = $files[$file_index];
      my $seqNum = $stat{$file_index}{seq_num};
      my $N50 = $stat{$file_index}{N50};
      my $Max = $stat{$file_index}{MAX};
      my $Min = $stat{$file_index}{MIN};
      my $mean = $stat{$file_index}{mean};
      my $median = $stat{$file_index}{median};
      my $totalLen = $stat{$file_index}{total};
      print Rscript
      "filePrint<-paste($file_index+1,\": \", \"$file\")
stats<-paste(paste(prettyNum($seqNum,big.mark=\",\"), \"Contigs;\"),\"N50 $N50;\",paste(\"Max\",prettyNum($Max,big.mark=\",\"),\";\"),\"Min $Min;\",\"Median $median;\",\"Mean $mean;\",paste(\"Total Base\",prettyNum($totalLen,big.mark=\",\")))
fileListtxt<-c(fileListtxt,filePrint,stats,\"\")
";
   }
   print Rscript
"
for (i in 1:length(fileListtxt))
{
  text(1,($numOfFile*30-i*5),fileListtxt[i],adj=0,font=2)
}
";
}

foreach my $file_index (0..$#files)
{
  my $file = $files[$file_index];
  my $filename = basename($file);
  my $seqNum = $stat{$file_index}{seq_num};
  my $N50 = $stat{$file_index}{N50};
  my $Max = $stat{$file_index}{MAX};
  my $Min = $stat{$file_index}{MIN};
  my $mean = $stat{$file_index}{mean};
  my $median = $stat{$file_index}{median};
  my $totalLen = $stat{$file_index}{total};
  if ($seqNum > $maxSeqNum) {$maxSeqNum = $seqNum; $maxSeqNumIndex=$file_index;} 
  if ($totalLen > $maxTotalLen){$maxTotalLen = $totalLen; $maxTotalLenIndex=$file_index;} 
  print Rscript "
file$file_index<-read.table(file=\"${file}_len_gc_cov.txt\",header=TRUE);
# length histogram plot
h$file_index<-hist(file$file_index\$Len,col=\'blue\',breaks=max(file$file_index\$Len)/100, main=\"Contig Length Distribution\",xlab=\"Length (bp)\",sub=\"$filename\");
leg.txt<-c(paste(prettyNum($seqNum,big.mark=\",\"), \"Contigs\"),\"N50 $N50\",paste(\"Max\",prettyNum($Max,big.mark=\",\")),\"Min $Min\",\"Median $median\",paste(\"Total Base\",prettyNum($totalLen,big.mark=\",\")))
legend(\"topright\",legend=leg.txt,inset=0.02);

# GC histogram plot
GC<-file$file_index\$GC*100;
hist(GC,breaks=c(0:100),xlab=\"GC (%)\",ylab=\"# of contigs\",main=\"Contig GC Histogram\",sub=\"$filename\");
GCmean<-sprintf(\"Mean %.2f %%\", mean(GC))
GCsd<-sprintf(\"Std %.2f\",sd(GC))
legend(\"topright\",legend=c(GCmean,GCsd),inset=0.02)
";

  if ($assembly_tool =~ /velvet|newbler/i){
        my $kmer =  ($assembly_tool =~ /velvet/i)?"(kmer)":"";
        print Rscript
"plot(GC,file$file_index\$NumReads,xlab=\"GC (%)\",ylab=\"# of reads in contig\", main=\"Num of Reads In Contigs vs. %GC\",pch=4,log='y');

plot(file$file_index\$Len,file$file_index\$NumReads,xlab=\"Contig Length\",ylab=\"# of reads in contig\", main=\"Num of Reads In Contigs vs. Contig Length\",pch=4,log=\"y\")

plot(file$file_index\$Cov,file$file_index\$GC,ylab=\"GC (%)\",xlab=\"Mean Depth Coverage $kmer\", main=\"Mean Depth Coverages $kmer vs. %GC\",pch=20,cex=log(file$file_index\$Len/mean(file$file_index\$Len)),log=\"x\")
";
  }
print Rscript
" 
LenSort$file_index<-(sort(file$file_index\$Len,decreasing=TRUE))
LenCumSum$file_index<-cumsum(LenSort$file_index)
LenCumSumList[[$file_index+1]]<-LenCumSum$file_index
NxVector$file_index<- rep(0,101)
NxCutValue$file_index<-sum(file$file_index\$Len)*seq(0,1,0.01)
for (i in 0:100)
{   
   NxVector$file_index\[i\] <- LenSort$file_index\[LenCumSum$file_index >= NxCutValue$file_index\[i\]\][1]
}
NxVectorList[[$file_index+1]]<-NxVector$file_index

";

}

if ($numOfFile == 1) 
{
  print Rscript 
"
# Cumulative Length plot
plot(LenCumSum0/1000000,xlab=\"Number of Contig\",ylab=\"Size (Mbp)\",type=\"l\",main=\"Contig Cumulative Length\",lwd=2)
mtext(side=3,\"Contigs are ordered from largest to smallest.\",adj=0,cex=0.7)
grid(col=\"grey\")

# Nx plot
plot(0:100,NxVector0/1000,type=\"l\",xlab=\"Nx\",ylab=\"Size (Kbp)\", main=\"Contig Nx Length\",lwd=2)
grid(col=\"grey\")
arrows(50+2,(NxVector0/1000)[51]+4,50,(NxVector0/1000)[51],length=0.1,angle=20)
text(50+4,(NxVector0/1000)[51]+4,paste(\"N50:\",(NxVector0/1000)[51]),pos=3,adj=1,cex=0.8)
";
}
if ($numOfFile >1)
{
  print Rscript
# Length Cumulative length
"
# Comparison plots
# Cumulative Length plot
plot(LenCumSum$maxSeqNumIndex/1000000,xlab=\"Number of Contig\",ylab=\"Size (Mbp)\",type=\"n\",main=\"Contig Cumulative Length\",ylim=c(0,$maxTotalLen/1000000))
mtext(side=3,\"Contigs are ordered from largest to smallest.\",adj=0,cex=0.7)
grid(col=\"grey\")
#legend(\"bottomright\",leg.txt,inset=0.02)
for (i in 1:$numOfFile)
{
    lines(1:length(LenCumSumList[[i]]),LenCumSumList[[i]]/1000000,col=myColor[i],lwd=3)
}
legend('bottomright',legend=1:$numOfFile,lwd=3,col=myColor,bty=\"n\",inset=0.02)

# Nx plot 
plot(0:100,NxVector$maxTotalLenIndex/1000,type=\"n\",xlab=\"Nx\",ylab=\"Size (Kbp)\", main=\"Contig Nx Length\")
grid(col=\"grey\")
for (i in 1:$numOfFile)
{
   lines(0:100,NxVectorList[[i]]/1000,col=myColor[i],lwd=3)
}
abline(v=50,col = \"gray60\",lty=2)
legend('topright',legend=1:$numOfFile,lwd=3,col=myColor,bty=\"n\",inset=0.02)
#arrows(50+4,(NxVector/1000)[51]+4,50,(NxVector/1000)[51],length=0.1,angle=20)
#text(50+4,(NxVector/1000)[51]+4,paste(\"N50:\",(NxVector/1000)[51]),pos=3,adj=1,cex=0.8)
";
 }



  print Rscript "quit();";
  #system ("R --vanilla --slave --silent < Rscript$$ 2>/dev/null");
  system ("R --vanilla --slave --silent < Rscript$$");
  unlink "Rscript$$";

  #unlink "${file}_len_gc_cov.txt";
}

    sub read_velvet_stats_and_log
    {
        my $path=shift;
        my @array;
        my %id_to_reads_num;
        # read stats.txt pile which contains numReads info each contigs.
        open (IN,"$path/stats.txt") || die "cannot open velvet stats.txt file";
        while (<IN>)
        {
           chomp;
           if ($_=~/^\d+/)
           {
               @array= split /\t/, $_;
               $id_to_reads_num{$array[0]}=$array[12];
           }
        }
        close IN;
        # raed  Laafile which contains 
        open (IN,"$path/Log") || die "cannot open velvet Log file";
        my ($used_percent,$singleton,$exp_cov);
        while(<IN>)
        {
           chomp;
            #Median coverage depth = 4.343750
            if ($_=~/Median coverage depth/)
            {
               ($exp_cov)= $_=~/(\d+\.\d+)/;
            }
            if ($_=~/^Final/){
               my ($used,$total)=$_=~/using\s+(\d+)\/(\d+)/;
               $used_percent = $used/$total*100;
               $singleton = $total - $used;
            }
        }
        close IN;
        return (\%id_to_reads_num,$used_percent,$singleton,$exp_cov);
    }
    sub read_Newbler_ContigGraph
    {
        my $path=shift;
        my $isFasta=shift;
        my @array;
        my %id_to_cov;
        open (IN,"$path/454ContigGraph.txt") ||die "cannot open 454ContigGraph.txt file";
        while (<IN>)
        {
           chomp;
           @array=split /\t/,$_;
           $id_to_cov{$array[1]}=$array[3];
        }
        close IN;
        return (\%id_to_cov);
    }
sub isFastq {
    #given a file it will return 1 for fastq, 0 for others.
    my $file=shift;
    my $fastq=0;
    open (my $fh, "< $file") or die "$!\n";
    my $header=<$fh>;
    close $fh;
    $fastq=1 if ($header =~ /^@/);
    return $fastq;
}

sub isFasta {
    #given a file it will return 1 for fasta, 0 for others.
    my $file=shift;
    my $fasta=0;
    open (my $fh, "< $file") or die "$!\n";
    my $header=<$fh>;
    close $fh;
    $fasta=1 if ($header =~ /^>/);
    return $fasta;
}
