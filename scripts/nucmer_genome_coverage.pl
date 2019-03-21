#! /usr/bin/env perl
# 20130125
use strict;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";



my $identity_cutoff=85;
my $end_cutoff=0.05;
my $query_cov_cutoff=0;
my $ref_cov_cutoff=0;
my $dotplot=0;
my $mincluster=65;
my $prefix="Output";
my $delta_file;
GetOptions('e=f', \$end_cutoff,
           'i=f', \$identity_cutoff,
 #          'q=f',  \$query_cov_cutoff,
 #          'r=f',  \$ref_cov_cutoff
           'c=i',  \$mincluster,
           'p=s'  , \$prefix,
           'd'    ,  \$dotplot,
           'delta=s', \$delta_file  
          );
          
my $referenceFile=$ARGV[0];
my $queryFile=$ARGV[1];
#my ($file_name, $file_path, $file_suffix)=fileparse("$file1", qr/\.[^.]*/);

if (scalar(@ARGV)<1){              
   print  "Usage:\nperl $0 ReferenceFasta QueryFasta\n" ;
   print  "       Following options are for defining mis-assembled contigs\n";
   print  "      -e      ratio of the query length on each end is ignored (default 0.05, >=90% query should align)\n";
#   print  "      -r     at least this percent alignment coverage in the reference sequence [0..100]\n";
#   print  "      -q     at least this percent alignment coverage in the query sequence [0..100]\n";
   print  "      -c      mincluster   Sets the minimum length of a cluster of matches (default 65)\n";
   print  "      -i      at least this percent identifity between alignement [0..100] (default 85)\n";
   print  "      -d      dotplot for first 50 reference if gnuplot existed\n";
   print  "      -p      output prefix\n";
   print  "      -delta  provide delta file to skip nucmer alignment step.\n";
   exit;
}
#$identity_cutoff=$identity_cutoff*100;
system("show-coords -h 2> /dev/null") == 0
        || die "show-coords is not in your PATH; please add it and re-run";
       
my $logFile="$prefix.log";
open (my $log_fh,">$logFile") or die "$logFile $!\n"; 
if ($delta_file)
{
   system ("show-coords -clTr $delta_file > $prefix.coords");
   system ("delta-filter -1 -i $identity_cutoff $delta_file > $prefix.filterforSNP");
   system ("show-snps -CT $prefix.filterforSNP > $prefix.snps"); 
}
else
{
   system ("nucmer --maxmatch -c $mincluster --prefix=$prefix $referenceFile $queryFile 2>/dev/null");
   system ("show-coords -clTr $prefix.delta > $prefix.coords");
   system ("delta-filter -1 -i $identity_cutoff $prefix.delta > $prefix.filterforSNP");
   system ("show-snps -CT $prefix.filterforSNP > $prefix.snps");
}

my ($numSNPs,$numINDELs,$variantCount_r)=&SNP_INDEL_COUNT("$prefix.snps");

print $log_fh "Mapping criteria\n";
#$end_cutoff=0.5 if ($end_cutoff>0.5);
#printf (">= %.1f%% query length in the middle should align\n",100-$end_cutoff*100*2);
printf $log_fh ("nucmer -c the minimum length of a cluster of matches %d \n",$mincluster);
printf $log_fh ("Aligned portion should be >= %.2f %% identity\n",$identity_cutoff);

my %query_used;
my %ref_used;
my %correctAssembled;
open (COORDS,"$prefix.coords") or die "Can't open $prefix.coords:$!";
my $coord_header= <COORDS>;
$coord_header=~ s/\n//;
#my ($referenceFile,$queryFile)=split (/\s+/,$coord_header);

# get reference info
print "Loading reference $referenceFile ...\n";
my %reference=&getSeqInfo($referenceFile,0);
# get query info and store sequence in hash.
print "Loading query $queryFile...\n";
my %query=&getSeqInfo($queryFile,1);
my $total_reads_num = scalar (keys (%query));
print $log_fh  "Total_contigs:\t$total_reads_num\n";

my ($total_identity,$total_hit)=(0,0);
open (my $Rtable,">${prefix}.coordsR");
print $Rtable "hit_start\thit_end\tsimilarity\ttemplate_name\tquery_name\tstrand\n";
while (<COORDS>)
{
    if ($_=~ /^\d+/)
    {
       my @array=split;
       my $r_start=$array[0];
       my $r_end=$array[1];
       my $q_start=$array[2];
       my $q_end=$array[3];
       my $r_align_len=$array[4];
       my $idy=$array[6];
       my $ref_len=$array[7];
       my $query_len=$array[8];
       my $ref_cov=$array[9];
       my $query_cov=$array[10];
       my $ref_tag=$array[11];
       my $query_tag=$array[12];
       my $r_loopstart=$r_start;
       my $r_loopend=$r_end;
       my $q_loopstart=$q_start;
       my $q_loopend=$q_end;
       my $strand = "1";
       if ($r_start>$r_end)
       {
           $r_loopstart=$r_end;
           $r_loopend=$r_start;
       }
       if ($q_start>$q_end)
       {
           $q_loopstart=$q_end;
           $q_loopend=$q_start;
           $strand="-1" 
       }
       my $allowed_error_len = $query_len * $end_cutoff;
       #print $query_len,"\t",$allowed_error_len,"\n";
      # next if ($ref_cov < $ref_cov_cutoff); 
          if (   #define aligned
              $idy>=$identity_cutoff  &&
              $q_loopstart <= $allowed_error_len &&
              $q_loopend >= ($query_len - $allowed_error_len)
          )
          {
              for (my $i=$r_loopstart; $i<= $r_loopend; $i++)
              {
                  $reference{$ref_tag}->{$i}++; 
              }
              for (my $i=$q_loopstart; $i<= $q_loopend; $i++)
              {
                  $query{$query_tag}->{$i}++;
              }
              $ref_used{$ref_tag}->{$query_tag}=1;
              $correctAssembled{$query_tag}=1;
              $total_identity += $idy;
	      $total_hit++;
              print $Rtable "$r_loopstart\t$r_loopend\t$idy\t$ref_tag\t$query_tag\t$strand\n" if ($r_align_len>=500);
          }
          $query_used{$query_tag}=1;
    }
    
}
close COORDS;

my $total_ref_len;
my $total_ref_cov;
my $total_fold_cov;
my $avg_cov_fold;
my $coverage;
my $used_read_num;
my $mis_assembled_num;
open (OUT, ">${prefix}_avg_coverage.table");
#aopen (OUT2,">${prefix}_base_coverage.txt");
#open (OUT3,">${prefix}_contig_coord.txt");
open (OUT2,">${prefix}_ref_zero_cov_coord.txt");
#open (UNUSEDREF,">${prefix}_ref_unUsed.fasta");

print OUT "ID\tLength\tGC(%)\tMapped_contigs\tBase_Coverage(%)\tAvg_fold(x)\tNumer_of_Gap\tTotal_Gap_bases\tNum_of_SNPs\tNum_of_INDELs\n";
my $unUsedRef=0;
#foreach my $ref (sort {$reference{$b}->{len}<=>$reference{$a}->{len}} keys %reference)
open (Rplot,">${prefix}.R") || die ": $!";
print Rplot "pdf(file = \"${prefix}_plot.pdf\",width = 10, height = 8)\n
par(mar=c(5,6,4,2))
# read the data
data.hits.all  <- read.table(file=\"${prefix}.coordsR\", header=TRUE, dec=\".\")
data.gaps.all  <- read.table(file=\"${prefix}_ref_zero_cov_coord.txt\")\n";
 
my $ref_count=0;
foreach my $ref (sort { $reference{$b}->{len} <=> $reference{$a}->{len} } keys  %reference){
  my $mapped_contig_count = scalar( keys %{$ref_used{$ref}});
  $reference{$ref}->{mappedCount} = $mapped_contig_count;
  $variantCount_r->{$ref}->{INDELs} = $variantCount_r->{$ref}->{INDELs} || 0;
  $variantCount_r->{$ref}->{SNPs} = $variantCount_r->{$ref}->{SNPs} || 0;
  $ref_count++;
  if ($ref_count <= 50  && $dotplot){
        (my $output_prefix = $ref) =~  s/\W/\_/g;
        $output_prefix = ${prefix}."_".$output_prefix."_dotplot";
  	system("mummerplot --fat -png $prefix.delta -r \"$ref\" --large  --prefix $output_prefix 1>/dev/null 2>/dev/null") if (`which gnuplot 2>/dev/null`);
  	system("rm -f $output_prefix.*plot $output_prefix.gp $output_prefix.filter");
  }
}

system("mummerplot --fat -png $prefix.delta --large --prefix $prefix.dotplot 1>/dev/null 2>/dev/null") if (`which gnuplot 2>/dev/null`);
system("rm -f $prefix.*plot $prefix.gp $prefix.filter");

foreach my $ref (sort {$reference{$b}->{mappedCount} <=> $reference{$a}->{mappedCount}} keys %reference)
{
  my $ref_len=$reference{$ref}->{len};
  my $ref_desc=$reference{$ref}->{desc};
  my $mapped_contig_count =  $reference{$ref}->{mappedCount};
  my ($genome_covered_pos, $gaps, $gap_num,$cov,$total_cov) = (0,0,0,0,0);
  my @gap_array;
  if ($ref_used{$ref})
  {   
      for (1..$ref_len){
        if ($reference{$ref}->{$_})
        {
            if (@gap_array){
             $gap_num++;
             print OUT2 $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t$ref","\n";
             @gap_array=();
            }
            $genome_covered_pos++;
            $cov=$reference{$ref}->{$_};
            $total_cov += $cov;
        }
        else
        {
           push @gap_array, $_;
           $gaps++;
        }
      }

      if (@gap_array){
             $gap_num++;
             print OUT2 $gap_array[0],"\t",$gap_array[-1],"\t",$gap_array[-1] - $gap_array[0]+1,"\t$ref","\n";
      }
      
      print OUT $ref,"\t",$reference{$ref}->{len},"\t",$reference{$ref}->{GC},"\t",
	    $mapped_contig_count,"\t",
            $genome_covered_pos/$ref_len*100,"\t",
            $total_cov/$ref_len,"\t",
	    $gap_num,"\t",
	    $gaps,"\t",
	    $variantCount_r->{$ref}->{SNPs},"\t",
	    $variantCount_r->{$ref}->{INDELs},"\n";
	    
	     
      $total_ref_cov += $genome_covered_pos;
      $total_fold_cov += $total_cov;
      
  }
  else
  {
      $unUsedRef++;
      print OUT $ref,"\t",$reference{$ref}->{len},"\t",$reference{$ref}->{GC},"\t",
            0,"\t",
            0,"\t",
            0,"\t",
            1,"\t",
            $reference{$ref}->{len},"\t",
            0,"\t",
            0,"\n";

      #print UNUSEDREF ">$ref\n".$reference{$ref}->{seq},"\n";
  }
  $total_ref_len += $ref_len;
  
  my $xlab= ($ref_desc)?$ref_desc:$ref;
  print Rplot <<RCODE;

data.hits <- subset(data.hits.all,template_name=="$ref")
data.gaps <- subset(data.gaps.all,V4=="$ref")
# variable
data.dimx <- $ref_len

ymin <- min(data.hits[,3])-20
plot(0,0,col="#00000000",pch=19,cex=2, las=1, xlim=c(0,data.dimx), ylim=c(ymin,100), xlab='', ylab='', xaxt ="n")
mtext('Identity (%)', side=2, line=3, cex=1.0)
mtext("$xlab", side=1, line=2, cex=1)
	
# set the x-axis ticks and labels
axisXS.at <- seq(0,data.dimx,data.dimx/10)
axisXS.text <- format(axisXS.at,digit=2,scientific=TRUE)
axis(1, at = axisXS.at, labels=FALSE, tcl=0.2, tck=-.01)
mtext(axisXS.text, side=1,at=axisXS.at, cex=1)

# plot all the hits
for(i in 1: dim(data.hits)[1]) {
 x<-c(data.hits[i,1],data.hits[i,2])
 y<-c(data.hits[i,3],data.hits[i,3])
 strand<-data.hits[i,6]
 color <- "red"
 if (strand == "-1")
 {
   color <- "green"
 }
 #points(x[1], y[1], type="p", cex=0.3, pch=16, col="black")
 rect(x[1],y[1]-0.2,x[2],y[1],border=color,col=color)
}

# get margin coordiates
pa<-par('usr');
# plot gap regions
if (dim(data.gaps)[1] > 0){
	for(i in 1:dim(data.gaps)[1]){
		rect(data.gaps[i,1],round(pa[3]),data.gaps[i,2],round(pa[3])+0.5,col="black",border=NA)	
	}
}

# add Legend
Coverage <- sprintf ("Coverage: %.2f %%", $genome_covered_pos/$ref_len*100)
legend('bottomright',Coverage,inset=0.05)

# add Title
title("Reference-based Analysis: Contig Mapped to Reference $ref")

RCODE

  
}
# Overall fold and coverage
$avg_cov_fold=$total_fold_cov/$total_ref_len;
$coverage=$total_ref_cov/$total_ref_len * 100;
close OUT;
close OUT2;
close UNUSEDREF;
print Rplot "garbage <- dev.off();\nq()\n";
close Rplot;
system ("R --vanilla --slave --silent < ${prefix}.R 2>/dev/null"); 
unlink "$prefix.R";
unlink "${prefix}.coordsR";


# Parse novel region in query.
open (OUT3,">${prefix}_query_novel_region_coord.txt");
#open (OUT5,">${prefix}_query_novel_region_coord.fasta");
open (OUT4,">${prefix}_query_novel_region_30bpUP.fasta");
#open (OUT4,">${prefix}_query_unUsed.fasta");
open (OUT6,">${prefix}_query_misAssembled.fasta");
my $misAssemblyNum=0;
my $unUsedContig=0;
my $allowQueryEndRegionForUnMap=20;
foreach my $q (sort keys %query)
{
  my $query_len=$query{$q}->{len};
  if ($correctAssembled{$q})
  {
     my ($genome_covered_pos, $gaps, @gap_array, $gap_num,$cov,$total_cov);
     my ($gap_start,$gap_end,$gap_len);
     for (1..$query_len){
       if ($query{$q}->{$_})
       {
          if (@gap_array){
             $gap_num++;
             $gap_start=$gap_array[0];
             $gap_end=$gap_array[-1];
             $gap_len=$gap_array[-1] - $gap_array[0]+1;
             unless (($gap_start > ($query_len - $allowQueryEndRegionForUnMap)) or ($gap_end < $allowQueryEndRegionForUnMap) )
             {
               print OUT3 $gap_start,"\t",$gap_end,"\t",$gap_len,"\t$q","\n";
               my $gap_seq=substr($query{$q}->{seq},$gap_start-1,$gap_len);
               print OUT4 ">$q"."_".$gap_start."_".$gap_end."\n".$gap_seq."\n" if ($gap_len >=30);
             }
             @gap_array=();
          }
          $genome_covered_pos++;
          $cov=$query{$q}->{$_};
          $total_cov += $cov;
       }
       else
       {
          push @gap_array, $_;
          $gaps++;
       }
     }
     if (@gap_array){
            $gap_num++;
             $gap_start=$gap_array[0];
             $gap_end=$gap_array[-1];
             $gap_len=$gap_array[-1] - $gap_array[0]+1;
             unless (($gap_start > ($query_len - $allowQueryEndRegionForUnMap)) or ($gap_end < $allowQueryEndRegionForUnMap) )
             {
               print OUT3 $gap_start,"\t",$gap_end,"\t",$gap_len,"\t$q","\n";
               my $gap_seq=substr($query{$q}->{seq},$gap_start-1,$gap_len);
               print OUT4 ">$q"."_".$gap_start."_".$gap_end."\n".$gap_seq."\n" if ($gap_len >=30);
             }
     }

    # print OUT $q,"\t",$query{$q}->{len},"\t",$query{$q}->{GC},"\t",
    #     $total_cov/$query_len,"\t",
    #    $genome_covered_pos/$query_len*100,"\n";
  #$total_query_cov += $genome_covered_pos;
  #$total_query_fold_cov += $total_cov;
  }
  elsif($query_used{$q})
  {
    $misAssemblyNum++;
    print OUT6 ">$q\n".$query{$q}->{seq}."\n";   
  }
  else
  {
    $unUsedContig++;
    print OUT4 ">$q\n".$query{$q}->{seq}."\n";
     
  }
  #$total_query_len += $query_len;
}
close OUT3;
close OUT4;
#close OUT5;
close OUT6;
  my $average_idn = ($total_hit)?$total_identity/$total_hit:0; 
  #printf  $log_fh ("MisAssembled Contigs#:\t%d (%.2f %%)\n",$misAssemblyNum,$misAssemblyNum/$total_reads_num*100);
  printf  $log_fh ("Unused Contigs#:\t%d (%.2f %%)\n",$unUsedContig, $unUsedContig/$total_reads_num*100);
  printf  $log_fh ("Avg_coverage_fold:\t%.4f\n",$avg_cov_fold);
  printf  $log_fh ("Reference_Coverage:\t%.2f%%\n",$coverage);
  printf  $log_fh ("Avg_Identity:\t%.2f%%\n",$average_idn);
  printf  $log_fh ("Number of SNPs:\t%d\n",$numSNPs);
  printf  $log_fh ("Number of INDELs:\t%d\n",$numINDELs);

close $log_fh;
system ("cat $logFile");
#clean
if ( -z "${prefix}_query_unUsed.fasta") { unlink "${prefix}_query_unUsed.fasta";}
#if ( -z "${prefix}_ref_unUsed.fasta") { unlink "${prefix}_ref_unUsed.fasta";}
if ( -z "${prefix}_query_novel_region_coord.txt") { unlink "${prefix}_query_novel_region_coord.txt";}
if ( -z "${prefix}_query_misAssembled.fasta") {unlink "${prefix}_query_misAssembled.fasta"; }
unlink "$prefix.filterforSNP";

  

sub getSeqInfo 
{
    my ($reference_file,$load_seq) = @_;
    my ($id, $desc, $seq, %hash, $GC_num, $GC_content,$seq_len);
    open (REF,"$reference_file") or die "Can't open $reference_file:$!";
    while (<REF>)
    {
        chomp;
        if(/^>(\S+)\s*(.*)/)
        {
           if ($seq){
               $seq_len= length($seq);
               $GC_num = $seq=~ tr/GCgc/GCgc/; 
               $GC_content = $GC_num/$seq_len;
               $hash{$id}->{len}= $seq_len;
               $hash{$id}->{desc}= $desc;
               $hash{$id}->{GC}=sprintf("%.2f",$GC_content*100);
               $hash{$id}->{seq}=$seq if ($load_seq==1);
               for my $pos (1..$seq_len){
                 $hash{$id}->{$pos}=0;
               }
           }
           $id =$1;
           $desc=$2;
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
       $hash{$id}->{len}= $seq_len;
       $hash{$id}->{desc}= $desc;
       $hash{$id}->{GC}=sprintf("%.2f",$GC_content*100);
       $hash{$id}->{seq}=$seq if ($load_seq==1);
       for my $pos (1..$seq_len){
           $hash{$id}->{$pos}=0;
       }
    }
    close REF;
    return %hash;
}

sub SNP_INDEL_COUNT 
{
    my $snp_file=shift;
    my $SNPs_count=0;
    my $Indels_count=0;
    my %count;
    my %indel;
    my $delet_base=0;
    open (my $fh, $snp_file) or die "$! $snp_file\n";
    while (my $line=<$fh>)
    {
        chomp $line;
        next if ($line !~ /^\d+/);
        my @array=split /\s+/, $line;
        my $r_position = $array[0];
        my $r_base = $array[1];
        my $q_base = $array[2];
        my $q_position = $array[3];
        my $r_tag = $array[8];
        my $q_tag = $array[9];
        
        if ($r_base =~ /\w/ and $q_base =~ /\w/)  # SNPs
        {
            $count{$r_tag}->{SNPs}++;
            $SNPs_count++;
        }
        else  # Indels 
        {
		if ($r_base eq '.'){ #insertion
                        $indel{$r_tag}->{$r_position}->{seq} .= $q_base;
			#$count{$r_tag}->{INDELs}++;
                }
                if ($q_base eq '.'){ # deletion
                        if ($indel{$r_tag}->{$r_position-$delet_base})
                        {
                                $r_position = $r_position - $delet_base;
                                $indel{$r_tag}->{$r_position}->{seq}.= $r_base;
                                $delet_base++;
			#	$count{$r_tag}->{INDELs}++;
                        }else {
                                $indel{$r_tag}->{$r_position}->{seq}.= $r_base;
                                $delet_base=1;
                        }

                }
        }
        
    }
    close $fh;
    foreach my $id (keys %indel){
        my %pos = %{$indel{$id}};
	$count{$id}->{INDELs} = scalar (keys %pos);
        $Indels_count += scalar (keys %pos);
    }
    return ($SNPs_count, $Indels_count,\%count);
}
