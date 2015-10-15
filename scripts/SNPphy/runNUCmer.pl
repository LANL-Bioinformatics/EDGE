#!/usr/bin/perl -w

######################################################
# Written by Sanaa Ahmed
# Nov. 30, 2012

# Given a directory containing fasta files, runs nucmer
#  Runs nucmer only on files ending in .fna 
#  Runs nucmer on all against all

######################################################

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use File::Basename;
use Parallel::ForkManager;

my $breaklen=200;
my $mincluster=65;
my $diagfactor=0.12;
my $maxgap=90;
my $minmatch=20;
my $identity=0;
my $gap_cutoff=0;
my $repeat_identity=97;
my $len_cutoff= 100;
my ($query_dir,$thread,$list,$reference);
my @query;
my $outdir=`pwd`;
   $outdir =~ s/\n//;

GetOptions(
   'q|querydir=s'   => \$query_dir,
   'd|outdir=s'     => \$outdir,
   'r|reference=s'  => \$reference,
   't|thread=i'     => \$thread,
   'l|list=s'       => \$list, 
   'h|help'         => sub{usage()}
);

&usage() unless ($query_dir);

my $max_thread = ($^O =~ /darwin/)?  `sysctl hw.ncpu | awk '{print \$2}'`:`grep -c ^processor /proc/cpuinfo`;
if ($thread < 1 || $thread>$max_thread){
   die("-thread value must be between 1 and $max_thread.\n");
}

chdir $outdir;
if (! -e "gaps"){mkdir "gaps";}
if (! -e "snps"){mkdir "snps";}
if (! -e "stats"){mkdir "stats";}
if (! -e "temp"){mkdir "temp";}
my $options= "--maxmatch -b $breaklen -c $mincluster -d $diagfactor -g $maxgap -l $minmatch ";
$gap_cutoff= "-l $gap_cutoff";
$identity= "-i $identity";

my %snp_hash;
my %indel_hash;
my %gap_hash;
my $snp_indel;
my $snp_n;
my $indel_n;
my $gaps;
my $ref_gaps;
my $query_gaps;

my $bindir= getBinDirectory();
read_directory($query_dir);
run_self_nucmer(@query);
if ($reference)
{
    #print "Test\n";
    run_nucmer_with_ref($reference,@query);
}
else
{
    run_nucmer(@query);
}
cleanup();

sub read_directory
{
my $dir=shift; 
my $query_file;
my %queries;

if ($dir=~ /.+\/$/){my $temp= chop($dir);}
$dir = $dir.'/files';

open (IN,$list);
while (<IN>){
   chomp;
   $queries{$_}++;
}
close IN;

opendir(PARENT,$dir);
while (my $files= readdir(PARENT)){  
   next if ($files=~ /^..?$/);
   if ($files=~ /(.+)\.fna$/){
 #     next if ($files=~ /contig/);
      if (exists $queries{$1}){
         $query_file=$dir.'/'.$files;
         push(@query,$query_file);
      }
   }
}
closedir(PARENT);
}

sub run_self_nucmer
{
my $ident=95;
my $len_cutoff= 75;

my $pm= new Parallel::ForkManager($thread);

foreach my $reference(@query){
   $pm->start($reference) and next;
   my ($name,$path,$suffix) = fileparse("$reference", qr/\.[^.]*/);
   my $coords=$outdir.'/'.$name.'_repeat_coords.txt';
   my $stat=$outdir.'/'.$name.'_repeat_stats.txt';
   my $fasta=$outdir.'/'.$name.'_norepeats.fna';
   
   my $command= "$bindir/get_repeat_coords.pl -l $len_cutoff -i $repeat_identity -o $coords -s $stat $reference";
   if (system ($command)){die "Error running $command.\n";}
   my $remove_repeats="$bindir/removeRepeats.pl -f $reference -c $coords -o $fasta";
   if (system ($remove_repeats)){die "Error running $remove_repeats.\n";}
   $pm->finish(0);
}
$pm->wait_all_children;
print "\nRepeat coords for all references found.\n";
print "Self-NUCmer complete\n\n";
}

sub run_nucmer
{
my $iteration= combo(2,@query);
my $pm= new Parallel::ForkManager($thread);

$pm->run_on_finish(
   sub{
      my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$hash_r)=@_;
# retrieve data structure from child
      if (defined($hash_r)){  # children are not forced to send anything
         foreach my $key1 (keys %{$hash_r}){
            foreach my $key2 (keys %{$hash_r}){
               $snp_hash{$key1}->{$key2} = $hash_r->{$key1}->{$key2}->{SNP};
               $indel_hash{$key1}->{$key2} = $hash_r->{$key1}->{$key2}->{indel};
               $gap_hash{$key1}->{$key2} = $hash_r->{$key1}->{$key2}->{gap};
            }
         }
      } 
# problems occuring during storage or retrieval will throw a warning
      else {print qq|No message received from child process $pid!\n|;}
   }
);

while (my @combo= $iteration->()){
   my ($first_name,$first_path,$first_suffix) = fileparse("$combo[0]", qr/\.[^.]*/);
   my ($second_name,$second_path,$second_suffix) = fileparse("$combo[1]", qr/\.[^.]*/);
   next if ( -e "$outdir/snps/${first_name}_${second_name}.snps");
   $pm->start(@combo) and next;
   my %hash;
   my $reference= "$first_path$first_name$first_suffix";
   my $query= "$second_path$second_name$second_suffix";
   $first_name =~ s/\.fna//;
   $second_name =~ s/\.fna//;
   $hash{$first_name}->{$first_name}->{SNP}="-";
   $hash{$first_name}->{$first_name}->{indel}="-";
   $hash{$first_name}->{$first_name}->{gap}="-";
   $hash{$second_name}->{$second_name}->{SNP}="-";
   $hash{$second_name}->{$second_name}->{indel}="-";
   $hash{$second_name}->{$second_name}->{gap}="-";
   my $prefix1= $first_name.'_'.$second_name;
   my $prefix2= $second_name.'_'.$first_name;

   my $first_fasta=$outdir.'/'.$first_name.'_norepeats.fna';
   my $second_fasta=$outdir.'/'.$second_name.'_norepeats.fna';

   print "Running nucmer on $prefix1\n";
   my $nucmer_command1= "nucmer $options -p $prefix1 $first_fasta $second_fasta";
   if (system ($nucmer_command1)){die "Error running nucmer_command1 $nucmer_command1.\n";}

   my $filter_command1= "delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.snpfilter";
   if (system ($filter_command1)){die "Error running filter_command1 $filter_command1.\n";}

   my $snp_command1= "show-snps -CT $outdir/$prefix1.snpfilter > $outdir/$prefix1.snps";
   if (system ($snp_command1)){die "Error running snp_command1 $snp_command1.\n";}
   $snp_indel= `$bindir/SNP_INDEL_count.pl $outdir/$prefix1.snps`;
   $snp_indel=~ s/\n//;
   ($snp_n,$indel_n)= split /\t/,$snp_indel;    
   $hash{$first_name}->{$second_name}->{SNP}=$snp_n;
   $hash{$first_name}->{$second_name}->{indel}=$indel_n;
   
   my $filter_command2= "delta-filter -1 $identity $outdir/$prefix1.delta > $outdir/$prefix1.gapfilter";
   if (system ($filter_command2)){die "Error running filter_command2 $filter_command2.\n";}

   my $coords_command1= "show-coords -clTr $outdir/$prefix1.gapfilter > $outdir/$prefix1.coords";
   if (system ($coords_command1)){die "Error running coords_command1 $coords_command1.\n";}
   my $gaps1= `$bindir/parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix1.coords 2>/dev/null`;
   ($ref_gaps,$query_gaps,undef)= split /\n/,$gaps1;
   $hash{$first_name}->{$second_name}->{gap}=$ref_gaps;
   $hash{$second_name}->{$first_name}->{gap}=$query_gaps;

   print "Running nucmer on $prefix2\n";
   my $nucmer_command2= "nucmer $options -p $prefix2 $second_fasta $first_fasta";
   if (system ($nucmer_command2)){die "Error running nucmer_command2 $nucmer_command2.\n";}

   my $filter_command3= "delta-filter -1 $identity $outdir/$prefix2.delta > $outdir/$prefix2.snpfilter";
   if (system ($filter_command3)){die "Error running filter_command3 $filter_command3.\n";}

   my $snp_command2= "show-snps -CT $outdir/$prefix2.snpfilter > $outdir/$prefix2.snps";
   if (system ($snp_command2)){die "Error running snp_command2 $snp_command2.\n";}
   $snp_indel= `$bindir/SNP_INDEL_count.pl $outdir/$prefix2.snps`;
   $snp_indel=~ s/\n//;
   ($snp_n,$indel_n)= split /\t/,$snp_indel;    
   $hash{$second_name}->{$first_name}->{SNP}=$snp_n;
   $hash{$second_name}->{$first_name}->{indel}=$indel_n;

   my $filter_command4= "delta-filter $identity $outdir/$prefix2.delta > $outdir/$prefix2.gapfilter";
   if (system ($filter_command4)){die "Error running filter_command4 $filter_command4.\n";}

   my $coords_command2= "show-coords -clTr $outdir/$prefix2.gapfilter > $outdir/$prefix2.coords";
   if (system ($coords_command2)){die "Error running coords_command2 $coords_command2.\n";}
   
   my $gaps2= `$bindir/parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix2.coords 2>/dev/null`;

   my $check= `$bindir/checkNUCmer.pl -i $outdir/$first_name\_norepeats_$second_name\_norepeats.gaps -r $reference`;
   if ($check==1){print "$query aligned < 25% of the $reference genome\n";}
   $check= `$bindir/checkNUCmer.pl -i $outdir/$second_name\_norepeats_$first_name\_norepeats.gaps -r $query`;
   if ($check==1){print "$reference aligned < 25% of the $query genome\n";} 

   ($ref_gaps,$query_gaps,undef)= split /\n/,$gaps2;
    $hash{$first_name}->{$second_name}->{gap}=$ref_gaps;
    $hash{$second_name}->{$first_name}->{gap}=$query_gaps;

   $pm->finish(0,\%hash);
}
$pm-> wait_all_children;

&PrintMatrix(\%snp_hash,"SNPs");
&PrintMatrix(\%indel_hash,"INDELs");
&PrintMatrix(\%gap_hash,"Gaps");

print "NUCmer on all reference genomes complete.\n";
}

sub run_nucmer_with_ref
{
my $reference=shift;

my $pm= new Parallel::ForkManager($thread);
$pm->run_on_finish(sub{my ($pid,$exit_code,$ident,$exit_signal,$core_dump)=@_;});

foreach my $query(@query){
   $pm->start($query) and next;
   my ($first_name,$first_path,$first_suffix) = fileparse("$reference", qr/\.[^.]*/);
   my ($second_name,$second_path,$second_suffix) = fileparse("$query", qr/\.[^.]*/);
   $first_name =~ s/\.fna//;
   $second_name =~ s/\.fna//;
   next if ($first_name eq $second_name);

   my $prefix= $first_name.'_'.$second_name;

   my $first_fasta=$outdir.'/'.$first_name.'_norepeats.fna';
   my $second_fasta=$outdir.'/'.$second_name.'_norepeats.fna';

   print "Running nucmer on $prefix\n";
   my $nucmer_command1= "nucmer $options -p $prefix $first_fasta $second_fasta";
   if (system ($nucmer_command1)){die "Error running nucmer_command1 $nucmer_command1.\n";}

   my $filter_command1= "delta-filter -1 $identity $outdir/$prefix.delta > $outdir/$prefix.snpfilter";
   if (system ($filter_command1)){die "Error running filter_command1 $filter_command1.\n";}

   my $snp_command1= "show-snps -CT $outdir/$prefix.snpfilter > $outdir/$prefix.snps";
   if (system ($snp_command1)){die "Error running snp_command1 $snp_command1.\n";}
   $snp_indel= `$bindir/SNP_INDEL_count.pl $outdir/$prefix.snps`;
   $snp_indel=~ s/\n//;
   ($snp_n,$indel_n)= split /\t/,$snp_indel;    
#   print "$snp_indel\t$snp_n\t$indel_n\n";
   
   my $filter_command2= "delta-filter -1 $identity $outdir/$prefix.delta > $outdir/$prefix.gapfilter";
   if (system ($filter_command2)){die "Error running filter_command2 $filter_command2.\n";}

   my $coords_command1= "show-coords -clTr $outdir/$prefix.gapfilter > $outdir/$prefix.coords";
   if (system ($coords_command1)){die "Error running coords_command1 $coords_command1.\n";}
   my $gaps= `$bindir/parseGapsNUCmer.pl $gap_cutoff $outdir/$prefix.coords 2>/dev/null`;
   ($ref_gaps,$query_gaps,undef)= split /\n/,$gaps;

   my $check= `$bindir/checkNUCmer.pl -i $outdir/$first_name\_norepeats_$second_name\_norepeats.gaps -r $reference`;
   if ($check==1){print "$query aligned < 25% of the $reference genome\n";}
   #cleanup($second_name,$prefix);
   $pm->finish(0);
}
$pm-> wait_all_children;

print "NUCmer on all reference genomes complete.\n";
}


sub combo
{
my $by= shift;
return sub {()} if ! $by || $by =~ /\D/ || @_ < $by;
my @list=@_;

my @position=(0 .. $by - 2, $by - 2);
my @stop    =@list - $by .. $#list;
my $end_pos =$#position;
my $done    = undef;

return sub{
   return ()if $done;
   my $cur=$end_pos;
   {
      if (++$position[$cur] > $stop[$cur]){
         $position[--$cur]++;
         redo if $position[$cur] > $stop[$cur];
         my $new_pos=$position[$cur];
         @position[$cur .. $end_pos]=$new_pos .. $new_pos + $by;
      }
   }
   $done=1 if $position[0]==$stop[0];
   return @list[@position];
}
}

sub PrintMatrix
{
my ($hash_r, $title)=@_;
my %hash= %{$hash_r};

open (OUT, ">$outdir/$title.txt")||die "$!";
#print $title,"\n";
print OUT ("Ref\\Query\t",join("\t",sort keys %hash),"\n");
foreach my $ref(sort keys %hash){
   print OUT $ref,"\t";
   my $query_count=0;
   foreach my $query(sort keys %{$hash{$ref}}){
      $query_count++;
      print OUT $hash{$ref}->{$query};
      print OUT "\t" if ($query_count != scalar(@ARGV));
   }
   print OUT "\n";
}
#print "\n";
}

sub cleanup
{
#if (-e "$outdir/*.coords"){`rm $outdir/*.coords`};
if (-e "$outdir/*.mgaps"){unlink "$outdir/*.mgaps"};
if (-e "$outdir/*.ntref"){unlink "$outdir/*.ntref"};
`cat $outdir/*_repeat_stats.txt > $outdir/repeat_stats.txt`;
`mv $outdir/*.snps $outdir/snps`;
`mv $outdir/*.gaps $outdir/gaps`;
`mv $outdir/*.txt $outdir/stats`;
`mv $outdir/*.coords $outdir/stats`;
`mv $outdir/*norepeats.fna $outdir/temp`;
}

sub getBinDirectory
{
my @t = split '/', "$FindBin::RealBin";
my $path = join '/', @t;
return ($path);
}

sub usage
{
print <<Usage;
$0 -q <query_fasta> -d <output_directory> -t <# threads>

Usage
exit;
}
