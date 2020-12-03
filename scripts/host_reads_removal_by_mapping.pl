#!/usr/bin/env perl
#Purpose: remove host reads
#     Human genome index for bwa > 0.6 is at /users/218819/scratch/data/databases/human_chromosomes/all_chromoaome.fasta
# required bwa > 0.7

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$Bin/../bin/:$ENV{PATH}";

umask 000;
my (@pairedReadsFile, @unpairedReadsFile, $outDir, $ref);
my $bwamemOpts="-T 50"; 
my $numCPU=4; #number of threads [4]
my $prefix="host_clean";
my $outputHost=0;
my $ref="/users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta";
my $outFasta=0;
my $similarity_cutoff=90;

GetOptions( 
            'p=s{,}'   => \@pairedReadsFile,  
            'u=s{,}'   => \@unpairedReadsFile,
            'prefix=s' => \$prefix,
            'ref=s' => \$ref,
            'host'  => \$outputHost,
            'fasta' => \$outFasta,
            'similarity|s=f'  => \$similarity_cutoff,
            'cpu=i'   => \$numCPU, # bwa option
            'bwaMemOptions=s'   => \$bwamemOpts,    # bwa mem options
            'o=s'   => \$outDir,
            'help|?'   => sub{&Usage("",1)}
);
unless (@pairedReadsFile || @unpairedReadsFile) {&Usage("Missing Input Files")}; 
unless ($ref && $outDir) {&Usage("Missing and Output directory")};
@pairedReadsFile = &checkFiles();
sub Usage
{
     my $msg=shift;
     my $printBwaMem=shift;
     print "\n     ".$msg."\n\n" if $msg;
     print <<"END";
 Usage: perl $0 [options] -p reads1.fastq reads2.fastq -ref host.fasta -o out_directory
        Input File:
        -ref          Host sequences in fasta  [default: human, if not provide, it will use
                      /users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta]

        -u            Unpaired reads, Single end reads
       
        -p            Paired reads in two files and separate by space
                      Or one file in interleaved format

        Output:
        -o            Output directory.
 
        -prefix       Output File Name Prefix [host_clean] 
 
        -fasta        <boolean> Output in fasta format instead of fastq.
 
        -host         <boolean> Output Host reads file.   
  
        Options:
        -similarity   <NUM > alignment similarity percentage (=~ matched/read_length*100) [default 90]
        
        -bwaMemOptions   <String in quote> see "bwa mem" options with -h flag
                         ex: '-T 60' 
                         
        -cpu          number of CPUs [4]
 
        -h            print "bwa mem" options
END
 
if ($printBwaMem)
{ 
print <<"END";

bwa mem algorithm options:

       -k INT     minimum seed length [19]
       -w INT     band width for banded alignment [100]
       -d INT     off-diagonal X-dropoff [100]
       -r FLOAT   look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
       -c INT     skip seeds with more than INT occurrences [10000]
       -P         skip pairing; perform mate SW only
       -A INT     score for a sequence match [1]
       -B INT     penalty for a mismatch [4]
       -O INT     gap open penalty [6]
       -E INT     gap extension penalty; a gap of size k cost {-O} + {-E}*k [1]
       -L INT     penalty for clipping [5]
       -U INT     penalty for an unpaired read pair [9]
       -T INT     minimum score to output [30]
       -H         hard clipping

END
}
exit;
}

mkdir $outDir if ( ! -e $outDir);

my ($numTotalReads,$numUnalignedReads)=&runMapping($ref,$outDir);

exit(0);

sub checkFiles
{
    unless (-s "$ref.bwt")  # bwa index is existed. no need to check ref fasta file.
    {
        if (is_file_empty($ref)) { die "$ref is empty";}
    }
 
    my %file;
    my @make_paired_paired_files;
    if (@pairedReadsFile)
    {
        if (scalar(@pairedReadsFile) % 2) { Usage("Please check paired data input are even file numbers\n") ;}
        map { if(is_file_empty($_)){ Usage("Please check paired data input at flag -p.\n    $_ doesn't not exist or empty.");} $file{basename($_)}=1; } @pairedReadsFile;
        #make pair in a new array 'read1_1 read1_2', 'read2_1 read2_2' ...
        for(my$i=0;$i<=$#pairedReadsFile;$i=$i+2)
        {
            if (&is_paired($pairedReadsFile[$i], $pairedReadsFile[$i+1]))
            {
                push @make_paired_paired_files, "$pairedReadsFile[$i] $pairedReadsFile[$i+1]";
            }
            else
            {
                print "The sequence names of the paired end reads in $pairedReadsFile[$i],$pairedReadsFile[$i+1] are not matching.\nWill use them as single end reads\n";
                push @unpairedReadsFile, $pairedReadsFile[$i],$pairedReadsFile[$i+1];
                delete $file{basename($pairedReadsFile[$i])};
                delete $file{basename($pairedReadsFile[$i+1])};
            }
        }
    }
    if (@unpairedReadsFile)
    {
        map { if(is_file_empty($_))
              { 
                  Usage("Please check unpaired data input at flag -u.\n    $_ doesn't not exist or empty.");
              } 
              if ($file{basename($_)}) 
              {
                  Usage("The single end file, $_,has been used in the paired end data.")
              }
            } @unpairedReadsFile;
    }
    return (@make_paired_paired_files);
}

sub runMapping 
{
    my $refFile=shift;
    my $outputDir=shift;
    my $mappingLogFile="$outputDir/$prefix.mappingPE.log";
    my $mappingLogFile2="$outputDir/$prefix.mappingSE.log";
    my $mappingStatsFile="$outputDir/$prefix.stats.txt";
    my $unalignedNonPairedFile="$outputDir/$prefix.unpaired.fastq";
    my $unalignedMate1File = "$outputDir/$prefix.1.fastq";
    my $unalignedMate2File = "$outputDir/$prefix.2.fastq";
    my $host_fastq= "$outputDir/host.fastq";
    my $numUnmappedPairedFile=0; # non host reads number
    my $numUnmappedUnpairedFile=0; # non host reads number
    my $numTotalReadsPaired=0;
    my $numTotalReadsUnpaired=0;
    my $numTotalReads=0;
    unlink $mappingLogFile;
    unlink $mappingLogFile2;
    unlink $unalignedNonPairedFile;
    unlink $unalignedMate1File;
    unlink $unalignedMate2File;
    unlink $host_fastq;
    open (my $stat_fh, ">$mappingStatsFile") or die "$! $mappingStatsFile\n";
    my $print_string;
    print $stat_fh "Remove Reads from Host ". basename($ref)."\n";
    my $command = "bwa index $refFile";
    if (! -e "$refFile.bwt")
    {
        #print colored ("Indexing the host sequences",'yellow'), "\n"; 
       &executeCommand($command);
    }

    #print colored ("Running reads mapping to host sequence ... and log file: $mappingLogFile",'yellow'),"\n";
    if (@pairedReadsFile)
    {
      foreach my $queryPairedFile(@pairedReadsFile)
      {
        #my ($mate1,$mate2) = split /s+/,$queryPairedFile;
        #if ( -z $mate1 || !$mate1) { die "One of the paired file is empty";}
        #if ( -z $mate2 || !$mate2) { die "One of the paired file is empty";}
        print $stat_fh $queryPairedFile,"\n";
        if ( -s $queryPairedFile)
        {
           $bwamemOpts .= " -p ";
        }
        $command = "bwa mem $bwamemOpts -t $numCPU $refFile $queryPairedFile 2>>$mappingLogFile| ";
        print " $command\n"; 
        open (my $fh, "$command") or die "$! bwa mem command failed\n";
        open (my $unalignedMate1_fh, ">>$unalignedMate1File") or die "$! $unalignedMate1File";
        open (my $unalignedMate2_fh, ">>$unalignedMate2File") or die "$! $unalignedMate2File";
        open (my $unalignedNonPaired_fh, ">>$unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
        open (my $host_fh, ">$host_fastq") or die "$host_fastq $!\n" if ($outputHost);
        while (<$fh>)
        {
            chomp;
            next if (/^\@/);
            my @samFields=split /\t/,$_;
            next if ($samFields[1] & 256 or $samFields[1] & 2048); # skip secondary alignment
            my $CIGAR=$samFields[5];
            my (@CIGARnum)=  $CIGAR =~ /(\d+)/g;
            my (@CIGARtag)=  $CIGAR =~ /([A-Z])/g;
  
            if ($samFields[10] eq "*" and !$outFasta)
            {
                $samFields[10] = "f" x length($samFields[9]);
            }
              # bit operation [and] on the flag 
            if (($samFields[1] & 4) and ($samFields[1] & 8))  # both paired reads unmapped
            {
                $samFields[0] =~ s/\/\d$//;
                if ($samFields[1] & 64) # the read is the first read in a pair
                {   
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedMate1_fh);
                }
                if ($samFields[1] & 128) # the read is the second read in a pair
                {
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedMate2_fh);
                }
                $numUnmappedPairedFile++;
            }
          #  elsif($samFields[1] & 4)  # query is unmapped
          #  {
          #      &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedNonPaired_fh);
          #      $numUnmappedPairedFile++;
          #  }
            else  #mapped reads  Host reads
            {
              
                my @samFields2=split /\t/,<$fh>;
                while ($samFields2[1] & 256 or $samFields2[1] & 2048)
                {
                    @samFields2=split /\t/,<$fh>;
                }
                my $CIGAR2=$samFields2[5];
                my (@CIGARnum2)=  $CIGAR2 =~ /(\d+)/g;
                my (@CIGARtag2)=  $CIGAR2 =~ /([A-Z])/g;
                my $read_similarity1=&get_similarity(\@CIGARnum,\@CIGARtag,$samFields[9]);
                my $read_similarity2=&get_similarity(\@CIGARnum2,\@CIGARtag2,$samFields[9]);
                
                if ($read_similarity1 < $similarity_cutoff and $read_similarity2 < $similarity_cutoff)
                {
                    if ($samFields[1] & 64) # the read is the first read in a pair
                    {   
                        &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedMate1_fh);
                    }
                    if ($samFields[1] & 128) # the read is the second read in a pair
                    {
                        &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedMate2_fh);
                    }
                    if ($samFields2[1] & 64) # the read is the first read in a pair
                    {   
                        &print_fastq($samFields2[1],$samFields2[0],$samFields2[9],$samFields2[10],$unalignedMate1_fh);
                    }
                    if ($samFields2[1] & 128) # the read is the second read in a pair
                    {
                        &print_fastq($samFields2[1],$samFields2[0],$samFields2[9],$samFields2[10],$unalignedMate2_fh);
                    }
                    $numUnmappedPairedFile += 2;
                }
                elsif ($read_similarity1 < $similarity_cutoff and $read_similarity2 >= $similarity_cutoff )
                {
                   # print $samFields[0],"\t",$read_similarity1, " \t2 ", $read_similarity2,"\n";
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedNonPaired_fh);
                    &print_fastq($samFields2[1],$samFields2[0],$samFields2[9],$samFields2[10],$host_fh) if ($outputHost);
                    $numUnmappedPairedFile++;
                    #$numMapped++;
                }
                elsif ($read_similarity2 < $similarity_cutoff and $read_similarity1 >= $similarity_cutoff)
                {
                    #print $samFields[0],"\t",$read_similarity1, " 1\t ", $read_similarity2,"\n";
                    &print_fastq($samFields2[1],$samFields2[0],$samFields2[9],$samFields2[10],$unalignedNonPaired_fh);
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$host_fh) if ($outputHost);
                    $numUnmappedPairedFile++;
                    #$numMapped++;
                }
                else   ## host reads 
                {
                    #print $samFields[0],"\t",$read_similarity1, " 1\t2 ", $read_similarity2,"\n";
                    #$numMapped+2;
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$host_fh) if ($outputHost); 
                    &print_fastq($samFields2[1],$samFields2[0],$samFields2[9],$samFields2[10],$host_fh) if ($outputHost); 
                }
            }
        }
        close $fh;
        close $unalignedMate1_fh;
        close $unalignedMate2_fh;
        close $unalignedNonPaired_fh;
        close $host_fh if ($outputHost);
      }# foreach (@pairedReadsFile)
        $numTotalReadsPaired= &parseMappingLog($mappingLogFile);
        $print_string .= " Total paired reads: $numTotalReadsPaired\n";
        if ($numTotalReadsPaired)
        {
          $print_string .= sprintf  (" Non host reads: %d (%.2f %%)\n", $numUnmappedPairedFile, $numUnmappedPairedFile/$numTotalReadsPaired *100);
          $print_string .= sprintf  (" Host reads: %d (%.2f %%)\n", $numTotalReadsPaired - $numUnmappedPairedFile, ($numTotalReadsPaired-$numUnmappedPairedFile)/$numTotalReadsPaired *100);
        }
    }

    if (@unpairedReadsFile)
    {
      foreach my $queryUnpairedFile (@unpairedReadsFile)
      {
        print $stat_fh $queryUnpairedFile,"\n";
        $command = "bwa mem $bwamemOpts -t $numCPU $refFile $queryUnpairedFile 2>$mappingLogFile2| ";
        print " $command\n"; 
        open (my $fh, "$command") or die "$! bwa mem command failed\n";
        open (my $unalignedNonPaired_fh, ">> $unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
        open (my $host_fh, ">>$host_fastq") or die "$host_fastq $!\n" if ($outputHost);
        while (<$fh>)
        {
            chomp;
            next if (/^\@/);
            my @samFields=split /\t/,$_;
            next if ($samFields[1] & 256 or $samFields[1] & 2048); # skip secondary alignment
            my $CIGAR=$samFields[5];
            my (@CIGARnum)=  $CIGAR =~ /(\d+)/g;
            my (@CIGARtag)=  $CIGAR =~ /([A-Z])/g;
            if ($samFields[10] eq "*" and !$outFasta)
            {
                $samFields[10] = "f" x length($samFields[9]);
            }

            if ($samFields[1] & 4) # query is unmapped
            {
                &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedNonPaired_fh);
                $numUnmappedUnpairedFile++;
            }
            else  # mapped reads
            {
                my $read_similarity=&get_similarity(\@CIGARnum,\@CIGARtag,$samFields[9]);
                if ($read_similarity >= $similarity_cutoff)
                {
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$host_fh) if ($outputHost); 	
                    #$numMapped++;
                }
                else
                {
                    &print_fastq($samFields[1],$samFields[0],$samFields[9],$samFields[10],$unalignedNonPaired_fh);
                    $numUnmappedUnpairedFile++;
                }
            }
        }
        close $fh;
        close $unalignedNonPaired_fh;
        close $host_fh if ($outputHost);
      } # foreach (@unpairedReadsFile)
        $numTotalReadsUnpaired=&parseMappingLog($mappingLogFile2);
        $print_string .= " Total Unpaired reads: $numTotalReadsUnpaired\n";
        if ($numTotalReadsUnpaired)
        {
            $print_string .= sprintf (" Non host reads: %d (%.2f %%)\n", $numUnmappedUnpairedFile, $numUnmappedUnpairedFile/$numTotalReadsUnpaired *100);
            $print_string .= sprintf (" Host reads: %d (%.2f %%)\n", $numTotalReadsUnpaired - $numUnmappedUnpairedFile, ($numTotalReadsUnpaired-$numUnmappedUnpairedFile)/$numTotalReadsUnpaired *100);
        }
    }
    if (@unpairedReadsFile && @pairedReadsFile)
    {
        print $stat_fh $print_string;
    }
    
    
        my $totalReads = $numTotalReadsUnpaired + $numTotalReadsPaired; 
        my $totalUnmapped = $numUnmappedUnpairedFile + $numUnmappedPairedFile;
        printf $stat_fh ("\nTotal reads: %d\n",  $totalReads);
        printf $stat_fh (" Total Non host reads: %d (%.2f %%)\n", $totalUnmapped, $totalUnmapped/$totalReads *100);
        printf $stat_fh (" Total Host reads: %d (%.2f %%)\n", $totalReads - $totalUnmapped, ($totalReads-$totalUnmapped)/$totalReads *100);
    
    close $stat_fh;
    system "cat $mappingStatsFile";
    return ($numTotalReadsUnpaired+$numTotalReadsPaired , $numTotalReadsUnpaired + $numTotalReadsPaired);
}

sub get_similarity
{
     my $CIGARnum_r=shift;
     my $CIGARtag_r=shift;
     my $seq=shift;
     my @CIGARnum=@{$CIGARnum_r};
     my @CIGARtag=@{$CIGARtag_r};
     my $aligned_pos=0;
     my $total_pos=0;
     foreach my $index (0..$#CIGARnum)
     {
        $aligned_pos += $CIGARnum[$index] if ($CIGARtag[$index] eq "M");
        $total_pos += $CIGARnum[$index];
     }
     $total_pos = length($seq) if (!$total_pos); #avoid 0
     $total_pos = 1 if (!$total_pos); #avoid 0
     return $aligned_pos/$total_pos*100;
}
sub print_fastq 
{
    my $flag=shift;
    my $id=shift;
    my $seq=shift;
    my $qual=shift;
    my $fh =shift;
    $id=~ s/\/\d$//;
    $id .= ($flag & 128)? "/2" : "/1";
    
    ($outFasta)? 
    print $fh  ">".$id."\n".$seq."\n":
    print $fh  "@".$id."\n".$seq."\n+\n".$qual."\n";
}

sub parseMappingLog 
{
    my $log=shift;
    my $numReads;
    open (my $fh, $log) or die "$! open $log failed\n";
    while (<$fh>)
    {  
        if ($_=~ /read (\d+) sequences/)
        {
           $numReads += $1;
        }
    }
    close $fh;
    return $numReads;
}

sub executeCommand 
{
    my $command = shift;
    system($command) == 0
         || die "the command $command failed\n";
}
sub open_file
{
    my ($file) = @_;
    my $fh;
    my $pid;
    if ( $file=~/\.gz$/i ) { $pid=open($fh, "gunzip -c $file |") or die ("gunzip -c $file: $!"); }
    else { $pid=open($fh,'<',$file) or die("$file: $!"); }
    $SIG{'PIPE'} = 'IGNORE';
    return ($fh,$pid);
}
sub is_paired
{
    my $paired_1=shift;
    my $paired_2=shift;
    my ($fh1,$pid1)=open_file($paired_1);
    my ($fh2,$pid2)=open_file($paired_2);
    my $count=0;
    my $check_num=1000;
    my $is_paired = 1;
    ## check top 1000 sequences paired by matching names.
    for ($count..$check_num)
    {
        my $id1=<$fh1>;
        my $seq1=<$fh1>;
        my $q_id1=<$fh1>;
        my $q_seq1=<$fh1>;
        my $id2=<$fh2>;
        my $seq2=<$fh2>;
        my $q_id2=<$fh2>;
        my $q_seq2=<$fh2>;
        my ($name1) = $id1 =~ /(\S+)/; 
        $name1 =~ s/\.\d$//;
        $name1 =~ s/\/\d$//;
        my ($name2) = $id2 =~ /(\S+)/; 
        $name2 =~ s/\.\d$//;
        $name2 =~ s/\/\d$//;
        if ($name1 ne $name2)
        {
            $is_paired=0;
        }
    }
    close $fh1;
    close $fh2;
    $SIG{'PIPE'} = 'DEFAULT';
    kill 9, $pid1; # avoid gunzip broken pipe
    kill 9, $pid2; # avoid gunzip broken pipe
    return $is_paired;
}
  
sub is_file_empty 
{
    #check file exist and non zero size
    my $file=shift;
    my $empty=1;
    if (-e $file) {$empty=0};
    if (-z $file) {$empty=1};
    return $empty;
}
