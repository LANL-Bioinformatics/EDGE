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
my ($pairedReadsFile, $unpairedReadsFile, $outDir, $ref);
my $bwamemOpts="-T 50"; 
my $numCPU=4; #number of threads [4]
my $prefix="host_clean";
my $outputHost=0;
my $ref="/users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta";
my $outFasta=0;

GetOptions( 
            'p=s'   => \$pairedReadsFile,  
            'u=s'   => \$unpairedReadsFile,
            'prefix=s' => \$prefix,
            'ref=s' => \$ref,
            'host'  => \$outputHost,
            'fasta' => \$outFasta,
            'cpu=i'   => \$numCPU, # bwa option
            'bwaMemOptions=s'   => \$bwamemOpts,    # bwa mem options
            'o=s'   => \$outDir,
            'help|?'   => sub{&Usage(1)}
);
unless ($pairedReadsFile || $unpairedReadsFile) {&Usage}; 
unless ($ref && $outDir) {&Usage};
&checkFiles($ref,$pairedReadsFile,$unpairedReadsFile);
sub Usage
{
     my $printBwaMem=shift;
     print <<"END";
 Usage: perl $0 [options] -p 'reads1.fastq reads2.fastq' -ref host.fasta -o out_directory
        Input File:
        -ref          Host sequences in fasta  [default: human, if not provide, it will use
                      /users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta]

        -u            Unpaired reads, Single end reads
       
        -p            Paired reads in two files and separate by space in quote
                      Or one file in interleaved format

        Output:
        -o            Output directory.
 
        -prefix       Output File Name Prefix [host_clean] 
 
        -fasta        <boolean> Output in fasta format instead of fastq.
 
        -host         <boolean> Output Host id list.   
  
        Options:
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

my ($numTotalReads,$numUnalignedReads)=&runMapping($ref,$pairedReadsFile,$unpairedReadsFile,$outDir);

exit(0);

sub checkFiles
{
    my $refFile=shift;
    my $queryPairedFile=shift;
    my $queryUnpairedFile=shift;
    if (defined $queryPairedFile)
    {
        my ($mate1,$mate2) = split /\s+/,$queryPairedFile;
        if ( -z $mate1 || ! -e $mate1) { die "$mate1 file is empty";}
        if (defined $mate2)
        {
            if ( -z $mate2 || ! -e $mate2) { die "$mate2 file is empty";}
        }
    }
    if ( -z $refFile || ! -e $refFile) { die "$refFile file is empty";}
    if (defined $queryUnpairedFile)
    {
       if ( -z $queryUnpairedFile || ! -e $queryUnpairedFile) { die "$queryUnpairedFile file is empty";}
    }
    
}

sub runMapping 
{
    my $refFile=shift;
    my $queryPairedFile=shift;
    my $queryUnpairedFile=shift;
    my $outputDir=shift;
    my $mappingLogFile="$outputDir/$prefix.mapping.log";
    my $mappingStatsFile="$outputDir/$prefix.stats.txt";
    my $unalignedNonPairedFile="$outputDir/$prefix.unpaired.fastq";
    my $unalignedMate1File = "$outputDir/$prefix.1.fastq";
    my $unalignedMate2File = "$outputDir/$prefix.2.fastq";
    my $hostIdList= "$outputDir/$prefix.hostId.txt";
    my $numUnmappedPairedFile=0; # non host reads number
    my $numUnmappedUnpairedFile=0; # non host reads number
    my $numTotalReadsPaired=0;
    my $numTotalReadsUnpaired=0;
    my $numTotalReads=0;
    unlink $unalignedNonPairedFile if ( -s $unalignedNonPairedFile);
    open (my $stat_fh, ">$mappingStatsFile") or die "$! $mappingStatsFile\n";
    my $command = "bwa index $refFile";
    if (! -e "$refFile.bwt")
    {
        #print colored ("Indexing the host sequences",'yellow'), "\n"; 
       &executeCommand($command);
    }

    #print colored ("Running reads mapping to host sequence ... and log file: $mappingLogFile",'yellow'),"\n";
    if ($queryPairedFile)
    {
        #my ($mate1,$mate2) = split /s+/,$queryPairedFile;
        #if ( -z $mate1 || !$mate1) { die "One of the paired file is empty";}
        #if ( -z $mate2 || !$mate2) { die "One of the paired file is empty";}
        print $stat_fh $queryPairedFile,"\n";
        if ( -s $queryPairedFile)
        {
           $bwamemOpts .= " -p ";
        }
        $command = "bwa mem $bwamemOpts -t $numCPU $refFile $queryPairedFile 2>$mappingLogFile| ";
        print " $command\n"; 
        open (my $fh, "$command") or die "$! bwa mem command failed\n";
        open (my $unalignedMate1_fh, ">$unalignedMate1File") or die "$! $unalignedMate1File";
        open (my $unalignedMate2_fh, ">$unalignedMate2File") or die "$! $unalignedMate2File";
        open (my $unalignedNonPaired_fh, ">$unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
        open (my $hostId_fh, ">$hostIdList") or die "$hostIdList $!\n" if ($outputHost);
        while (<$fh>)
        {
            chomp;
            next if (/^\@/);
            my @samFields=split /\t/,$_;
            if ($samFields[10] eq "*" and !$outFasta)
            {
                $samFields[10] = "f" x length($samFields[9]);
            }
              # bit operation [and] on the flag 
            if (($samFields[1] & 4) and ($samFields[1] & 8))  # both paired reads unmapped
            {
                if ($samFields[1] & 64) # the read is the first read in a pair
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedMate1_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedMate1_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                if ($samFields[1] & 128) # the read is the second read in a pair
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedMate2_fh  ">".$samFields[0]."/2\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedMate2_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                $numUnmappedPairedFile++;
            }
            elsif($samFields[1] & 4)  # query is unmapped
            {
                if ($samFields[1] & 64)
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedNonPaired_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                if ($samFields[1] & 128)
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedNonPaired_fh  ">".$samFields[0]."/2\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedNonPaired_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                $numUnmappedPairedFile++;
            }
            else  #mapped reads
            {
                #$numMapped++;
                print $hostId_fh $samFields[0],"\t",$samFields[2],"\t",$samFields[5],"\n" if ($outputHost); 
                #print $hostId_fh $samFields[0],"\n" if ($outputHost); 
            }
        }
        close $fh;
        close $unalignedMate1_fh;
        close $unalignedMate2_fh;
        close $unalignedNonPaired_fh;
        close $hostIdList if ($outputHost);
        $numTotalReadsPaired= &parseMappingLog($mappingLogFile);
        print $stat_fh " Total reads: $numTotalReadsPaired\n";
        if ($numTotalReadsPaired)
        {
          printf $stat_fh (" Non host reads: %d (%.2f %%)\n", $numUnmappedPairedFile, $numUnmappedPairedFile/$numTotalReadsPaired *100);
          printf $stat_fh (" Host reads: %d (%.2f %%)\n", $numTotalReadsPaired - $numUnmappedPairedFile, ($numTotalReadsPaired-$numUnmappedPairedFile)/$numTotalReadsPaired *100);
        }
    }

    if ( -s $queryUnpairedFile)
    {
        print $stat_fh $queryUnpairedFile,"\n";
        $command = "bwa mem $bwamemOpts -t $numCPU $refFile $queryUnpairedFile 2>$mappingLogFile| ";
        print " $command\n"; 
        open (my $fh, "$command") or die "$! bwa mem command failed\n";
        open (my $unalignedNonPaired_fh, ">> $unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
        open (my $hostId_fh, ">>$hostIdList") or die "$hostIdList $!\n" if ($outputHost);
        while (<$fh>)
        {
            chomp;
            next if (/^\@/);
            my @samFields=split /\t/,$_;
            if ($samFields[10] eq "*" and !$outFasta)
            {
                $samFields[10] = "f" x length($samFields[9]);
            }

            if ($samFields[1] & 4)
            {
                if ($samFields[1] & 64)
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedNonPaired_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                if ($samFields[1] & 128)
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedNonPaired_fh  ">".$samFields[0]."/2\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedNonPaired_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                else 
                {
                    $samFields[0] =~ s/\/\d$//;
                    if ($outFasta)
                    {
                        print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                    }
                    else
                    {
                        print $unalignedNonPaired_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                }
                $numUnmappedUnpairedFile++;

            }
            else  # mapped reads
            {
                #$numMapped++;
                print $hostId_fh $samFields[0],"\t",$samFields[2],"\t",$samFields[5],"\t",$samFields[11],"\t",$samFields[12],"\n" if ($outputHost); 
            }
        }
        close $fh;
        close $unalignedNonPaired_fh;
        close $hostId_fh if ($outputHost);
        $numTotalReadsUnpaired=&parseMappingLog($mappingLogFile);
        print $stat_fh " Total reads: $numTotalReadsUnpaired\n";
        if ($numTotalReadsUnpaired)
        {
            printf $stat_fh (" Non host reads: %d (%.2f %%)\n", $numUnmappedUnpairedFile, $numUnmappedUnpairedFile/$numTotalReadsUnpaired *100);
            printf $stat_fh (" Host reads: %d (%.2f %%)\n", $numTotalReadsUnpaired - $numUnmappedUnpairedFile, ($numTotalReadsUnpaired-$numUnmappedUnpairedFile)/$numTotalReadsUnpaired *100);
        }
    }
    close $stat_fh;
    system "cat $mappingStatsFile";
    return ($numTotalReadsUnpaired+$numTotalReadsPaired , $numTotalReadsUnpaired + $numTotalReadsPaired);
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


  
