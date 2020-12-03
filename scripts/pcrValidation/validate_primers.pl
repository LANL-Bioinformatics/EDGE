#!/usr/bin/perl -w 
#This script validates primers used to amplify a target
#Author: Chien-Chi Lo at lanl dot gov  20130918


use strict;
use FindBin qw($Bin);
use Getopt::Long;
use Cwd;
$ENV{PATH} = "$Bin:$Bin/../../bin/:$Bin/../../scripts:$ENV{PATH}";

my $working_dir = Cwd::getcwd;

my $outDir;
my $ref_file;
my $primer_file;
my $numMAX_MISMATCHES=1;
my $numCPUs=4;
my $version=0.1;
my $prefix="primerMapping";
GetOptions(
            "output|o=s"       => \$outDir,
            "ref|r=s"          => \$ref_file,
            "primer|p=s"       => \$primer_file,
            "mismatch|m=i"     => \$numMAX_MISMATCHES,
            "threads|t=i"      => \$numCPUs,
            "prefix=s"         => \$prefix,
            "version"          => sub{print "Version $version\n";exit;},
            "help|?"           => sub{Usage()} );   
     
if (! $ref_file or ! $primer_file) {&Usage();}

if (! $outDir) {$outDir=$working_dir;}
system ("mkdir -p $outDir");

my $pRef= &readFastaSeq($ref_file);  # reference pointer  ->{seq}    ->{len}
     #read the primer file
my ($pPrimer,$primerRenameFile)= &readPrimerSeq($primer_file);

# index reference 
if ( ! -e "$ref_file.bwt")
{
    system("bwa index $ref_file 1>/dev/null");
}

# run bwa mapping
system("bwa mem -t $numCPUs -p -k 5 -A 1 -B 1 -O 3 -E 1 -T 5 -v 1 $ref_file $primerRenameFile > $outDir/$prefix.sam 2>/dev/null");
system("samtools view -HS $outDir/$prefix.sam > $outDir/$prefix.name.sam" );
system("rm $ref_file.*") if (-e "$ref_file.bwt");
        
open (my $sam_fh, "$outDir/$prefix.sam") or die "Cannot read $outDir/$prefix.sam\n";
open (my $sam_ofh, ">>", "$outDir/$prefix.name.sam") or die "Cannot write $outDir/$prefix.name.sam\n";
#p001    97      gi|49183039|ref|NC_005945.1|    4524337 22      20M     =       4525115 798     TATCGTGCACCAACTCCACC    *       NM:i:0  AS:i:20 XS:i:15
#p001    145     gi|49183039|ref|NC_005945.1|    4525115 36      20M     =       4524337 -798    TCGTCGCCTTATCAGCACTC    *       NM:i:0  AS:i:20 XS:i:12

my $reason;
my %result;
my $fail_flag=0;
while(<$sam_fh>)
{
    chomp;
    next if (/^\@/);
    my $primerId;
    my $ref_seq;
    my $primer_seq;
    my @samFields =split /\t/,$_;
    my @samFields_name_replace = @samFields;
    next if ( $samFields[1] & 2048 ); #Chimera
    next if ( $samFields[1] & 256 ); #not primary alignment
    # not primary alignment
  #  my $next=<$sam_fh>;
   # chomp $next;
   # my @samFields2 =split /\t/,$next;
   # my ($numMisMatch) = $samFields1[11] =~ m/NM:i:(\d+)/;
    if ($samFields[1] & 64)  # forward primer
    {
       $reason="";
       $fail_flag=0;
       $primerId=$pPrimer->{"$samFields[0]/1"}->{oldID};
       $ref_seq=substr($pRef->{$samFields[2]}->{seq},$samFields[3]-1,$pPrimer->{"$samFields[0]/1"}->{len});
       $primer_seq=$pPrimer->{"$samFields[0]/1"}->{seq};
       $primer_seq=&ReverseComplement($primer_seq) if ($samFields[1] &16);
       $samFields_name_replace[0] = $primerId;
       print $sam_ofh join("\t",@samFields_name_replace),"\n";
    }
    if ($samFields[1] & 128) # reverse primer
    {
       $primerId=$pPrimer->{"$samFields[0]/2"}->{oldID};
       $ref_seq=substr($pRef->{$samFields[2]}->{seq},$samFields[3]-1,$pPrimer->{"$samFields[0]/2"}->{len});
       $primer_seq=$pPrimer->{"$samFields[0]/2"}->{seq};
       $primer_seq=&ReverseComplement($primer_seq) if ($samFields[1] &16);
       $samFields_name_replace[0] = $primerId;
       print $sam_ofh join("\t",@samFields_name_replace),"\n";
    }
    # primer not match
    if ($samFields[1] &4)
    {  
        $fail_flag=1;
        $reason .="Primer $primerId fail to align to input sequences\n";
    }
    elsif ($samFields[5] =~ /I|D/)
    {
        $fail_flag=1;
        $reason .="Primer $primerId has INDELs ($samFields[5]) when aligning to input sequence ($samFields[2])\n";
        $reason .=" r    $ref_seq\n p    $primer_seq\n";
    }
    else
    {
        my ($alignment,$numMis)=&checkAlign($primer_seq,$ref_seq); 

       # if ($numMisMatch1 > $numMAX_MISMATCHES)
       # {
            if ($numMis > $numMAX_MISMATCHES)
            {
                $fail_flag=1;
                $reason .="Primer $primerId has $numMis mismatches to the input sequence ($samFields[2])\n";
                $reason .=" r    $ref_seq\n      $alignment\n p    $primer_seq\n\n";
            }
            else
            {
                $reason .="Primer $primerId Alignment to ($samFields[2])\n" . " r    $ref_seq\n      $alignment\n p    $primer_seq\n";
            }
       # }
    }
   # if (!$fail_flag) # Success
   # {
   #         $reason .="Primer $primerId Alignment to the reference ($samFields[2])\n" . " r    $ref_seq\n      $alignment\n p    $primer_seq\n";
   # }

    if ($samFields[1] & 128) # reverse primer
    {
        $result{$samFields[0]} = "Primer Pair ". $pPrimer->{"$samFields[0]/1"}->{oldID} . " and " . $pPrimer->{"$samFields[0]/2"}->{oldID}. "\n";
        if($fail_flag)
        {
            $result{$samFields[0]} .= "PCR failure! Because ... \n" . $reason ;
        }
        else
        {
            if ($samFields[6] ne "=")
            {
            	$result{$samFields[0]} .= "Alignment success! \n" . $reason;
                $result{$samFields[0]} .= " WARNINGS:  The primers pair match to different input sequences. ($samFields[2] vs $samFields[6])\n";
            }
            else
            {
            	$result{$samFields[0]} .= "PCR success! \n" . $reason;
                my ($start,$end,$amplicon_size,$direction);
                if ( $samFields[8] < 0 )
                {
                    $start=$samFields[7];
                    $end=$samFields[3]+$pPrimer->{"$samFields[0]/2"}->{len}-1;
                    $amplicon_size = $end - $start + 1;
                    $direction="Same Strand";
                }
                else
                {
                    $start=$samFields[3];
                    $end=$samFields[7]+$pPrimer->{"$samFields[0]/1"}->{len}-1;
                    $amplicon_size = $end - $start + 1;
                    $direction="Reverse Strand";
                }
                $result{$samFields[0]} .="The primers amplify $samFields[2] from $start to $end, with size $amplicon_size\n";
                if ($amplicon_size > 2000)
                {
                    $result{$samFields[0]} .= " WARNINGS:  The amplicon size is > 2000 bp.\n";
                }
            }
        }
    }
    
}    
close $sam_fh;
close $sam_ofh;

system("samtools view -huS $outDir/$prefix.name.sam | samtools sort -T $outDir -O BAM -o $outDir/$prefix.bam -" );
print $result{$_},"----------\n" foreach (sort keys %result);
unlink "$outDir/$prefix.sam";
#unlink "$outDir/$prefix.name.sam";
unlink $primerRenameFile;
exit(0);

## SUBs ##
sub checkAlign
{
    my $primer_seq=shift;
    my $ref_seq=shift;
    my @primer_seq= split //, $primer_seq;
    my @ref_seq= split //, $ref_seq;
    my $numMis=0;
    my $alignment;
    for my $i(0..$#primer_seq)
    {
       if ($primer_seq[$i] eq $ref_seq[$i])
       {
           $alignment.="|";
       }
       else
       {   
           if (&checkDegen($ref_seq[$i],$primer_seq[$i]))
           {
               $alignment.="+";
           }
           else
           {
               $alignment.=" ";
               $numMis++;
           }
       }
    }
    return ($alignment,$numMis);
}

sub checkDegen 
{
    my $ref_base=shift;
    my $query_base=shift;
    my %DEG=(U=>"T",
             M=>"AC",
             R=>"AG",
             W=>"AT",
             S=>"CG",
             Y=>"CT",
             K=>"GT",
             V=>"ACG",
             H=>"ACT",
             D=>"AGT",
             B=>"CGT",
             N=>"AGCT" 
            );
    my $match=0;
    if ($ref_base !~ /[ATCG]/i and $query_base !~ /[ATCG]/i)
    {
       foreach my $base1(split //,$DEG{$ref_base})
       {
           foreach my $base2(split //,$DEG{$query_base})
           {
                 $match =1 if ($base1 eq $base2);
           }
       }
    }
    elsif($ref_base !~ /[ATCG]/i)
    {
        $match=1 if ($DEG{$ref_base} =~ /$query_base/);
    }
    elsif($query_base !~ /[ATCG]/i)
    {
        $match=1 if ($DEG{$query_base} =~ /$ref_base/); 
    }
    return $match;
}


exit 0;

sub readPrimerSeq
{
    my $seqFile=shift;
    my %hash;
    my $output="$outDir/primers_rename.fa";
    open (my $ofh, ">$output") or die "Cannot write $output\n"; 
    open (my $fh, $seqFile) or die "$! $seqFile";
    $/ = ">";
    my $primer_id="p000";
    while (<$fh>)
    { 
        $_ =~ s/\>//g;
        my ($id1, $seq1) = split /\n/, $_;
        next if (!$id1);
        ($id1) =~ s/^(\S+).*/$1/;
        $primer_id++;
        #my $len = length($seq);
        print $ofh ">$primer_id/1\n$seq1\n";
        my $next_seq=<$fh>;
        my ($id2, $seq2) = split /\n/, $next_seq;
        ($id2) =~ s/^(\S+).*/$1/;
        print $ofh ">$primer_id/2\n$seq2\n";
        $hash{"$primer_id/1"}->{oldID}=$id1;
        $hash{"$primer_id/2"}->{oldID}=$id2;
        $hash{"$primer_id/1"}->{len}=length $seq1;
        $hash{"$primer_id/2"}->{len}=length $seq2;
        $hash{"$primer_id/1"}->{seq}=$seq1;
        $hash{"$primer_id/2"}->{seq}=$seq2;
    }
    close $fh;
    close $ofh;
    $/="\n";
    return (\%hash,$output);
}

sub readFastaSeq
{
    my $seqFile=shift;
    my %hash;
    open (my $fh, $seqFile) or die "$! $seqFile";
    $/ = ">";
    while (<$fh>)
    { 
        $_ =~ s/\>//g;
        my ($id, @seq) = split /\n/, $_;
        next if (!$id);
        ($id) =~ s/^(\S+).*/$1/;
        my $seq = join "", @seq;
        my $len = length($seq);
        $hash{$id}->{seq}=$seq;
        $hash{$id}->{len}=$len;
    }
    close $fh;
    $/="\n";
    return \%hash;
}

sub Usage
{
    print  <<USAGE;
perl $0 -ref <reference.fasta> -primer <primers.fa> -mismatch <num_MAX_MISMATCHES> -threads <num_of_CPU>

   -ref        the assembled genome i.e. a set of contigs
   -primer     the primers file in fasta format
   -mismatch   the max number of Mismatch between assembled genome and primer (default: 1)
   -prefix     Output prefix
   -threads    num_of_CPU:  bwa cpu usage
   -output     Output Directory 
 
USAGE
   #-prefix     Output file prefix (default: primerMapping)
exit;
}
 
sub ReverseComplement
{
        my $dna = $_[0];
        my $ReverseCompSeq = reverse ($dna);
        $ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
        return($ReverseCompSeq);
}
