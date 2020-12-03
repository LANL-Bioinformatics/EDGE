#!/usr/bin/env perl
# pcrUniquePrimer.pl is used to design unique primer pairs for input contigs.
#
# Required third-party tools:
# - Primer3 >=2.0
# - BWA >=0.7
# - BLAST+
#
# Po-E Li 20130820
#
# Update log:
# 20140519 - Fixing the bug if soft clip is prior to the start of the reference
# 20131002 - Showing the predicted Tm of the primer on background sequence
#

use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long;
use Tm_calculate;
use strict;

my $edge_home = $ENV{"EDGE_HOME"};
$ENV{PATH} = "$Bin:$edge_home/bin:$edge_home/scripts:$ENV{PATH}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'threads|t=s',
    'gff3|g=s',
    'cutoff|c=s',
    'no_tm_chk|n',
    'bg_tm_diff|f=s',
    'tm_max=i',
    'tm_opt=i',
    'tm_min=i',
    'len_max=i',
    'len_opt=i',
    'len_min=i',
    'top|p=s',
    'tmp=s',
    'debug|d',
    'help|?') || &Usage();

if ( $opt{help} ) { &Usage(); }

# DEBUG flag
my $DEBUG = defined $opt{debug} ? 1 : 0;

# tmp directory
my $tmp =  $opt{tmp} || '/tmp';

# PCR reagent salt concentration (Molar)
# Following parameters are the default values of Primer3-Web v4.0.0
my $PRIMER_SALT_MONOVALENT = 0.05;
my $PRIMER_SALT_DIVALENT   = 0.0015;
my $PRIMER_DNTP_CONC       = 0.0006;
my $PRIMER_DNA_CONC        = 5e-8;

my $PRIMER_MAX_TM          = $opt{tm_max};
my $PRIMER_OPT_TM          = $opt{tm_opt};
my $PRIMER_MIN_TM          = $opt{tm_min};
my $PRIMER_MAX_SIZE        = $opt{len_max};
my $PRIMER_OPT_SIZE        = $opt{len_opt};
my $PRIMER_MIN_SIZE        = $opt{len_min};

$PRIMER_MAX_TM           ||= 63;
$PRIMER_OPT_TM           ||= 60;
$PRIMER_MIN_TM           ||= 57;
$PRIMER_MAX_SIZE         ||= 27;
$PRIMER_OPT_SIZE         ||= 20;
$PRIMER_MIN_SIZE         ||= 18;

# Path of the libraries required by Primer3
my $PRIMER_THERMODYNAMIC_PARAMETERS_PATH = "$edge_home/lib/primer3_config/";
#my $PRIMER_THERMODYNAMIC_PARAMETERS_PATH = "/users/218817/scratch/opt/apps/bin/primer3_config/";

# Filter out the primer pairs if both primers have alignments < $BG_MISMATCH_CUTOFF mismatches
# 0: keep all primer pairs
my $BG_MISMATCH_CUTOFF = defined $opt{cutoff} ? $opt{cutoff} : 2;

# Background Tm
my $BG_TM_TOLERATE = defined $opt{bg_tm_diff} ? $opt{bg_tm_diff} : 5;

# Display top # qualified primer pairs. 
# Set DISPLAY_TOP_PRIMER to 0 to display all qualified primer pairs.
my $DISPLAY_TOP_PRIMER = defined $opt{top} ? $opt{top} : 5;

# Background database (NCBI blast/BWA indexed db required)
my $BLASTN_DB  = "$edge_home/database/blast_nucl_db/NCBI-Bacteria-Virus.fna";
my $BWA_DB     = "$edge_home/database/bwa_index/NCBI-Bacteria-Virus.fna";
my $BWA_DB_MAP = "$edge_home/database/bwa_index/id_mapping.txt";

my $chunkSize=300000;

my $id_mapping_ref = &readIdMapping( $BWA_DB_MAP );

Tm_calculate->Tm_parameters(
    -oligo => $PRIMER_DNA_CONC,
    -salt  => $PRIMER_SALT_MONOVALENT,
    -mg    => $PRIMER_SALT_DIVALENT,
    -bound => 0.5,
    -dntp  => $PRIMER_DNTP_CONC,
    -sc    => 2 #santalucia correction
);

&pick_primer( $opt{input}, $opt{threads}, "bwa" );

sub pick_primer {
    ## File input
	my ($fasta, $numCPU, $mode) = @_;
    $numCPU ||= 2;

	&Usage() if ( !defined $fasta || !-r $fasta);
    
    # parsing fasta file
    print STDERR " Parsing input file...";
    my ($fasta_r, $contigCnt, $chunkCnt)= &readFastaSeq($fasta);
    print STDERR "$contigCnt contig(s) found. Splitted to $chunkCnt chunks.\n";

    # Identifying shared regions in background sequences
	print STDERR " Identifying shared regions in background sequences...";
	my ($exclude_region, $etime) = &exclude_region($fasta_r, $numCPU);
	print STDERR "done. [$etime]\n";


	# design primer for each contigs
	print STDERR " Picking primer pairs by PRIMER3...";
	my ($primer, $tol_num, $time) = &design_primers($fasta_r, $exclude_region, $numCPU);
	print STDERR "done. $tol_num primer pairs picked. [$time]\n";

	if( $tol_num == 0 ){
		print STDERR "\nNo primer pair is picked.\n";
		exit 0;
	}

    # Checking primer pairs 
	print STDERR " Checking similarity with background sequences...";
	($primer, $tol_num, $time) = &check_primers($primer, $BG_MISMATCH_CUTOFF, $mode, $numCPU);
	print STDERR "done. $tol_num primer pairs passed. [$time]\n";

	if( $tol_num == 0 ){
		print STDERR "\nNo primer pair is qualified.\n";
		exit 0;
	}

    # checking background Tm
	print STDERR " Checking background Tm...";

    unless( $opt{no_tm_chk} ){
	    ($primer, $tol_num, $time) = &check_background_tm($primer, $BG_TM_TOLERATE, $mode, $numCPU);
	    print STDERR "done. $tol_num primer pairs passed after heuristic searching. [$time]\n";

	    if( $tol_num == 0 ){
	    	print STDERR "\nNo primer pair is qualified.\n";
	    	exit 0;
	    }
    }
    else{
        print STDERR "skipped.\n";
    }

	# sort and print suggested primer pairs
	my $count=0;
    my ($p,$ori_fasta_r);

    foreach my $id ( sort
            {
		    	$primer->{$a}->{MIN_BG_TM} <=> $primer->{$b}->{MIN_BG_TM} ||
		    	$primer->{$b}->{TOL_MIN_BG_MISMATCH} <=> $primer->{$a}->{TOL_MIN_BG_MISMATCH} ||
		    	$a cmp $b
		    }  keys %$primer )
	{
        my ($contig,$num) = $id =~ /^(.*)-(\d+)$/;
        $p->{$contig}->{$num}->{COUNT} = $count++;
        $p->{$contig}->{$num}->{ID} = $id;
    }

    foreach my $id ( sort {$fasta_r->{$a}->{order} <=> $fasta_r->{$b}->{order}} keys %$fasta_r ) {
        my $pid = $id;
        $pid =~ s/%part%\d+$//;
        $ori_fasta_r->{$pid}->{len} = 0 unless defined $ori_fasta_r->{$pid}->{len};
        $ori_fasta_r->{$pid}->{len} += $fasta_r->{$id}->{len};
        $ori_fasta_r->{$pid}->{order} = $fasta_r->{$id}->{order};
    }

    $opt{gff3} = "/dev/null" unless defined $opt{gff3};
    open GFF3OUT, ">$opt{gff3}" or die "ERROR: Can't create $opt{gff3} file: $!\n";

    foreach my $contig ( sort {$ori_fasta_r->{$a}->{order} <=> $ori_fasta_r->{$b}->{order}} keys %$ori_fasta_r ){
        my $cnt=1;

        my $len = $ori_fasta_r->{$contig}->{len};
        print "\n[$contig]";

        if( defined $p->{$contig} ){
            print "\n#Length: $len bp\n" if defined $DEBUG && $DEBUG;
            #print "\nPrimer Name\tForward Primer\tForward Tm\tReverse Primer\tReverse Tm\tProduct Length\tUnmatches\tBackground [Tm]\n";
            print "\n\nPrimer Name\tForward Primer\tForward Tm\tReverse Primer\tReverse Tm\tProduct Length\tBackground distance\t[Tm] Background\tAmplicon Location\n";
        }
        else{
            print " -- no qualified primers --\n";
            print "#Length: $len bp\n" if defined $DEBUG && $DEBUG;
        }

        foreach my $num ( sort {$p->{$contig}->{$a}->{COUNT} <=> $p->{$contig}->{$b}->{COUNT}} keys %{$p->{$contig}} ){
            my $id = "$contig-$num";
    		my $out_text = "F".($primer->{$id}->{LEFT}->{MIN_BG_MISMATCH}  == 100 ? "-" : $primer->{$id}->{LEFT}->{MIN_BG_MISMATCH} );
    		$out_text .= "R".($primer->{$id}->{RIGHT}->{MIN_BG_MISMATCH} == 100 ? "-" : $primer->{$id}->{RIGHT}->{MIN_BG_MISMATCH} );
      
            my $out_bg_text = "No background found";

            $out_bg_text = sprintf "[%.2f C] %s", $primer->{$id}->{MIN_BG_TM}, $primer->{$id}->{BG_ORG} if $primer->{$id}->{BG_ORG} ne "NA";

            printf "%s\t%s\t%.1f C\t%s\t%.1f C\t%s bp\t%s\t%s\t%s\n",
    						"$contig-$cnt",
    	                    $primer->{$id}->{LEFT}->{SEQUENCE},
    	                    $primer->{$id}->{LEFT}->{TM},
    	                    $primer->{$id}->{RIGHT}->{SEQUENCE},
    	                    $primer->{$id}->{RIGHT}->{TM},
    	                    $primer->{$id}->{PRODUCT_SIZE},
                            $out_text,
                            $out_bg_text,
                            "$primer->{$id}->{LEFT}->{START}..$primer->{$id}->{RIGHT}->{START}"
                            ;
            
            printf GFF3OUT "$contig\tEDGE\tpcr_amplicon\t%d\t%d\t.\t+\t.\tID=$contig-$cnt;product_length=%s;background=%s\n",
                            $primer->{$id}->{LEFT}->{START},
                            $primer->{$id}->{RIGHT}->{START},
                            $primer->{$id}->{PRODUCT_SIZE},
                            $out_bg_text;

            printf GFF3OUT "$contig\tEDGE\tpcr_primer\t%d\t%d\t.\t+\t.\tID=$contig-$cnt-LEFT;parent=$contig-$cnt;Tm=%.2f;seq=%s;length=%dbp\n",
                            $primer->{$id}->{LEFT}->{START},
                            $primer->{$id}->{LEFT}->{START}+$primer->{$id}->{LEFT}->{LEN}-1,
                            $primer->{$id}->{LEFT}->{TM},
                            $primer->{$id}->{LEFT}->{SEQUENCE},
                            $primer->{$id}->{LEFT}->{LEN};

            printf GFF3OUT "$contig\tEDGE\tpcr_primer\t%d\t%d\t.\t-\t.\tID=$contig-$cnt-RIGHT;parent=$contig-$cnt;Tm=%.2f;seq=%s;length=%dbp\n",
                            $primer->{$id}->{RIGHT}->{START}-$primer->{$id}->{LEFT}->{LEN}+1,
                            $primer->{$id}->{RIGHT}->{START},
                            $primer->{$id}->{RIGHT}->{TM},
                            $primer->{$id}->{RIGHT}->{SEQUENCE},
                            $primer->{$id}->{RIGHT}->{LEN};

    		last if ( $cnt++ == $DISPLAY_TOP_PRIMER );
        }
    }
    close GFF3OUT;
}

sub exclude_region {
	my ($fasta_r, $numCPU) = @_;
	my $exclude_region;
	my $now = time;

    # write primer sets to a fasta file
    open CHUNK, ">$tmp/FASTA_CHUNK_$$.fa" or die "failed to write primers: $!n";
    foreach my $contig ( keys %$fasta_r ){
		print CHUNK ">$contig\n".$fasta_r->{$contig}->{seq}."\n";
	}
	close CHUNK;

    my $cmd = "blastn -db $BLASTN_DB -perc_identity 100 -word_size 30 -ungapped -penalty -1000 -query $tmp/FASTA_CHUNK_$$.fa -num_threads $numCPU -outfmt '6 qseqid qstart length' -out $tmp/PRIMER_EXR_$$.txt";
    print STDERR "\n\t[COMMAND] $cmd\n" if defined $DEBUG && $DEBUG;
    
    &execCmd($cmd, 1);
    
	# run blast & parse result
    open BLAST_OUT, "cat $tmp/PRIMER_EXR_$$.txt | awk -F\\\\t '{if(\$3 >= 30) print \$0}' | sort -u -k1,1 -k2n,2 |" or die "Fail to running blastn: $!\n";

	my($cur_contig, $start, $length) = ("",0,0);

	while(<BLAST_OUT>){
		chomp;
		my @temp = split /\t/, $_;

		($cur_contig, $start, $length) = @temp if( $temp[0] ne $cur_contig );

		if( $start+$length < $temp[1] ){
			$exclude_region->{$cur_contig}->{$start+10}=$length-10;
			($start, $length) = @temp[1..2];
		}
		else{
			$length = $temp[1]+$temp[2]-$start+1;
		}
	}
	close BLAST_OUT;
	
    my $time = &timeInterval($now);

	return ($exclude_region, $time);
}

sub design_primers {
	my ($fasta_r, $exclude_r, $numCPU) = @_;
    my $totalNum=0;
	my $now = time;
    
    # primer3 configuration
    foreach my $id ( keys %$fasta_r ) {
        my $seq = $fasta_r->{$id}->{seq};
        $seq =~ s/[^ATCG]/N/gi;

        $totalNum++;
    
        # preparing the excluded regions
		my @ex;
        my $cnt=1;
		foreach my $start ( sort { $exclude_r->{$id}->{$b} <=> $exclude_r->{$id}->{$a} } keys %{$exclude_r->{$id}} ) {
            # allow max 200 elements in SEQUENCE_EXCLUDED_REGION
            last if $cnt++ == 200;
            push @ex, "$start,".$exclude_r->{$id}->{$start}; 
		}
		my $ex_txt = "";
		$ex_txt = join " ", @ex;
		
        my $primer_num_return = int($fasta_r->{$id}->{len}/150);
        $primer_num_return = $primer_num_return < 50 ? 50 : $primer_num_return;

        my $monov_mm = $PRIMER_SALT_MONOVALENT*(10**3);
        my $div_mm   = $PRIMER_SALT_DIVALENT*(10**3);
        my $dntp_mm  = $PRIMER_DNTP_CONC*(10**3);
        my $dna_nm   = $PRIMER_DNA_CONC*(10**9);

        if( defined $DEBUG && $DEBUG ){
            my $cid = $id;
            $cid =~ s/%part%/ - Chunk /;
            print STDERR "\t[PRIMER3] $cid - ".length($seq)." bp\n";
            print STDERR "\t[PRIMER3] $cid - PRIMER_NUM_RETURN: $primer_num_return\n";
            print STDERR "\t[PRIMER3] $cid - EXCLUDE_REGION: $ex_txt\n";
        } 

        my $p3_config = "SEQUENCE_ID=$id
SEQUENCE_TEMPLATE=$seq
SEQUENCE_EXCLUDED_REGION=$ex_txt
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$PRIMER_THERMODYNAMIC_PARAMETERS_PATH
PRIMER_TASK=pick_detection_primers
PRIMER_PAIR_MAX_DIFF_TM=5.0
PRIMER_PRODUCT_SIZE_RANGE=200-500 501-1000
PRIMER_NUM_RETURN=$primer_num_return
PRIMER_SALT_MONOVALENT=$monov_mm
PRIMER_SALT_DIVALENT=$div_mm
PRIMER_DNTP_CONC=$dntp_mm
PRIMER_DNA_CONC=$dna_nm
PRIMER_SALT_CORRECTIONS=1
PRIMER_TM_FORMULA=1
PRIMER_MIN_GC=30.0
PRIMER_MAX_GC=70.0
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MAX_POLY_X=4
PRIMER_MAX_TM=$PRIMER_MAX_TM
PRIMER_OPT_TM=$PRIMER_OPT_TM
PRIMER_MIN_TM=$PRIMER_MIN_TM
PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE
PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE
PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE
=
";
 
        # save primer config to p3_config file
        open (my $configFH, ">$tmp/p3-$id-$$.config" ) or die "write p3_config fail $!";
        print $configFH $p3_config;
        close $configFH;
   }

    # prepare command
    my $cmd = "cd $tmp; parallel -j $numCPU primer3_core {} ::: p3-*-$$.config > $tmp/PRIMER_P3_$$.out";
    print STDERR "\t[COMMAND] $cmd\n" if defined $DEBUG && $DEBUG;
    &execCmd($cmd,1);

    # run primer3
    open PRIMER3, "$tmp/PRIMER_P3_$$.out" or die "ERROR: can't read Primer3 result.\n";
	
    # parse primer3 result
    my $primer;
    my $id;
    my $p=0;
    my $numWithPrimer3=0;
    while(<PRIMER3>){
        chomp;
        if( /^SEQUENCE_ID=(.*)$/ ){
            $id = $1;
            ($id, $p) = $id =~ /^(.*)%part%(\d+)$/ if ( $id =~ /%part%/ );
        }
        elsif( /^PRIMER_(LEFT|RIGHT)_(\d+)=(\d+),(\d+)$/ ){
            my $pos = $chunkSize*$p;
            $primer->{"$id-$p$2"}->{$1}->{START} = $pos+$3;
            $primer->{"$id-$p$2"}->{$1}->{LEN}   = $4;
        }
        elsif( /^PRIMER_(LEFT|RIGHT)_(\d+)_(TM|SEQUENCE)=([\w\d\.]+)$/ ){
            $primer->{"$id-$p$2"}->{$1}->{$3} = $4;
	    	#initiate minium background mismatch value 100 (do not hit bg)
            $primer->{"$id-$p$2"}->{$1}->{MIN_BG_MISMATCH} = 100;
        }
        elsif( /^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(\d+)/ ){
            $primer->{"$id-$p$1"}->{PRODUCT_SIZE} = $2;
	    	$numWithPrimer3++;
        }
        elsif( /^PRIMER_ERROR=(.*)/ ){
            print STDERR "[PRIMER3 ERROR] $id: $1\n" if defined $DEBUG && $DEBUG;
        }
        elsif( /^=$/ ){
            ($id, $p) = ("",0);
        }
    }
    close PRIMER3;

    my $time = &timeInterval($now);

	return ($primer, $numWithPrimer3, $time);
}

sub check_primers {
	my ($primer,$cutoff,$mode,$numCPU) = @_;
	my $now = time;
	my $qualifiedPairCount=0;
    my %ctg_primer_cnt;

	# write primer sets to a fasta file
	open PRI, ">$tmp/PRIMER_$$.fa" or die "failed to write primers: $!n";
	foreach my $id ( keys %$primer ){
		printf PRI ">$id-LEFT\n".$primer->{$id}->{LEFT}->{SEQUENCE}."\n";
		printf PRI ">$id-RIGHT\n".$primer->{$id}->{RIGHT}->{SEQUENCE}."\n";
	}
	close PRI;

	# run bwa & parse result
    my $cmd = "bwa mem -k5 -A1 -B1 -O3 -E1 -T5 -t$numCPU $BWA_DB $tmp/PRIMER_$$.fa > $tmp/PRIMER_CHK_BWA_$$.sam 2>/dev/null";
    &execCmd($cmd,1);
    print STDERR "\t[COMMAND] $cmd\n" if defined $DEBUG && $DEBUG;

    $cmd = "samtools view -F 4 -S $tmp/PRIMER_CHK_BWA_$$.sam > $tmp/PRIMER_CHK_BWA_$$.mapped.sam 2>/dev/null";
    &execCmd($cmd,1);
    print STDERR "\t[COMMAND] $cmd\n" if defined $DEBUG && $DEBUG;
    
    open BWA_OUT, "$tmp/PRIMER_CHK_BWA_$$.mapped.sam" or die "failed to read BWA output: $!\n";
	while(<BWA_OUT>){
		chomp;
		my @temp = split /\t/, $_;
		my ($id,$fr) = $temp[0]  =~ /^(.*)-(LEFT|RIGHT)$/;
        
        my ($clip5, $clip3) = (0,0);
        $clip5 = $1 if $temp[5]  =~ /^(\d+)S/;
        $clip3 = $1 if $temp[5]  =~ /(\d+)S$/;
        my ($mm) = $temp[11] =~ /NM:i:(\d+)/;
        $mm = $mm+$clip5+$clip3;

        $primer->{$id}->{$fr}->{MIN_BG_MISMATCH} = $mm;
	}
	close BWA_OUT;

	# remove non-qualified pairs - background mismatch excess the cutoff on both primers
	foreach my $id ( keys %$primer ) {
        my $tol = $primer->{$id}->{LEFT}->{MIN_BG_MISMATCH} + $primer->{$id}->{RIGHT}->{MIN_BG_MISMATCH};
		
        if( $tol < $cutoff ){
			delete $primer->{$id};
		}
        else{
            $primer->{$id}->{TOL_MIN_BG_MISMATCH} = $tol;
			$qualifiedPairCount++;
        }
	}

    # [HEURISTIC] only keep triple display number of primers for improving speed
    if( $DISPLAY_TOP_PRIMER > 0 ){
        foreach my $id ( sort { 
                $primer->{$b}->{TOL_MIN_BG_MISMATCH} <=> $primer->{$a}->{TOL_MIN_BG_MISMATCH} || 
                $a cmp $b 
            }  keys %$primer )
    	{
            my ($contig,$num) = $id =~ /^(.*)-(\d+)$/;
            $ctg_primer_cnt{$contig}=0 unless defined $ctg_primer_cnt{$contig};
            $ctg_primer_cnt{$contig}++;
            
            if( $ctg_primer_cnt{$contig} > $DISPLAY_TOP_PRIMER*3 ){
                delete $primer->{$id};
            }
        }
    }

    my $time = &timeInterval($now);

    return ($primer, $qualifiedPairCount, $time);
}

sub check_background_tm {
    my ($primer,$cutoff,$mode,$numCPU) = @_;
    my $now = time;
    my $qualifiedPairCount=0;

    open BWA_OUT, "$tmp/PRIMER_CHK_BWA_$$.mapped.sam" or die "failed to read BWA output: $!\n";

    while(<BWA_OUT>){
        chomp;
	my @temp = split /\t/, $_;
	my ($id,$fr) = $temp[0]  =~ /^(.*)-(LEFT|RIGHT)$/;
 
        #check current available primers
        next unless defined $primer->{$id};
        
        if( $primer->{$id}->{$fr}->{MIN_BG_MISMATCH} < 100 ){
            my $org = $id_mapping_ref->{$temp[2]};
            $org =~ s/,.*$// if defined $org;
            $primer->{$id}->{$fr}->{BG_ORG} = defined $org ? $org : $temp[2];

            my $tm = $primer->{$id}->{$fr}->{TM};
            $tm = &samToTm(@temp) if( $primer->{$id}->{$fr}->{MIN_BG_MISMATCH} > 0 );
            $primer->{$id}->{$fr}->{MIN_BG_TM} = $tm;
        }
        else{ # primer with no background hit
            $primer->{$id}->{$fr}->{BG_ORG} = "NA";
            $primer->{$id}->{$fr}->{MIN_BG_TM} = 0;
        }
    }
    close BWA_OUT;

    # remove non-qualified pairs - background mismatch excess the cutoff on both primers
    foreach my $id ( keys %$primer ) {
        my $tm = $primer->{$id}->{LEFT}->{TM};
        $tm = $primer->{$id}->{RIGHT}->{TM} if $primer->{$id}->{RIGHT}->{TM} < $tm;
        
        my $bgtm;
        $bgtm = $primer->{$id}->{RIGHT}->{MIN_BG_TM};
        $primer->{$id}->{BG_ORG} = $primer->{$id}->{RIGHT}->{BG_ORG};

        if( $primer->{$id}->{LEFT}->{MIN_BG_TM} < $bgtm ){
            $bgtm = $primer->{$id}->{LEFT}->{MIN_BG_TM};
            $primer->{$id}->{BG_ORG} = $primer->{$id}->{LEFT}->{BG_ORG};
        }

        $primer->{$id}->{BG_ORG} = $primer->{$id}->{LEFT}->{BG_ORG};
        
        if( $tm-$bgtm < $cutoff ){
            delete $primer->{$id};
        }
        else{
            $primer->{$id}->{MIN_BG_TM} = $bgtm;
            $qualifiedPairCount++;
        }
    }

    my $time = &timeInterval($now);

    return ($primer, $qualifiedPairCount, $time);
}

sub samToTm {
    my @temp = @_;

    my ($clip5, $clip3, $ins, $del, $len, $mm, $seq, $ref) = (0, 0, 0, 0, 0, "", "");
    $clip5 = $1 if $temp[5]  =~ /^(\d+)S/;
    $clip3 = $1 if $temp[5]  =~ /(\d+)S$/;
    $mm    = $1 if $temp[11] =~ /NM:i:(\d+)/;
    $len   = length($temp[9]);
    $seq   = $temp[9];

    #calculate indels
    $ins += $1 while( $temp[5] =~ /(\d+)I/g );
    $del += $1 while( $temp[5] =~ /(\d+)D/g );

    #get reference seq (note that the start position has been excluded clipping length)
    my $start = ($temp[3] - $clip5) > 0 ? ($temp[3] - $clip5) : 1; #the soft-clip starts prior the reference
    $ref = &getRefSeq($BWA_DB, $temp[2], $start, ($len-$ins+$del) );
    $ref = ('A'x(1-$temp[3]+$clip5)).$ref if ($temp[3] - $clip5) < 1; #faking the prior sequence
    
    #replace every bases other than ATCG to A
    $ref =~ s/[^ATCG]/A/g;

    my ($seq_aln, $ref_aln, $aln) = &cigar2alignment($temp[5], $seq, $ref);

    #dangling sequence of primer sequence
    my ($qs5, $qm, $qs3) = ("","","");
    
    if( $clip5 ){
        my $cnt = $clip5-1;
        ($qs5) = $seq_aln =~ /^\w{$cnt}(\w)/;
        $qs5 = "-$qs5";
    }
    if( $clip3 ){
        my $cnt = $clip3-1;
        ($qs3) = $seq_aln =~ /(\w)\w{$cnt}$/;
        $qs3 = "$qs3-";
    }

    ($qm) = $seq_aln =~ /^\w{$clip5}(.+)\w{$clip3}$/;

    $seq_aln = "$qs5$qm$qs3";
    
    #dangling sequence of reference sequence
    ($qs5, $qm, $qs3) = ("","","");
    
    if( $clip5 ){
        my $cnt = $clip5-1;
        ($qs5) = $ref_aln =~ /^\w{$cnt}(\w)/;
        $qs5 = "-$qs5";
    }
    if( $clip3 ){
        my $cnt = $clip3-1;
        ($qs3) = $ref_aln =~ /(\w)\w{$cnt}$/;
        $qs3 = "$qs3-";
    }

    ($qm) = $ref_aln =~ /^\w{$clip5}(.+)\w{$clip3}$/;

    $ref_aln = "$qs5$qm$qs3";

    #calculating Tm
    my $tm = Tm($seq_aln, $ref_aln);

    #return Tm, total number of unmatched bases
    return $tm;
}

sub cigar2alignment
{
    my ($cigar, $seq, $ref) = @_;
    my $cstr = "";
    while( $cigar =~ /(\d+)(\w)/g ){
        $cstr .= "$2"x$1;
    }

    my $seq_m = "";
    my $ref_m = "";
    my $aln   = "";

    my @c = split //, $cstr;
    my @s = split //, $seq;
    my @r = split //, $ref;

    for( my $i=0; $i<=$#c; $i++){
        if( $c[$i] eq "M" ){
            my $curs = shift @s;
            my $curr = shift @r;
            $seq_m .= $curs;
            $ref_m .= $curr;
            $aln   .= "|" if $curs eq $curr;;
            $aln   .= " " if $curs ne $curr;;
        }
        elsif( $c[$i] eq "S" ){
            $seq_m .= shift @s;
            $ref_m .= shift @r;
            $aln   .= ".";
        }
        elsif( $c[$i] eq "I" ){
            $seq_m .= shift @s;
            $ref_m .= "-";
            $aln   .= " ";
        }
        elsif( $c[$i] eq "D" ){
            $seq_m .= "-";
            $ref_m .= shift @r;
            $aln   .= " ";
        }
    }

    return ($seq_m, $ref_m, $aln);
}

sub getRefSeq
{
    my ($file, $refid, $start, $len) = @_;    
    my $end = $start+$len-1;
    my $cmd = "samtools faidx $file '$refid:$start-$end'";
    my $output = `$cmd`;
    my ($ref) = $output =~ />\S+\n(\w+)\n/;

    return $ref;
}

sub readIdMapping
{
    my ($file) = @_;
    my $ref;
    open MAP, $file or die "Can't open file: $!.\n";
    while(<MAP>){
        chomp;
        next if /^$/;
        my ($idx,$val) = $_ =~ /^(\S+) (.*)/;
        $ref->{$idx}=$val;
    }
    close MAP;

    return $ref;
}

sub readFastaSeq
{
    my $seqFile=shift;
    my %hash;
    my $count=0;
    my $contigCnt=0;
    my $chunkCnt=0;

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
        $contigCnt++;

        if( $len < $chunkSize ){
            $hash{$id}->{seq}=$seq;
            $hash{$id}->{len}=$len;
            $hash{$id}->{order}=$count++;
            $chunkCnt++;
        }
        if( $len >= $chunkSize ){
            my $cnt=1;
            for( my $start=0; $start<$len; $start+=$chunkSize ){
                my $l =  ($start+$chunkSize) > $len ? ($len-$start) : $chunkSize;
                my $p = "$id%part%$cnt";
                $hash{$p}->{seq} = substr($seq, $start, $l);
                $hash{$p}->{len} = $l;
                $hash{$p}->{order} = $count++;
                $cnt++;
                $chunkCnt++;
            }
        }
    }
    close $fh;
    $/="\n";
	
	return (\%hash, $contigCnt, $chunkCnt);
}

sub timeInterval{
    my $now = shift;
	$now = time - $now;
    return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}

sub execCmd {
    my $command = $_[0];
    my $fatal = $_[1];

    my $rv = system($command);
    if($rv != 0) {
        die "\nFAILED: $command - $!\n" if $fatal;
    }   
}

sub Usage
{
    print <<USAGE;
USAGE:
    perl $0 [OPTIONS] --input <FILE>

    --input | -i <FILE> : input FASTA file

OPTIONS:
    --threads    | -t INT  : number of threads (default: 2)
    --gff3       | -g FILE : output GFF3 file
    --cutoff     | -c INT  : total unmatches (clipping + mismatch) cutoff (default: 2)
    --no_tm_chk  | -n      : ignore background Tm checking
    --bg_tm_diff | -f INT  : reject primer having Tm < tm_diff difference with background Tm (default: 5)
    --top        | -t INT  : display # top results

USAGE
exit;
}

