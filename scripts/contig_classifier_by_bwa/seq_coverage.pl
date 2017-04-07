#!/usr/bin/perl
use Getopt::Long;
use FindBin qw($Bin);
use lib $Bin;
use gi2lineage;
use strict;

my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'type|t=s',
    'orig_seq|p=s',
    'help|?') || &usage();

&usage() if !defined $opt{input} || !-e $opt{input};

print STDERR "Input file: $opt{input}\n";

#preload taxonomy data
#loadTaxonomy("preload");
loadTaxonomy();
print STDERR "Done loading taxanomy data.\n";

my $input = $opt{input};
my $type  = $opt{type};
my $file  = $opt{orig_seq};

$type ||= "guess";

#############################################################################################################
## Support 3 alignment format:
#############################################################################################################
#
#  # LAST tab output    # BLAST M8       SAM Field  Description
#  0  score             0  qseqid        0   QNAME  Query (pair) NAME
#  1  name1 (ref)       1  sseqid        1   FLAG   bitwise FLAG
#  2  start1            2  pident        2   RNAME  Reference sequence NAME
#  3  alnSize1          3  length        3   POS    1-based leftmost POSition/coordinate of clipped sequence
#  4  strand1           4  mismatch      4   MAPQ   MAPping Quality (Phred-scaled)
#  5  seqSize1          5  gapopen       5   CIAGR  extended CIGAR string
#  6  name2 (query)     6  qstart        6   MRNM   Mate Reference sequence NaMe (‘=’ if same as RNAME)
#  7  start2            7  qend          7   MPOS   1-based Mate POSistion
#  8  alnSize2          8  sstart        8   ISIZE  Inferred insert SIZE
#  9  strand2           9  send          9   SEQ    query SEQuence on the same strand as the reference
#  10 seqSize2          10 evalue        10  QUAL   query QUALity (ASCII-33 gives the Phred base quality)
#  11 blocks            11 bitscore      11  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
#                                        **** OPT fields
#                                        11  AS:i:
#                                        12  XS:i:
#                                        13  XF:i:
#                                        14  XE:i:
#                                        15  NM:i:
#
#############################################################################################################

my $seq;
my $cov;
my $cnt=0;
my $length;

$length = &readFastaSeq($file) if -e $file;
print STDERR "Done loading original ".scalar(keys %$length)." sequences.\n" if -e $file;

open INPUT, $input or die "ERROR: can't open input file: $!\n";
while(<INPUT>)
{
    next if /^#/;
    chomp;
    my @temp = split /\t/, $_;

    # guess input format
    if( $type eq "guess" ){
        $type = "last"  if scalar @temp == 12 && $temp[0] =~ /^\d+$/;
        $type = "blast" if scalar @temp >= 12 && $temp[8] =~ /^\d+$/;
        $type = "sam" if scalar @temp > 13;
        die "ERROR: Can not recognize input format!\n" if $type eq "guess";
        print STDERR "Guess input format: $type.\n";
    }

    my ($sid, $qid, $qlen, $qstart, $qend, $dist) = ("","",0,0,0,1);
    
    if( $type =~ /^last/i ){
        my $len = defined $length->{$temp[6]} ? $length->{$temp[6]} : $temp[10];
        #LAST use 0-based position
        ($sid, $qid, $qlen, $qstart, $qend) = ($temp[1], $temp[6], $len, $temp[7]+1, ($temp[7]+$temp[8]) ); 
    }
    elsif( $type =~ /^blast/i ){
        my $len = defined $length->{$temp[0]} ? $length->{$temp[0]} : $temp[12];
        die "No query sequence found." unless defined $len;
        ($sid, $qid, $qlen, $qstart, $qend, $dist) = ($temp[1], $temp[0], $len, $temp[6], $temp[7], ($temp[4]+$temp[5]) );
    }
    elsif( $type =~ /^sam/i ){
        my ($nm) = $_ =~ /NM:i:(\d+)/;
        $length->{$temp[0]} ||= length($temp[9]);
        my $len = $length->{$temp[0]};
        my $clip5 = $1 if $temp[5] =~ /^(\d+)[SH]/;
        my $clip3 = $1 if $temp[5] =~ /(\d+)[SH]$/;
        
        if( $temp[5] =~ /r/ ){
            my $temp = $clip5;
            $clip5 = $clip3;
            $clip3 = $temp;
        }
        
        ($sid, $qid, $qlen, $qstart, $qend, $dist) = ($temp[2], $temp[0], $len, ($clip5+1), ($len-$clip3), $nm);
    }
    else{
        die "ERROR: unknown input format.\n";
    }

    #my ($gi) = $sid =~ /gi\|(\d+)/;
    my $acc = getAccFromSeqID($sid);
    my $taxid = acc2taxID($acc);
    $length->{$qid}=$qlen;
    $seq->{$qid}->{$taxid}->{"$qstart..$qend"}=$dist;
	$cov->{$qid}->{$cnt}->{$taxid}="$qstart..$qend";

    $cnt++;
}
close INPUT;
print STDERR "Done loading $cnt $type hits.\n";

print "##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tACC_COV_RGN\tACC_COV_LEN\n";
# 1  SEQ = query sequence name
# 2  RANK = rank name
# 3  ORGANISM = taxonomy name
# 4  TAX_ID = taxonomy id
# 5  PARENT = parent taxonomy name
# 6  LENGTH = query sequence name
# 7  NUM_HIT = number of hits
# 8  TOL_HIT_LEN = total bases of hits
# 9  TOL_MISM = total bases of mismatches
# 10 AVG_IDT = average hit identity
# 11 LINEAR_LEN = linear length of hits
# 12 RANK_LINEAR_LEN = total linear length of hits in the certain rank
# 13 COV = LINEAR_LEN/LENGTH
# 14 SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN

foreach my $pname ( sort keys %$seq )
{
    my $p;
    my $r;
    my $upper_level = "root";

    foreach my $rank (("superkingdom","phylum","class","order","family","genus","species","strain"))
    {
        foreach my $taxid ( keys %{$seq->{$pname}} )
        {
            my $name = taxid2rank($taxid, $rank);
            
            #upper taxa
            my $upname = taxid2rank($taxid, $upper_level);
            $upname = "NA" unless $upname;
            $name = "$upname $rank" unless $name;
            
            $p->{$rank}->{$name}->{UP_RANK} = $upname;
            $p->{$rank}->{$name}->{TAXO_ID} = taxid2rank_taxid($taxid,$rank);
            $p->{$rank}->{$name}->{TAXO_ID} = $taxid if $rank eq "strain";

            my $pcov;

            if( defined $p->{$rank}->{$name}->{LINEAR_LEN} ){
                $pcov = $p->{$rank}->{$name}->{LINEAR_LEN};
            }
            else{
                $pcov = "0"x$length->{$pname};
            }

            foreach my $region ( keys %{$seq->{$pname}->{$taxid}} ){
                my ($qs, $qe) = $region =~ /^(\d+)\.\.(\d+)$/;
                my $end = length($pcov);
                my $nm = $seq->{$pname}->{$taxid}->{$region};
                
                #update linear length
                my $str = "0"x($qs-1) . "1"x($qe-$qs+1) . "0"x($end-$qe);
                $pcov = $pcov | $str;
                #total mapped
                $p->{$rank}->{$name}->{TOL_HIT_LEN} ||= 0;
                $p->{$rank}->{$name}->{TOL_HIT_LEN} += $qe-$qs+1;
                #number of hits
                $p->{$rank}->{$name}->{NUM_HIT} ||= 0;
                $p->{$rank}->{$name}->{NUM_HIT}++;

                #distance
                $p->{$rank}->{$name}->{TOL_MISM} ||= 0;
                $p->{$rank}->{$name}->{TOL_MISM} += $nm;
            }
            $p->{$rank}->{$name}->{AVG_IDT} = ($p->{$rank}->{$name}->{TOL_HIT_LEN} - $p->{$rank}->{$name}->{TOL_MISM})/$p->{$rank}->{$name}->{TOL_HIT_LEN};
            $p->{$rank}->{$name}->{LINEAR_LEN} = $pcov;
        }

        # 
        foreach my $name ( keys %{$p->{$rank}} ){
            $p->{$rank}->{$name}->{LINEAR_LEN} =~ s/0//g;
            my $sum = length($p->{$rank}->{$name}->{LINEAR_LEN});

            $p->{$rank}->{$name}->{LINEAR_LEN} = $sum;
            $r->{$rank}->{TOL_LINEAR_LEN} ||= 0;
            $r->{$rank}->{TOL_LINEAR_LEN} += $sum;
        }
    	
		#accumulated coverage
    	my $acc_cov;
    	$acc_cov = "0"x$length->{$pname};
    	my $map;
		my $map->{48}="unclassified";
		my $mid_ascii = 49;

    	foreach my $cnt ( sort {$a<=>$b} keys %{$cov->{$pname}} )
    	{
    		foreach my $taxid ( keys %{$cov->{$pname}->{$cnt}} )
    		{
		        my $name = taxid2rank($taxid, $rank);    
    	        #upper taxa
		        my $upname = taxid2rank($taxid, $upper_level);
    	        $upname = "NA" unless $upname;
    		    $name = "$upname $rank" unless $name;

				unless( defined $map->{$name} ){
					$map->{$name} = $mid_ascii;
					$map->{$mid_ascii} = $name;
					$mid_ascii++;
				}

		    	my $region = $cov->{$pname}->{$cnt}->{$taxid};
		    	$acc_cov = &accCov($acc_cov, chr($map->{$name}), $region);
    		}
    	}
		my $csum = &accCovSummary($acc_cov, $map);
		foreach my $name ( keys %$csum )
		{
			$p->{$rank}->{$name}->{ACC_COV_RGN} = join ";", @{$csum->{$name}};
			my $len=0;
			foreach my $rgn ( @{$csum->{$name}} )
			{
				my ($qs,$qe) = $rgn =~ /(\d+)\.\.(\d+)/;
				$len += $qe-$qs+1;
			}
			$p->{$rank}->{$name}->{ACC_COV_LEN} = $len;
		}
    
        $upper_level = $rank;
    }

    foreach my $rank (("superkingdom","phylum","class","order","family","genus","species","strain"))
    {
        foreach my $name ( sort {
                $p->{$rank}->{$b}->{ACC_COV_LEN} <=> $p->{$rank}->{$a}->{ACC_COV_LEN} 
            } keys %{$p->{$rank}} )
        {
			next if $name eq "unclassified";
            my $pref = $p->{$rank}->{$name};
            #sequence rank orig taxid parent #hit tol_mapped_len linear_len tol_linear_len coverage prob
            printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%s\t%d\n",
                $pname,
                $rank,
                $name,
                $pref->{TAXO_ID},
                $pref->{UP_RANK},
                $length->{$pname},
                $pref->{NUM_HIT},
                $pref->{TOL_HIT_LEN},
                $pref->{TOL_MISM},
                $pref->{AVG_IDT},
                $pref->{LINEAR_LEN},
                $r->{$rank}->{TOL_LINEAR_LEN},
                $pref->{LINEAR_LEN}/$length->{$pname},
                $pref->{LINEAR_LEN}/$r->{$rank}->{TOL_LINEAR_LEN},
                $pref->{ACC_COV_RGN},
                $pref->{ACC_COV_LEN}
            ;
        }
    }

	my $pref = $p->{"superkingdom"}->{"unclassified"};
	printf "%s\t%s\t%s\t\t\t%d\t\t\t\t\t%d\t\t%.4f\t\t%s\t%d\n",
	    $pname,
	    "unclassified",
	    "unclassified",
	    $length->{$pname},
	    $pref->{ACC_COV_LEN},
	    $pref->{ACC_COV_LEN}/$length->{$pname},
	    $pref->{ACC_COV_RGN},
	    $pref->{ACC_COV_LEN}
	;
}

sub accCov {
	my ($acc_cov,$id,$region) = @_;
	my ($qs, $qe) = $region =~ /^(\d+)\.\.(\d+)$/;

	while( $acc_cov =~ /0+/g ){
		my ($us, $ue) = ($-[0]+1,$+[0]);
		last if $us > $qe;
		if( $qs>=$us && $qe<=$ue){ #whole overlapping
			my $len = $qe-$qs+1;
			substr $acc_cov, $qs-1, $len, ${id}x$len;
		}
		elsif( $us>=$qs && $qe>=$us && $ue>=$qe ){ #cov overlapping 3" 0s
			my $len = $ue-$qs+1;
			substr $acc_cov, $qs-1, $len, ${id}x$len;
		}
		elsif( $qs>=$us && $qs<=$ue && $qe>$ue ){ #overlapping 5"
			my $len = $qe-$us+1;
			substr $acc_cov, $us-1, $len, ${id}x$len;
		}
	}
	
	return $acc_cov;
}

sub accCovSummary {
	my ($acc_cov,$map) = @_;
	my $c;
	my $csum;
	while( $acc_cov =~ /(.)\1*/g ){
		my ($qs,$qe) = ($-[0]+1,$+[0]);
		my $tax = $map->{ord($1)};
		
		#dealing with an upper limit of '32766' on the MAX value of the regex {MIN,MAX} quantifier.
		my ($prev_end) = $csum->{$tax}[-1] =~ /\.\.(\d+)/;
		if( defined $prev_end && $prev_end+1 == $qs ){
			$csum->{$tax}[-1] =~ s/\.\.$prev_end/\.\.$qe/;
		}
		else{
			push @{$csum->{$tax}}, "$qs..$qe";
		}
	}
	return $csum;
}

sub readFastaSeq
{
    my $seqFile=shift;
    my $hash;
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
        $hash->{$id}=$len;
    }
    close $fh;
    $/="\n";
    
    return $hash;
}

sub usage {
    print <<__END__;

$0 [OPTIONS] --input <FILENAME>

    --input   | -i <STRING>  the sequence of contigs in FASTA format

[OPTIONS]

    --type     | -t <STR>   input format
    --orig_seq | -p <FILE>  original sequence file
    --help/h/?              display this help                   

__END__
exit();
}

