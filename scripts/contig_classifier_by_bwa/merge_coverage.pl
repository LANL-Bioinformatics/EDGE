#!/usr/bin/perl -w
use strict;

my $cov;
my $length;
my $old_contig;
my $input=$ARGV[0];
my $contig=$ARGV[1];

print "##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tNUM_MERGED\tACC_COV_LEN\n";

open (IN,$input) or die "Cannot read $input\n";
while(<IN>){
    next if /^#/;
    chomp;
    #SEQ    RANK    ORGANISM    TAX_ID  PARENT  LENGTH  NUM_HIT TOL_HIT_LEN TOL_MISM    AVG_IDT LINEAR_LEN  RANK_LINEAR_LEN COV SCALED_COV  ACC_COV_RGN  ACC_COV_LEN
    # 0       1        2           3      4       5        6         7         8           9        10              11       12     13           14           15
    my @temp = split /\t/, $_;
    #my ( $contig, $id ) = $temp[0] =~ /^(.*)_(\d+)$/;
    #if( !defined $contig || !defined $id ){
	#$contig = "assembly";
        my $id = $temp[0];
    #}
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{TAXA_ID}         = $temp[3];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{PARENT}          = $temp[4];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{LENGTH}          = $temp[5];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{NUM_HIT}         = $temp[6];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{TOL_HIT_LEN}     = $temp[7];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{TOL_MISM}        = $temp[8];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{LINEAR_LEN}      = $temp[10];
    $cov->{$contig}->{$id}->{$temp[1]}->{$temp[2]}->{ACC_COV_LEN}     = $temp[15];
    $length->{$contig}->{$id}->{LENGTH} = $temp[5];

    $old_contig ||= $contig;

    # merging parts when the classification of a contig is fully loaded
    if( $old_contig ne $contig ){
        covergeMerger($old_contig);
        delete $cov->{$old_contig};
        delete $length->{$contig};
        $old_contig = $contig;
    }
}
close IN;
covergeMerger($old_contig);

# 1   SEQ = query sequence name
# 2   RANK = rank name
# 3   ORGANISM = taxonomy name
# 4   TAX_ID = taxonomy id
# 5   PARENT = parent taxonomy name
# 6   LENGTH = query sequence name
# 7   NUM_HIT = number of hits
# 8   TOL_HIT_LEN = total bases of hits
# 9   TOL_MISM = total bases of mismatches
# 10  AVG_IDT = average hit identity
# 11  LINEAR_LEN = linear length of hits
# 12  RANK_LINEAR_LEN = total linear length of hits in the certain rank
# 13  COV = LINEAR_LEN/LENGTH
# 14  SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN
# 15  NUM_MERGED = number of sequences merged

sub covergeMerger
{
    my $contig = shift;
    my $cds_len = 0;
    my $ctg_cov;
    my $r;
    
    foreach my $id ( keys %{$cov->{$contig}} )
    {
        $cds_len += $length->{$contig}->{$id}->{LENGTH};

        foreach my $rank ( keys %{$cov->{$contig}->{$id}} )
        {
            foreach my $name ( keys %{$cov->{$contig}->{$id}->{$rank}} )
            {
                if( !defined $ctg_cov->{$contig}->{$rank}->{$name} )
                {
                    $ctg_cov->{$contig}->{$rank}->{$name}->{NUM_CDS}     = 0;
                    $ctg_cov->{$contig}->{$rank}->{$name}->{NUM_HIT}     = 0;
                    $ctg_cov->{$contig}->{$rank}->{$name}->{TOL_HIT_LEN} = 0;
                    $ctg_cov->{$contig}->{$rank}->{$name}->{LINEAR_LEN}  = 0;
                    $ctg_cov->{$contig}->{$rank}->{$name}->{TOL_MISM}    = 0;
                }
                
                my $cov_ctg_ref = $ctg_cov->{$contig}->{$rank}->{$name};
                $cov_ctg_ref->{NUM_CDS}++;
                $cov_ctg_ref->{PARENT}         = $cov->{$contig}->{$id}->{$rank}->{$name}->{PARENT};
                $cov_ctg_ref->{TAXA_ID}        = $cov->{$contig}->{$id}->{$rank}->{$name}->{TAXA_ID};
                $cov_ctg_ref->{NUM_HIT}       += $cov->{$contig}->{$id}->{$rank}->{$name}->{NUM_HIT} || 0;
                $cov_ctg_ref->{TOL_HIT_LEN}   += $cov->{$contig}->{$id}->{$rank}->{$name}->{TOL_HIT_LEN} || 0;
                $cov_ctg_ref->{LINEAR_LEN}    += $cov->{$contig}->{$id}->{$rank}->{$name}->{LINEAR_LEN};
                $cov_ctg_ref->{TOL_MISM}      += $cov->{$contig}->{$id}->{$rank}->{$name}->{TOL_MISM} || 0;
                $cov_ctg_ref->{ACC_COV_LEN}   += $cov->{$contig}->{$id}->{$rank}->{$name}->{ACC_COV_LEN};
                $r->{$rank} ||= 0;
                $r->{$rank} += $cov->{$contig}->{$id}->{$rank}->{$name}->{LINEAR_LEN};
            }
        }
    }

    foreach my $rank (("superkingdom","pylum","class","order","family","genus","species","strain"))
    {
        foreach my $name ( sort {
                $ctg_cov->{$contig}->{$rank}->{$b}->{ACC_COV_LEN} <=> $ctg_cov->{$contig}->{$rank}->{$a}->{ACC_COV_LEN}
				#$ctg_cov->{$contig}->{$rank}->{$a}->{TOL_MISM} <=> $ctg_cov->{$contig}->{$rank}->{$b}->{TOL_MISM}
                } keys %{ $ctg_cov->{$contig}->{$rank}} )
        {

            my $cov_ctg_ref = $ctg_cov->{$contig}->{$rank}->{$name};
            printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%d\t%d\n",
                $contig,                                                 # 1   SEQ = query sequence name                                        
                $rank,                                                   # 2   RANK = rank name                                                 
                $name,                                                   # 3   ORGANISM = taxonomy name                                         
                $cov_ctg_ref->{TAXA_ID},                                 # 4   TAX_ID = taxonomy id                                             
                $cov_ctg_ref->{PARENT},                                  # 5   PARENT = parent taxonomy name                                    
                $cds_len,                                                # 6   LENGTH = query sequence name                                     
                $cov_ctg_ref->{NUM_HIT},                                 # 7   NUM_HIT = number of hits                                         
                $cov_ctg_ref->{TOL_HIT_LEN},                             # 8   TOL_HIT_LEN = total bases of hits                                
                $cov_ctg_ref->{TOL_MISM},                                # 9   TOL_MISM = total bases of mismatches                             
                ($cov_ctg_ref->{TOL_HIT_LEN}-$cov_ctg_ref->{TOL_MISM})/$cov_ctg_ref->{TOL_HIT_LEN}, # 10  AVG_IDT = average hit identity                                   
                $cov_ctg_ref->{LINEAR_LEN},                              # 11  LINEAR_LEN = linear length of hits                               
                $r->{$rank},                                             # 12  RANK_LINEAR_LEN = total linear length of hits in the certain rank
                $cov_ctg_ref->{LINEAR_LEN}/$cds_len,                     # 13  COV = LINEAR_LEN/LENGTH
                $cov_ctg_ref->{LINEAR_LEN}/$r->{$rank},                  # 14  SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN
                $cov_ctg_ref->{NUM_CDS},                                 # 15  NUM_MERGED = number of sequences merged
                $cov_ctg_ref->{ACC_COV_LEN}                              # 16  ACC_COV_LEN = number of merged ACC_COV_LEN
            ;                                                            
        }                                                                
    }
    
	my $cov_ctg_ref = $ctg_cov->{$contig}->{"unclassified"}->{"unclassified"};
    printf "%s\t%s\t%s\t\t\t%d\t\t\t\t\t%d\t\t%.4f\t\t\t%d\n",
        $contig,                             
        "unclassified",                      
        "unclassified",                      
        $cds_len,                            
        $cov_ctg_ref->{LINEAR_LEN},          
        $cov_ctg_ref->{LINEAR_LEN}/$cds_len, 
        $cov_ctg_ref->{ACC_COV_LEN}          
    ;                               
}
