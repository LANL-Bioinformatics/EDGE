#!/usr/bin/perl -w
# prokka_gff3_stats.pl is used to generate annotation stats from PROKKA gff3 output
#
# Required third-party tools:
# - R
#
# Po-E Li 20130926

use FindBin qw($Bin);
use Getopt::Long;
use strict;

$ENV{PATH} = "$Bin:$Bin/../bin/:$Bin/../scripts:$ENV{PATH}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'title|t=s',
    'prefix|p=s',
    'outfmt|f=s',
    'help|?') || &usage();

if ( $opt{help} || !defined $opt{input} ) { &usage(); }

$opt{title}  ||= "Contigs";
$opt{prefix} ||= "prokka_gff3_stats";
$opt{outfmt} ||= "PDF";

my $stats;
my $ctg_len;
my $q_feat = "CDS"; #for prokka gff
my $q_attr;
my $ignore_contig = 1;
my $ctg_num=0; #total number of contigs


my $ec_num;
$ec_num->{1}="EC1: Dehydrogenase, oxidase";
$ec_num->{2}="EC2: Transaminase, kinase";
$ec_num->{3}="EC3: Lipase, amylase, peptidase";
$ec_num->{4}="EC4: Decarboxylase";
$ec_num->{5}="EC5: Isomerase, mutase";
$ec_num->{6}="EC6: Synthetase";

open GFF3, "$opt{input}" or die "Can't open $opt{input}: $!\n";
my $ctg;

while(<GFF3>){
    chomp;
    
    if ( /^##sequence-region (\S+) 1 (\d+)$/ ){
        $ctg = $ignore_contig ? $opt{title} : $1;
        $ctg_len->{$ctg} = 0 unless defined $ctg_len->{$ctg};
        $ctg_len->{$ctg} += $2;
        $ctg_num++;
        next;
    }
    
    my @temp = split /\t/, $_;
    next if scalar @temp != 9;

    $temp[0] = $ctg if $ignore_contig;

    #counting features
    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{COUNT} );
    $stats->{$temp[0]}->{$temp[2]}->{COUNT}++;

    # counting rRNA features
    if( $temp[2] eq 'rRNA' ){
        my ($rrna_type) = $temp[8] =~ /product=([^;]+)/;
        &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$rrna_type}->{COUNT} );
        $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$rrna_type}->{COUNT}++;
    }
    
    #counting GENE attributes
    if( $temp[2] eq $q_feat ){
        my @a_temp = split /;/, $temp[8];
        foreach my $attr ( @a_temp ){
            my ($tag, $val) = $attr =~ /^([^=]+)=(.*)$/;

            #skip tags
            next if $tag =~ /^(locus_tag|note)$/;

            $tag = uc($tag);

            &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{COUNT} );
            $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{COUNT}++;

            # inference
            if( $tag =~ /inference/i ){
                my @values = split /,/, $val;
                foreach my $value ( @values ){
                    my ($cate) = $value =~ /:([^:]+):/;
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{TYPE}->{$cate}->{COUNT} );
                    $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{TYPE}->{$cate}->{COUNT}++;
                }
            }

            # EC
            if( $tag =~ /eC_number/i ){
                my @nums = split /,/, $val;
                foreach my $ec ( @nums ){
                    my ($cate) = $ec =~ /^(\d+)\./;
                    $cate = $ec_num->{$cate};
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{TYPE}->{$cate}->{COUNT} );
                    $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{$tag}->{TYPE}->{$cate}->{COUNT}++;
                }
            }

            #hypothetical protein
            if( $tag =~ /product/i && $val =~ /(hypothetical|putative|proteins of unknown function)/i ){
               
                $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{COUNT}++;
                
                my $ref = $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE};

                if( $val =~ /conserved hypothetical protein/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"conserved hypothetical"}->{COUNT} );
                    $ref->{"conserved hypothetical"}->{COUNT}++;
                }
                elsif( $val =~ /conserved.*domain hypothetical protein/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"conserved domain hypothetical"}->{COUNT} );
                    $ref->{"conserved domain hypothetical"}->{COUNT}++;
                }
                elsif( $val =~ /domain protein$/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{COUNT} );
                    $ref->{"conserved domain hypothetical"}->{COUNT}++;
                }
                elsif( $val =~ /putative lipoprotein/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"putative lipoprotein"}->{COUNT} );
                    $ref->{"putative lipoprotein"}->{COUNT}++;
                }
                elsif( $val =~ /putative membrane/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"putative membrane"}->{COUNT} );
                    $ref->{"putative membrane"}->{COUNT}++;
                }
                elsif( $val =~ /putative/ ){
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"putative"}->{COUNT} );
                    $ref->{"putative"}->{COUNT}++;
                }
                else{
                    &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (hypothetical/putative)"}->{TYPE}->{"hypothetical"}->{COUNT} );
                    $ref->{"hypothetical"}->{COUNT}++;
                }
            }
            elsif( $tag =~ /product/i ){
                &init_ref( $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (Function assigned)"}->{COUNT} );
                $stats->{$temp[0]}->{$temp[2]}->{ATTR}->{"CDS (Function assigned)"}->{COUNT}++;
            }
        }
    }
}
close GFF3;

my ($text, $length) = ("",0);

foreach my $contig ( keys %$stats )
{
    print "\nAnnotation stats - $contig (".$ctg_len->{$contig}."bp)\n\n";
    $length = $ctg_len->{$contig};
    
    foreach my $feat ( sort keys %{$stats->{$contig}} )
    {
        ##skipping feat
        #next if $feat =~ /(polypeptide|exon|CDS)/;
        
        printf "%-40s%6s\n",
            $feat,
            $stats->{$contig}->{$feat}->{COUNT};

        $text .= "$feat: ".$stats->{$contig}->{$feat}->{COUNT}."\\n" if $feat eq "CDS";

        if( defined $stats->{$contig}->{$feat}->{ATTR} ){
           
            foreach my $tag ( sort keys %{$stats->{$contig}->{$feat}->{ATTR}} )
            {
                # skipping attributes
                next if $tag =~ /(ID|Name|description|gene_symbol_source|Parent|product)/i;
        
                printf "%-40s%6s\n", "  - $tag", $stats->{$contig}->{$feat}->{ATTR}->{$tag}->{COUNT};
               
                $text .= "  - $tag: ".$stats->{$contig}->{$feat}->{ATTR}->{$tag}->{COUNT}."\\n" if $feat eq "CDS";

                foreach my $type ( sort keys %{$stats->{$contig}->{$feat}->{ATTR}->{$tag}->{TYPE}} ){
                    my $value = $stats->{$contig}->{$feat}->{ATTR}->{$tag}->{TYPE}->{$type}->{COUNT};
                    $value = 0 unless defined $value;

                    printf "%-40s%6s\n", "    * $type", $value;
                    
                    $text .= "    * $type: $value\\n" if $feat eq "CDS";
                }
            }
        }
    }
    
    #plot
    print STDERR "\nGenerate plots in R...";
    $ctg_num = 1 unless $ignore_contig;
    Rplot($opt{input}, $opt{title}, $text, $opt{prefix}, $opt{outfmt}, $ctg_num, $length);
    print STDERR "Done.\n";
}

sub Rplot {
    my ($input, $title, $text, $prefix, $outfmt, $ctg_num, $length) = @_;

    my $cwd = `pwd`;
    chomp($cwd);

    my $outcmd = "pdf(file='${prefix}_plots.pdf',width=10,height=8);";
    $outcmd = "png('${prefix}_plots.png',height=800, width=1000, type='cairo');" if $outfmt =~ /png/i;

    my $Rscript = "
setwd('$cwd');
$outcmd
#page setup
par(mfrow=c(1,2), oma=c(2,1.5,2,0))
layout(matrix(c(1,2),nrow=1), widths=c(1,3))

#read gff3
gffRead <- function(gffFile, nrows = -1) {
  raw = read.table(gffFile, sep='\t', fill=TRUE, as.is=TRUE, quote='',
  header=FALSE, comment.char='#', nrows = nrows,
  colClasses=c('character', 'character', 'character', 'integer',
               'integer', 'character', 'character', 'character', 'character'))
  colnames(raw) = c('seqname', 'source', 'feature', 'start', 'end',
                    'score', 'strand', 'frame', 'attributes')
  gff = raw[!is.na(raw\$start),]
       
  #add length column
  gff['length']<-NA
  gff\$length<-gff\$end-gff\$start+1
  return(gff)
}

gff<-gffRead('$input')

barc = 'grey'
#feature plot
table.feat<-table(gff\$feature)
barplot(table.feat,ylab=expression('Count  log'[10]), ylim=c(1,10000), col=barc, log='y', las=3, main='Feature count')
#protein size distribution
dist<-hist(gff\$length/3,
         xlab='Protein size (aa)',
         ylab='Count',
         main='Distribution of protein size',
         breaks=seq(0, max(gff\$length)/3+100, 100), 
         xlim=c(0,max(gff\$length)/3+400),
         probability=FALSE, col=barc, border='black', cex=0.8 )
#detail info
text(max(gff\$length)/3/2,max(dist\$counts),adj=c(0,1),cex=0.8,'$text')
title('ANNOTATION STATS - $title ($ctg_num contigs, ${length}bp)', outer=TRUE)

tmp<-dev.off();

quit();
";

    open (Rscript, ">/tmp/Rscript$$") or die "Can't write to /tmp/Rscript$$: $!\n";
    print Rscript $Rscript;
    close Rscript;

    if (system ("R --vanilla --slave --silent < /tmp/Rscript$$ 2>/dev/null")) {warn "$!\n"};
    unlink "/tmp/Rscript$$";
    unlink "Rplots.pdf" if ( -e "Rplots.pdf");
}

sub init_ref {
    my ($ref) = @_;
    $ref = 0 unless defined $ref;
}

sub usage
{
    print  <<USAGE;
Usage:
    perl $0 --input <GFF3> --title <STRING> --prefix <STRING> --outfmt <PDF|PNG>

    --input|-i <GFF3> : input annotation 
    --title|-t <STRING> : organism name to show in R plot
    --prefix|-p <STRING> : output prefix
    --outfmt|-o <PDF|PNG>  : output format (default is pdf)

USAGE
exit;
}
