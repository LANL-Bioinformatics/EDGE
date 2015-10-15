#!/usr/bin/perl -w

#A program to translate a nucleotide sequence to an amino acid sequence

use strict;

my ($id,$values,$seq,$length,$dna,$start,$index,$codon,$codlen)=0;

my $file=shift @ARGV;

open(FILE, "$file")||die "$!";

my $protein;

while(<FILE>){
   chomp;
   if (/^(>\S+)/){
      $id=$1;
      print "$id\n";
   }
   if (! /^>/){
      $protein='';
      $seq=$_;
      $length=length $seq; 	               #length of sequence
      foreach ($index=0; $index <= $length; $index+=3){
         $codon=substr($seq,$index,3);
         $codlen=length $codon;           #avoiding errors when two bases left for last codon
         $protein .=&codon2aa($codon), unless $codlen <=2;
#         print $dict{$codon}, unless $codlen <=2;
      }
      print "$protein\n";
   }
}

#Codon Table
sub codon2aa
{
my($codon)=@_;
$codon=uc $codon;
my %dict=(
   'AAA'=>'K', 'AAC'=>'N', 'AAG'=>'K', 'AAT'=>'N', 'ACA'=>'T', 'ACC'=>'T', 'ACG'=>'T', 'ACT'=>'T', 'AGA'=>'R', 'AGC'=>'S', 'AGG'=>'R', 'AGT'=>'S', 'ATA'=>'I', 'ATC'=>'I', 'ATG'=>'M', 'ATT'=>'I', 'CAA'=>'Q', 'CAC'=>'H', 'CAG'=>'Q', 'CAT'=>'H', 'CCA'=>'P', 'CCC'=>'P', 'CCG'=>'P', 'CCT'=>'P', 'CGA'=>'R', 'CGC'=>'R', 'CGG'=>'R', 'CGT'=>'R', 'CTA'=>'L', 'CTC'=>'L', 'CTG'=>'L', 'CTT'=>'L', 'GAA'=>'E', 'GAC'=>'D', 'GAG'=>'E', 'GAT'=>'D', 'GCA'=>'A', 'GCC'=>'A', 'GCG'=>'A', 'GCT'=>'A', 'GAC'=>'D', 'GAG'=>'E', 'GAT'=>'D', 'GCA'=>'A', 'GCC'=>'A', 'GCG'=>'A', 'GCT'=>'A', 'GGA'=>'G', 'GGC'=>'G', 'GGG'=>'G', 'GGT'=>'G', 'GTA'=>'V', 'GTC'=>'V', 'GGC'=>'G', 'GGG'=>'G', 'GGT'=>'G', 'GTA'=>'V', 'GTC'=>'V', 'GTG'=>'V', 'GTT'=>'V', 'TAC'=>'Y', 'TAT'=>'Y', 'TCA'=>'S', 'TCC'=>'S', 'TCG'=>'S', 'TCT'=>'S', 'TGC'=>'C', 'TGG'=>'W', 'TGT'=>'C', 'TTA'=>'L', 'TTC'=>'F', 'TTG'=>'L', 'TTT'=>'F', 'TAA'=>'*', 'TAG'=>'*', 'TGA'=>'*',
'AAR'=>'K', 'AAS'=>'', 'AAY'=>'N', 'ACN'=>'T', 'ACR'=>'T', 'ACY'=>'T', 'AKC'=>'', 'AKT'=>'', 'AMT'=>'', 'ARA'=>'', 'ARC'=>'', 'ART'=>'', 'AST'=>'', 'ATR'=>'', 'ATY'=>'I', 'AYT'=>'', 'CAK'=>'', 'CAN'=>'', 'CAR'=>'Q', 'CCK'=>'P', 'CCN'=>'P', 'CGM'=>'R', 'CGN'=>'R', 'CGK'=>'R', 'CGY'=>'R', 'CKM'=>'', 'CMA'=>'', 'CMC'=>'', 'CTY'=>'L', 'CVC'=>'', 'CWG'=>'', 'CYG'=>'', 'CYT'=>'', 'DAG'=>'', 'GAK'=>'', 'GAN'=>'', 'GAR'=>'E', 'GAS'=>'', 'GAY'=>'D', 'GCN'=>'A', 'GCR'=>'A', 'GCS'=>'A', 'GCY'=>'A', 'GGK'=>'G', 'GGY'=>'G', 'GKA'=>'', 'GMT'=>'', 'GNT'=>'', 'GRC'=>'', 'GSG'=>'', 'GST'=>'', 'GTS'=>'V', 'GTW'=>'V', 'GTY'=>'V', 'GYG'=>'', 'GYT'=>'', 'KAA'=>'', 'KAG'=>'', 'KCC'=>'', 'KGC'=>'', 'MTC'=>'', 'MTG'=>'', 'MTT'=>'', 'NAC'=>'', 'NTG'=>'', 'RAA'=>'', 'RAC'=>'', 'RAG'=>'', 'RAT'=>'', 'RAW'=>'', 'RCG'=>'', 'RGC'=>'', 'RGT'=>'', 'RTA'=>'', 'RTC'=>'', 'RTG'=>'', 'RTT'=>'', 'SAG'=>'', 'SCC'=>'', 'SGT'=>'', 'SGS'=>'', 'SGY'=>'', 'STG'=>'', 'TCM'=>'S', 'TCR'=>'S', 'TGM'=>'', 'TGY'=>'C', 'TST'=>'', 'TTN'=>'', 'TTS'=>'', 'TTW'=>'', 'TTY'=>'', 'TWC'=>'', 'TWG'=>'', 'WAA'=>'', 'WTT'=>'', 'YGT'=>'', 'YTC'=>'', 'YTG'=>'', 'YTY'=>'',
);
if(exists $dict{$codon}){return $dict{$codon};}
else{
   print STDERR "Bad codon: $codon\n";
   return '';
}
}

=head


&translate($seq, 1);
}

exit;

sub translate {
 $dna = shift @_;
 $start = shift @_;

# print $start,"\t\t";
  foreach ($index=abs ($start)-1; $index <= $length; $index+=3){
    $codon = substr ($dna, $index, 3);
    $codlen = length $codon; #avoiding errors when two bases left for last codon
    print $dict{$codon}, unless $codlen <=2;
  }

  print "\n";
}
=cut

close FILE;

exit;
