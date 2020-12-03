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
'AAR'=>'K', 'AAS'=>'X', 'AAY'=>'N', 'ACN'=>'T', 'ACR'=>'T', 'ACY'=>'T', 'AKC'=>'X', 'AKT'=>'X', 'AMT'=>'X', 'ARA'=>'X', 'ARC'=>'X', 'ART'=>'X', 'AST'=>'X', 'ATR'=>'X', 'ATY'=>'I', 'AYT'=>'X', 'CAK'=>'X', 'CAN'=>'X', 'CAR'=>'Q', 'CCK'=>'P', 'CCN'=>'P', 'CGM'=>'R', 'CGN'=>'R', 'CGK'=>'R', 'CGY'=>'R', 'CKM'=>'X', 'CMA'=>'X', 'CMC'=>'X', 'CTY'=>'L', 'CVC'=>'X', 'CWG'=>'X', 'CYG'=>'X', 'CYT'=>'X', 'DAG'=>'X', 'GAK'=>'X', 'GAN'=>'X', 'GAR'=>'E', 'GAS'=>'X', 'GAY'=>'D', 'GCN'=>'A', 'GCR'=>'A', 'GCS'=>'A', 'GCY'=>'A', 'GGK'=>'G', 'GGY'=>'G', 'GKA'=>'X', 'GMT'=>'X', 'GNT'=>'X', 'GRC'=>'X', 'GSG'=>'X', 'GST'=>'X', 'GTS'=>'V', 'GTW'=>'V', 'GTY'=>'V', 'GYG'=>'X', 'GYT'=>'X', 'KAA'=>'X', 'KAG'=>'X', 'KCC'=>'X', 'KGC'=>'X', 'MTC'=>'X', 'MTG'=>'X', 'MTT'=>'X', 'NAC'=>'X', 'NTG'=>'X', 'RAA'=>'X', 'RAC'=>'X', 'RAG'=>'x', 'RAT'=>'X', 'RAW'=>'X', 'RCG'=>'X', 'RGC'=>'X', 'RGT'=>'X', 'RTA'=>'X', 'RTC'=>'X', 'RTG'=>'X', 'RTT'=>'X', 'SAG'=>'X', 'SCC'=>'X', 'SGT'=>'X', 'SGS'=>'X', 'SGY'=>'X', 'STG'=>'x', 'TCM'=>'S', 'TCR'=>'S', 'TGM'=>'X', 'TGY'=>'C', 'TST'=>'X', 'TTN'=>'X', 'TTS'=>'X', 'TTW'=>'X', 'TTY'=>'X', 'TWC'=>'X', 'TWG'=>'X', 'WAA'=>'X', 'WTT'=>'X', 'YGT'=>'X', 'YTC'=>'X', 'YTG'=>'X', 'YTY'=>'X',
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
