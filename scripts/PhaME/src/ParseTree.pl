#!/usr/bin/perl -w 

use strict;
use File::Basename; 

my $previous="";
my $count=1;
my $file=shift || die "Tree file (newick format) to parse";
my $size=0;
my $rawtree;
my $node=0;
my %child;
my $offset=0;

my ($name,$path,$suffix)=fileparse("$file",qr/\.[^.]*/);
open(IN, $file) || die "Can't open $file";
while (<IN>){
   chomp; 
   $rawtree=$_;
   while (length $_>0){
      my $char=chop($_);
      if($char=~/:/){$node++;}
   }
}
close IN;

$size=length $node;
#$node--;

FromStringWithDistance($rawtree);

$rawtree=~ s/:\d.\d+,/,/g;
$rawtree=~ s/:\d.\d+\)/)/g;
$rawtree=~ s/:0.0;/;/g;
$rawtree=~ s/:\d.\d+E-\d,/,/g;
$rawtree=~ s/:\d.\d+E-\d\)/)/g;

my $branch=$path.$name.'_BranchNumber'.$suffix;
open (OUT,">$branch");

$offset=length $rawtree;
my $comma=rindex($rawtree,",",$offset);
my $paren=rindex($rawtree,")",$offset);
my $number=1;
while ($comma!= -1 || $paren!=-1){
   if ($paren>$comma){
      substr($rawtree,$paren,0)=":$number";
      $offset=$paren-1;
   }
   else{
      substr($rawtree,$comma,0)=":$number";
      $offset=$comma-1;
   }
   $number++;
   $comma=rindex($rawtree,",",$offset);
   $paren=rindex($rawtree,")",$offset);
}
print OUT $rawtree;
close OUT;

sub FromStringWithDistance
{
my $in=shift;
my $distance;
my $rparen=rindex($in,")");
my $rcolon=rindex($in,":");
my $length=length $rawtree;
$length=$length-$rcolon;
if ($rcolon>1){
   my $outfile=$path.$name.'_';
   $outfile.=sprintf "%0${size}d", $count;
   $outfile.=$suffix;
   open (OUT,">$outfile");

   my $string=substr($rawtree,0,$rcolon); # adjust the string to exclude the distance
   $previous=substr($rawtree,$rcolon,$length);

   if ($previous!~/^:0\.0;$/){
      print OUT "$string#1$previous\n";
      close OUT;
      $count++;
   }
   FromStringWithDistance($string);
}
}


