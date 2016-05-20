#!/usr/bin/perl -w

use strict;
use FileHandle;
use diagnostics;

my $dir= shift @ARGV;
my $output=shift @ARGV;

opendir (DIR, $dir);

my $file;
my %sequences;
my ($header,@seq);
foreach my $files(sort{$a cmp $b} readdir(DIR)){
   next if ($files=~ /^..?$/);
   if ($files=~ /\.cdn$/){
      print "$files\n";
      $file= $dir.'/'.$files;
   }
   if (!-z $file){read_file($file);}
}
open (OUT, ">$output");
for my $keys (sort keys %sequences) {   
   my $alignment= join '', @{$sequences{$keys}};
   print OUT ">$keys\n$alignment\n";
}


sub read_file
{
my $alignment= shift;
my $name;
 
my $fh= FileHandle->new($alignment)|| die "$!";
if ($fh->open("< $alignment")){
   $/= ">";
   while (<$fh>){
      $_=~ s/\>//g;
      unless ($_){next};
      ($header,@seq) = split /\n/, $_;
      my $line= join "",@seq;
      push @{$sequences{$header}},$line;
   }
   $/="\n";
   $fh->close;
}
}
