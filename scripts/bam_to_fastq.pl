#! /usr/bin/perl
# require samtools 
use strict;
use warnings;
use Getopt::Long;
use Cwd;

my $workingDir = Cwd::getcwd();
my $in_offset=33;
my $out_offset=33;
my $opt_mapped;
my $opt_unmapped;
my $opt_passfilter;
my $prefix = "Reads";
my $mapped_ref_id="";
GetOptions(  
            'in_offset=i'  =>  \$in_offset,
            'out_offset=i' =>  \$out_offset,
 #           'type=s'       =>  \$type, 
            'mapped'       =>  \$opt_mapped,
	    'id=s'	   =>  \$mapped_ref_id, 
            'unmapped'     =>  \$opt_unmapped,
            'pf'           =>  \$opt_passfilter,
            'prefix=s'     =>  \$prefix,
            'help|?'       =>  sub{Usage()}
);

if (!@ARGV){&Usage}
if ($in_offset and !$out_offset) { die "Please provide output offset for converstion.\n";} 
if (!$in_offset and $out_offset) { die "Please provide iutput offset for converstion.\n";} 
if ($opt_mapped and $opt_unmapped) { die "Please provide either -mapped or -unmapped.\n";}
my $file=$ARGV[0];
my $samtools_flag="";
if ($opt_mapped) {$samtools_flag="-F 4"};
if ($opt_unmapped) {$samtools_flag="-f 4"};
open (my $pair1_fh, ">$prefix.1.fastq");
open (my $pair2_fh, ">$prefix.2.fastq");
open (my $se_fh, ">$prefix.se.fastq");
my $sortName_sam_file = "unmapped$$.sam";
$mapped_ref_id = "\"$mapped_ref_id\"" if ($mapped_ref_id);
system("samtools view $samtools_flag $file $mapped_ref_id | sort -T $workingDir -k 1,1 > $sortName_sam_file");
open (IN, $sortName_sam_file) or die $!;
my $p1_count=0;
my $p2_count=0;
my $se_count=0;
while(<IN>)
{
  chomp;
  my @array=split /\t/,$_;
  next if ($opt_passfilter and ($array[1] & 512));
  if ($in_offset != $out_offset)
  {
     $array[10]=&quality_conversion($array[10],$in_offset,$out_offset);
  }
  if($array[1] & 1)    
  {
      if ( ($opt_mapped and ($array[1] & 8))  or # the other mate unmapped reads
           ($opt_unmapped and !(($array[1] & 4) && ($array[1] & 8)))
      )  
      {
          $se_count++;
          $array[9]=ReverseComplement($array[9]) if ($array[1] & 16);
          $array[10]=reverse($array[10])if ($array[1] & 16);
          print $se_fh "@",$array[0],"\n";
          print $se_fh $array[9],"\n";
          print $se_fh "+\n";
          print $se_fh $array[10],"\n";
      }else{
          $array[9]=ReverseComplement($array[9]) if ($array[1] & 16);
          $array[10]=reverse($array[10])if ($array[1] & 16);
          if ($mapped_ref_id &&  $array[6] ne "="){
              $se_count++;
              print $se_fh "@",$array[0],"\n";
              print $se_fh $array[9],"\n";
              print $se_fh "+\n";
              print $se_fh $array[10],"\n";
              next;	
          }
          if ($array[1] & 64)
          {
              $p1_count++;
              print $pair1_fh "@",$array[0],"/1\n";
              print $pair1_fh $array[9],"\n";
              print $pair1_fh "+\n";
              print $pair1_fh $array[10],"\n";
          }
          elsif($array[1] & 128)
          {
              $p2_count++;
              print $pair2_fh "@",$array[0],"/2\n";
              print $pair2_fh $array[9],"\n";
              print $pair2_fh  "+\n";
              print $pair2_fh $array[10],"\n";
          }
      }
  }
  else
  { 
     $se_count++;
     $array[9]=ReverseComplement($array[9]) if ($array[1] & 16);
     $array[10]=reverse($array[10])if ($array[1] & 16);
     print $se_fh "@",$array[0],"\n";
     print $se_fh $array[9],"\n";
     print $se_fh "+\n";
     print $se_fh $array[10],"\n";
  }
}
close IN;
close $pair1_fh;
close $pair2_fh;
close $se_fh;


print "Paired 1: $p1_count\n";
print "Paired 2: $p2_count\n";
print "Single End: $se_count\n";

unlink "$prefix.1.fastq" if ( -z "$prefix.1.fastq");
unlink "$prefix.2.fastq" if ( -z "$prefix.2.fastq");
unlink "$prefix.se.fastq" if ( -z "$prefix.se.fastq");
unlink $sortName_sam_file;


sub quality_conversion
{
   my ($seq,$in_offset,$out_offset) = @_;
   my $seq_len=length $seq;
   my ($q, $cov_q);
   for (0..($seq_len-1)){
     $q = substr($seq,$_,1);
     $cov_q = chr (ord($q) - $in_offset + $out_offset);
     substr($seq,$_,1,$cov_q);
   }
   return ($seq);
}

sub ReverseComplement{
        my $dna = $_[0];
        my $ReverseCompSeq = reverse ($dna);
        $ReverseCompSeq =~ tr/atgcrywsmkATGCRYWSMK/tacgyrswkmTACGYRSWKM/;
        return($ReverseCompSeq);
}

sub Usage {
print <<END;
Usage: $0  <bam_file>
# require samtools in PATH.

Options:
-in_offset      Offset nubmer for ASCII encoding from input
-out_offset     Offset nubmer for ASCII encoding for output 
                Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)
                Use this option to convert the encoding offset.
                default: no conversion.
-mapped         For retrieving mapped reads from alignment bam file
		-id   retrieving reads mapped to this reference ID. (e.g. contig_00001)
 
-unmapped       For retrieving unmapped reads from alignment bam file
-pf             Output passFilter reads only
-prefix         Output file prefix (default: Reads)
-help           Show this usage

END
#-type           pe or se.  pe=paired end, se=single end. The pe type will output reads' name with /1 or /2.

exit;
}
