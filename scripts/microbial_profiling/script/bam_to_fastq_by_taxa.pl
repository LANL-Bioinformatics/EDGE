#! /usr/bin/perl
# require samtools 
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin qw($RealBin);
use lib $RealBin;
use gi2lineage;


my $workingDir = Cwd::getcwd();
my $in_offset=33;
my $out_offset=33;
my $opt_mapped;
my $opt_unmapped;
my $opt_passfilter;
my $prefix = "Reads";
my $mapped_ref_id=""; # this is the id in SAM 'RNAME' field (3rd column);
my $sam_format;
my $rank = "genus";
my $mapped_ref_name=""; # this is the name in a taxonomy rank.
my $preload;
my $single_end;
my $zip;
my $fastq_file;
my %reads_id;

GetOptions(  
            'in_offset=i'  =>  \$in_offset,
            'out_offset=i' =>  \$out_offset,
			'sam|S'        =>  \$sam_format,
 #           'type=s'       =>  \$type, 
            'mapped'       =>  \$opt_mapped,
            'id=s'         =>  \$mapped_ref_id,
            'rank=s'       =>  \$rank,
            'name=s'       =>  \$mapped_ref_name,
			'fastq=s'      =>  \$fastq_file,
            'preload'      =>  \$preload,
            'unmapped'     =>  \$opt_unmapped,
            'pf'           =>  \$opt_passfilter,
            'se'           =>  \$single_end,
            'zip'          =>  \$zip,
            'prefix=s'     =>  \$prefix,
            'help|?'       =>  sub{Usage()}
);
my $time = time;
if (!@ARGV){&Usage}
if ($in_offset and !$out_offset) { die "Please provide output offset for converstion.\n";} 
if (!$in_offset and $out_offset) { die "Please provide iutput offset for converstion.\n";} 
if ($opt_mapped and $opt_unmapped) { die "Please provide either -mapped or -unmapped.\n";}
my $file=$ARGV[0];
my $samtools_flag="";
if ($opt_mapped || $mapped_ref_name || $mapped_ref_id) {$samtools_flag="-F 4"};
if ($opt_unmapped) {$samtools_flag="-f 4"};
open (my $pair1_fh, ">$prefix.1.fastq");
open (my $pair2_fh, ">$prefix.2.fastq");
open (my $se_fh, ">$prefix.se.fastq");
my $sortName_sam_file = "unmapped$$.sam";
$mapped_ref_id = "\"$mapped_ref_id\"" if ($mapped_ref_id);
$samtools_flag .= " -S" if ($sam_format);

$preload = 'preload' if ($preload);
if ($mapped_ref_name){
  &print_timeInterval($time, "Load Tanonomy info ... ");
  loadTaxonomy($preload)
};

if ( ! $single_end){
  &print_timeInterval($time, "Sorting the SAM file ... ");
  system("samtools view $samtools_flag $file $mapped_ref_id | sort -T $workingDir -k 1,1 > $sortName_sam_file");
  &print_timeInterval($time, "Parsing the SAM file ...");
  open (IN, $sortName_sam_file) or die $!;
}else{ # skip sorting to speed up process 
  &print_timeInterval($time, "Parsing the SAM file ...");
  open (IN, "samtools view $samtools_flag $file $mapped_ref_id|") or die $!;
}

my $p1_count=0;
my $p2_count=0;
my $se_count=0;
while(<IN>)
{
  chomp;
  next if /^\@SQ/;
  my @array=split /\t/,$_;
  next if ($opt_passfilter and ($array[1] & 512));

  if ($mapped_ref_name){
    my $acc = getAccFromSeqID($array[2]);
    my $acc_to_name = acc2rank($acc,lc($rank));
    next if ( lc($acc_to_name) ne lc($mapped_ref_name));
  }
  if ($fastq_file){
    $reads_id{$array[0]}=1;
    ( my $tmp  = $array[0]) =~ s/_\d+$//;
    $reads_id{$tmp}=1;
    next;
  }

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

&extract_from_original_fastq($fastq_file,\%reads_id) if ($fastq_file);

close $pair1_fh;
close $pair2_fh;
close $se_fh;

print "$rank: " if ($rank);
print "$mapped_ref_name\n" if ($mapped_ref_name);
print "Paired 1: $p1_count; ";
print "Paired 2: $p2_count; ";
print "Single End: $se_count; \n";

unlink "$prefix.1.fastq" if ( -z "$prefix.1.fastq");
unlink "$prefix.2.fastq" if ( -z "$prefix.2.fastq");
unlink "$prefix.se.fastq" if ( -z "$prefix.se.fastq");
unlink  $sortName_sam_file if ( ! $single_end);
if ($zip){
  &print_timeInterval($time, "Compressed fastq files ... \n");
  my ($file_name, $file_path, $file_suffix)=fileparse("$prefix", qr/\.[^.]*/);
  chdir $file_path;
  `zip $file_name.fastq.zip $file_name*fastq`;
  `rm -f $file_name*fastq`;
}

&print_timeInterval($time, "All Done.\n");


sub extract_from_original_fastq {
	my $fastq = shift;
	my $id = shift;
	open (my $fh, $fastq) or die "Cannot read $fastq $!\n";
	while(<$fh>){
		my $s_id=$_;
		my $seq=<$fh>;
		my $q_id=<$fh>;
		my $q_seq=<$fh>;
		if($s_id=~ /^\@(\S+)/){
			if ($id->{$1}){
				print $se_fh $s_id.$seq.$q_id.$q_seq;
				$se_count++;
			}
		}
	}
	close $fh;
}

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
-sam|S          Input is sam format.
-in_offset      Offset nubmer for ASCII encoding from input
-out_offset     Offset nubmer for ASCII encoding for output 
                Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)
                Use this option to convert the encoding offset.
                default: no conversion.
-mapped         For retrieving mapped reads from alignment bam file
                -id   retrieving reads mapped to this reference ID. (e.g. contig_00001)
                  Or   retiveving reads mapped to a specified taxonomy name
                -rank  phylum, class, order, family, genus, species or strain (default: genus)
                -name  the name in a taxonomy rank (e.g. Pseudomonas)
 
-unmapped       For retrieving unmapped reads from alignment bam file
-pf             Output passFilter reads only
-se             The sam file is from single end reads, will bypass sorting to speed up.
-zip            Compress output files with tar gz.
-fastq          Original mapping fastq file. Will use this to extract reads for single end reads only
-pefix          Output file prefix (default: Reads)
-help           Show this usage

END
#-type           pe or se.  pe=paired end, se=single end. The pe type will output reads' name with /1 or /2.

exit;
}

sub print_timeInterval{
    my $now = shift;
    my $msg = shift;
    $now = time - $now;
    my $string=sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
    print STDERR "[$string]  $msg\n";
}
