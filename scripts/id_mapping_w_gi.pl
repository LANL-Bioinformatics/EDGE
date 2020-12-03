#! /usr/bin/env perl

my $coverage_table=$ARGV[0];
my $sort_method=$ARGV[1];
my $EDGE_HOME=$ENV{'EDGE_HOME'};
my $id_mapping_file="$EDGE_HOME/database/bwa_index/id_mapping.txt";
if (!@ARGV) {print "perl $0 coverageTable sort_method
    sort_method: fold std cov reads bases(default:cov)\n"; exit;}

if (!$sort_method){$sort_method="cov";}

my (%mapping,$mapping);
$mapping= &load_mapping_file($id_mapping_file,\%mapping);
%mapping=%{$mapping};
#$mapping=&load_mapping_file($id_mapping_file2,\%mapping);
#%mapping=%{$mapping};

sub load_mapping_file 
{
   my ($file,$mapping)= @_;
   my %hash=%{$mapping};
   open (IN,$file) or die "$!";

   while(<IN>)
   {
     chomp;
     my ($id,@name)=split /\s+/,$_;
     my $name = join (" ",@name);
     $hash{"$id"}=$name;
   }
   close IN;
   return \%hash;
}
#ID	Length	GC%	Avg_fold	Fold_std	Base_Coverage%	Mapped_reads
open (IN, $coverage_table) or die "$!";
my $head=<IN>;
my $head = "Accession\tLength\tGC%\tAvg_fold\tFold_std\tBase_Coverage%\tMapped_reads\tLinear_length\tID\n";
my %coverage;
while (<IN>)
{
   chomp;
   my @array = split /\t/,$_;
   $coverage{$array[0]}{len}=$array[1];
   $coverage{$array[0]}{cg}=$array[2];
   $coverage{$array[0]}{fold}=$array[3];
   $coverage{$array[0]}{std}=$array[4];
   $coverage{$array[0]}{cov}=$array[5];
   $coverage{$array[0]}{reads}=$array[6];
   $coverage{$array[0]}{bases}=int($array[5]*$array[1]/100);
}
close IN;

print $head;
foreach my $id (sort {$coverage{$b}{$sort_method}<=>$coverage{$a}{$sort_method}} keys %coverage)
{
     next if ($coverage{$id}{len}<1000);
     next if ($coverage{$id}{reads}==0);
     
     #my ($gi) = $id =~ /gi\|(\d+)\|/;
     my $acc = getAccFromSeqID($id);
     print $acc,"\t";
     #print $id,"\t";
     print $coverage{$id}{len},"\t";
     print $coverage{$id}{cg},"\t";
     print $coverage{$id}{fold},"\t";
     print $coverage{$id}{std},"\t";
     print $coverage{$id}{cov},"\t";
     print $coverage{$id}{reads},"\t";
     print $coverage{$id}{bases},"\t";
     print $mapping{"$id"},"\n";
}

sub getAccFromSeqID
{
	my ($seqID) = @_;
	
	$seqID =~ /^>?(\S+)/;
	
	my $acc = $1;
	
	if ( $acc =~ /\|/ )
	{
		$acc = (split /\|/, $acc)[3];
	}
	
	if ( $acc !~ /^\d+$/ && $acc !~ /^[A-Z]+_?[A-Z]*\d+(\.\d+)?$/ )
	{
		$invalidAccs{$acc} = 1;
		return undef;
	}
	
	return $acc;
}
