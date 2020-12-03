#! /usr/bin/env perl

if (!@ARGV) {print "perl $0 <fasta>\n";exit;}
my %chr;
my $len;
my %all_chr;
while(<>)
{
   chomp;
   if (/>(\S+)/)
   {
	   $id=$1;
   }
   else
   {
      @letter = split //,$_;
      foreach my $nu (@letter)
      {
           $chr{$id}{$nu}++;
	   $all_chr{$nu}++;
      }
   }
}

&print_by_id(\%chr);
if(scalar(keys %chr)>1){
  &print_all(\%all_chr);
}

sub print_by_id{
  my $hash = shift;
  my %chr = %{$hash};
  foreach my $id(sort keys %chr){
    print "##$id\n";
    my $total_GC=0;
    my $total_len_no_gap=0;
    my $len=0;
    foreach my $nu ("A","T","C","G","N","a","t","c","g","n"){
        if ($chr{$id}{$nu}){
          my $count = $chr{$id}{$nu};
          print $nu, "\t", $count,"\n";
          $total_GC += $count if ($nu eq "C" or $nu eq "G" or $nu eq "c" or $nu eq "g" );
          $total_len_no_gap += $count;
	  $len += $count;
          delete $chr{$id}{$nu};
        }
    }
    my %chr2 = %{$chr{$id}}; 
    foreach my $nu (sort {$chr2{$b} <=> $chr2{$a} } keys %chr2)
    {
        print $nu, "\t", $chr2{$nu},"\n";
	$len += $chr2{$nu};
        $total_len_no_gap += $chr2{$nu} if ($nu ne "-");
    
    }
    print "GC_content\t". $total_GC/$len ."\n";
    print "Total length\t$len\n";
    print "Total no gap bases\t$total_len_no_gap\n\n";
  }
}

sub print_all{
  my $hash = shift;
  my %chr = %{$hash};
  my $total_len_no_gap=0;
  my $total_GC=0;
  my $len=0;
  print "##ALL\n";
  foreach my $nu ("A","T","C","G","N","a","t","c","g","n"){
    if ($chr{$nu}){
	    my $count = $chr{$nu};
	    print $nu, "\t", $count,"\n";
	    $total_GC += $count if ($nu eq "C" or $nu eq "G" or $nu eq "c" or $nu eq "g" );
	    $total_len_no_gap += $count;
	    $len += $count;
	    delete $chr{$nu};
    }
  }
  foreach my $nu (sort {$chr{$b} <=> $chr{$a} } keys %chr){
	  print $nu, "\t", $chr{$nu},"\n";
	  $len += $chr{$nu};
	  $total_len_no_gap += $chr{$nu} if ($nu ne "-");
  }
  print "GC_content\t". $total_GC/$len ."\n";
  print "Total length\t$len\n";
  print "Total no gap bases\t$total_len_no_gap\n";
}
