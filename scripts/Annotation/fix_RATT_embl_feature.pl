#!/usr/bin/env perl
#
#
open (IN, $ARGV[0]);
while(<IN>){
	if (/order|complement|join\(/){ 
		$line = $_;
		$line =~ s/join\(order/order/;
		$line =~ s/order\(join/join/;
		#my $extra;
		#INNER: while($nline=<IN>){
	#		if ($nline =~ /\//){
	#			$extra = $nline;
	#			last;
	#		}
	#		$line .= $nline;
	#	}
		($op)=$line=~ tr/\(/\(/;
		($cp)=$line=~ tr/\)/\)/;  
		if ($op > $cp){ 
			chomp $line;
			print $line. ")" x ($op - $cp) . "\n";
	#		print $extra if ($extra);
		}else{
			print $line;
		}
	}else{
		print $_;
	}
}
close IN;
