#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($MEGAN_CSV, $OUTDIR, $PREFIX, $CWD);

GetOptions(
	'cwd|c=s',    \$CWD,
	'outdir|o=s', \$OUTDIR,
	'prefix|p=s', \$PREFIX,
	'megan|m=s',  \$MEGAN_CSV
);

my $x = "
set dir=$CWD;
import csv=summary separator=tab file=$MEGAN_CSV;
show window=mainviewer;
set nodedrawer=piechart;
nodelabels assigned=true;
update;
exportimage file=$OUTDIR/$PREFIX.megan.pdf format=pdf visibleOnly=false replace=true title=$PREFIX;
save file='$OUTDIR/$PREFIX.megan' summary=true;
quit;
";

print STDERR "$x\n";

$x =~ s/\n//mg;

if( -s $MEGAN_CSV ){
	my $cmd = "MEGAN +g -x \"$x\"";
	system($cmd);
}
else{
	print STDERR "\n\n";
	print STDERR "********************************************************************\n";
	print STDERR "** WARRNING: FILE IS EMPTY!!\n";
	print STDERR "** FILE: $MEGAN_CSV \n";
	print STDERR "********************************************************************\n\n";
}
