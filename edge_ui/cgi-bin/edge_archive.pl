#!/usr/bin/env perl
use File::Basename;
use Getopt::Long;
use Term::ANSIColor;
use POSIX qw(strftime);
use Cwd;
use FindBin qw($RealBin);
use lib "$RealBin/lib";

my ($proj, $adir) = @ARGV;
my $time = time;
exit if !@ARGV;
&stringSanitization($proj, $adir);
chdir($proj);
my $command = "tar cf - ./ | ( cd $adir; tar xfp -) >> $proj/error.log 2>&1";

if( -d $adir ){
	`echo "Archiving command failed. The directory is already exist." |tee -a $proj/process.log >> $proj/process_current.log`;
	exit;
}
else{
	`mkdir -p $adir`;
}

if ( system($command) != 0 ){
	`echo "Archiving command failed." |tee -a $proj/process.log >> $proj/process_current.log`;
}
else{
	`rm -rf $proj`;
	`ln -s $adir $proj`;
	my $runningtime = &runTime($time);
	`echo "$runningtime" |tee -a $proj/process.log >> $proj/process_current.log`;
	`echo "All Done." |tee -a $proj/process.log >> $proj/process_current.log`;
	`echo "\n*** EDGE_UI: project archived ***" |tee -a $proj/process.log >> $proj/process_current.log`;
}

sub runTime {
  my $time=shift;
  my $runTime = time() - $time;
  my $time_string = sprintf(" Running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime % 3600) / 60), 
  int($runTime % 60));
  return $time_string;
}


sub stringSanitization{
	my @values=@_;
	my $dirtybit=0;
	foreach my $str (@values){
		$dirtybit=1 if ($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/);
		if ($dirtybit){
			print "Content-Type: text/html\n\n", "Invalid characters detected \'$str\'.\n\n";
			exit;
		}
	}
}
