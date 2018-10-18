#!/usr/bin/env perl
# Los Alamos National Lab.
# 2016/10/25
#

use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use POSIX qw(strftime);
use Getopt::Long;
######################################################################################
# DATA STRUCTURE:
#
#     $list->{1}->{NAME}    // project name
#               ->{TIME}    // submission time
#               ->{STATUS}  // status [finished|running|unstarted]
#          ->{2}...
#
######################################################################################

exit if ( $ENV{"REQUEST_METHOD"} );

my $debug=0;
GetOptions("debug" => \$debug);

# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $out_dir     = $ARGV[0] || $sys->{edgeui_output};
my $edge_total_cpu = $sys->{"edgeui_tol_cpu"};
my $max_num_jobs = $sys->{"max_num_jobs"};

#cluster
my $cluster 	= $sys->{cluster};

my $list; # ref for project list

my ($memUsage, $cpuUsage, $diskUsage) = &getSystemUsage();
print join("\n",$memUsage, $cpuUsage, $diskUsage),"\n" if ($debug);
if ( ($memUsage > 99 or $cpuUsage > 99 or $diskUsage > 98) and !$cluster){
	print STDERR "No enough computation resource to run CPU/MEM/DISK\n";
	exit;
}


#check projects vital
my ($vital, $name2pid, $error);
($vital, $name2pid) = &checkProjVital();

&scanNewProjToList();

#autorun
if( scalar keys %$list && $sys->{edgeui_auto_run} ){
	my ( $progs, $proj, $projCode, $p_status, $proj_start, $proj_dir, $log, $config , $domain);
	my $num_cpu_used = 0;
	foreach my $i ( sort {$list->{$a}->{TIME} cmp $list->{$b}->{TIME}} keys %$list ) {
		$proj     = $list->{$i}->{NAME};
		$projCode = $list->{$i}->{PROJCODE};
		$proj_dir = $list->{$i}->{PROJDIR};
		$domain   = $list->{$i}->{PROJDOMAIN}; 
		my $run=0;
		print $list->{$i}->{PROJNAME},"\t$list->{$i}->{TIME}\t$list->{$i}->{STATUS}\n" if  $list->{$i}->{STATUS} ne "finished" && $debug;
		if ($list->{$i}->{STATUS} eq "unstarted"){
			$run = &availableToRun($list->{$i}->{CPU}, $num_cpu_used ) if $list->{$i}->{STATUS} eq "unstarted";
			if($run or $cluster){
				#my $cmd="$RealBin/edge_action.cgi $proj rerun '' '' '' '' $domain 0 2>> $proj_dir/error.log";
				chdir $RealBin;
				my $json = `$RealBin/edge_action.cgi $proj rerun "" "" "" "" $domain 0 2>> $proj_dir/error.log`;
				print STDERR "$json\n" if ($debug);
				#print STDERR "$cmd\n" if ($debug);
				$num_cpu_used += $list->{$i}->{CPU};
			}else{
				
			}
		}
	}
}

######################################################


sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			$flag=1;
			while(<CONF>){
				chomp;
				last if /^\[/;
				if ( /^([^=]+)=([^=]+)/ ){
					$sys->{$1}=$2;
				}
			}
		}
		last;
	}
	close CONF;
	die "Incorrect system file\n" if (!$flag);
	return $sys;
}

sub scanNewProjToList {
	my $cnt = 1;
	
	opendir(BIN, $out_dir) or die "Can't open $out_dir: $!";
	my @dirfiles = readdir(BIN);
	foreach my $file (@dirfiles)  {
		next if ($file eq '.' || $file eq '..' || ! -d "$out_dir/$file");
		my $config = "$out_dir/$file/config.txt";
		my $processLog = "$out_dir/$file/process_current.log";
		$cnt++;
		if (-r "$config"){
			$list->{$cnt}->{TIME} =  strftime "%F %X",localtime((stat("$processLog"))[9]); 
			$list->{$cnt}->{STATUS} = "running" if $name2pid->{$file};
			if ( -r "$processLog"){
				open (my $fh, $processLog);
				while(<$fh>){
					if (/queued/){
						$list->{$cnt}->{STATUS} = "unstarted";
					}
					if (/All Done/){
						$list->{$cnt}->{STATUS} = "finished";
					}
					if (/failed/i){
						$list->{$cnt}->{STATUS} = "failed";
					}
					if (/interrupted/){
						$list->{$cnt}->{STATUS} = "interrupted";
					}
				}
				close $fh;
			}
			my ($projname, $projid, $projcode,$projCPU,$runhost);
			# if the system change from User management on to off. will need parse the project name from config file
			#  using grep is slow than open file and regrex
			open (my $fh, $config);
			while(<$fh>){
				chomp;
				$runhost=$1 if (/projrunhost=(.*)/);
				$projname=$1 if (/projname=(.*)/);
				$projid=$1 if (/projid=(.*)/);
				$projcode=$1 if (/projcode=(.*)/);
				$projCPU=$1 if (/cpu=(.*)/);
			}
			close $fh;
			
			chomp $projname;
			(my $domain = $runhost) =~ s/https?:\/\///;
			$list->{$cnt}->{NAME} = $file ;
			$list->{$cnt}->{PROJNAME} = $projname;
			$list->{$cnt}->{PROJCODE} = $projcode;
			$list->{$cnt}->{PROJDOMAIN} = $domain;
			$list->{$cnt}->{PROJDIR} = "$out_dir/$file";
			$list->{$cnt}->{CPU} = $projCPU;
		}

	}
	closedir(BIN);
}

sub availableToRun {
	my ($num_cpu, $cpu_been_used) = @_;
	my $count=0;
	return 0 if $cpu_been_used + $num_cpu >= $sys->{edgeui_tol_cpu};
	if( $sys->{edgeui_auto_queue} && $sys->{edgeui_tol_cpu} ){
		foreach my $i ( keys %$list ){
			if( $list->{$i}->{STATUS} eq "running" ){
				$count++;
				$cpu_been_used += $list->{$i}->{CPU};
				return 0 if $cpu_been_used + $num_cpu >= $sys->{edgeui_tol_cpu};
				return 0 if ($count >= $max_num_jobs);
			}
		}
		return 1;
	}
	else{
		return 0;
	}
}

sub getSystemUsage {
	my $mem = `vmstat -s | awk  '\$0 ~/total memory/ {total=\$1 } \$0 ~/free memory/ {free=\$1} \$0 ~/buffer memory/ {buffer=\$1} \$0 ~/cache/ {cache=\$1} END{print (total-free-buffer-cache)/total*100}'`;
	my $cpu = `top -bn1 | grep load | awk '{printf "%.1f", \$(NF-2)}'`;
	my $disk = `df -h $out_dir | tail -1 | awk '{print \$5}'`;
	$disk= `df -h $out_dir | tail -1 | awk '{print \$4}'` if ($disk !~ /\%/);
	$cpu = $cpu/$sys->{edgeui_tol_cpu}*100;
	$disk =~ s/\%//;
	if( $mem || $cpu || $disk ){
		$mem = sprintf "%.1f", $mem;
		$cpu = sprintf "%.1f", $cpu;
		$disk = sprintf "%.1f", $disk;
		return ($mem,$cpu,$disk);
	}
	else{
		return (0,0,0);
	}
}

sub checkProjVital {
	my $ps = `ps aux | grep run[P]`;
	my $vital;
	my $name2pid;
	my @line = split "\n", $ps;
	foreach my $line ( @line ){
		my ( $pid, $proj, $numcpu) = $line =~ /^\S+\s+(\d+) .*\/(\S+) -cpu (\d+)/;
		$vital->{$pid}->{PROJ} = $proj;
		$vital->{$pid}->{CPU} = $numcpu;
		$name2pid->{$proj} = $pid;
	}
	return ($vital,$name2pid);
}
