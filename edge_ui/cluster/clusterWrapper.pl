#!/usr/bin/perl
use XML::Simple;
use strict;

my $DEBUG = 0; 

my $BIN="/opt/uge/bin/lx-amd64";
$ENV{SGE_ROOT} ="/opt/uge";
$ENV{SGE_CELL} ="seqprod";
$ENV{SGE_CLUSTER_NAME} ="seqclust";

sub clusterSubmitJob {#submit job
    my ($in, $opt) = @_;
    my $error;

    if(! -e $in) {
        $error = "Job script file $in does not exist.\n";
	return ('', $error);
    }

    #submit job
    print "$BIN/qsub $opt $in\n" if $DEBUG;
    my $job = `$BIN/qsub $opt $in`;
    print "$job\n" if $DEBUG;

    my $job_id;
    #uge format
    #Your job 15504 ("Sleeper") has been submitted
    if($job =~ /Your job (\d+).*/){
        $job_id = $1;
    }

    unless ($job_id) {
        $error = "Queuing system error.\n";
    }
    
    return ($job_id, $error);

} 

sub clusterGetJobs {#get all jobs
    my $error;
    my $stats = `$BIN/qstat -xml`;
   	my $joblist;

	if($stats !~ /<queue_info>/) {
		$error = "Queuing system error.\n";
    } else {
		$joblist = XMLin($stats);
	}
	return ($joblist,$error);
} 

sub checkProjVital_cluster {
	my $cluster_job_prefix = shift;
	my $vital;
	my $name2pid;
	my ($joblist, $error) = clusterGetJobs();
	my ($pid, $proj, $stat, $jobName,$numcpu);
	if(!$error) {
		foreach my $job (@{$joblist->{queue_info}->{job_list}}) {
			$pid = $job->{JB_job_number};
   			$jobName = $job->{JB_name};
   			$stat = $job->{state}[1];
			if($stat eq 'Eqw') {
				$stat = "Cluster Job Error.";
			} 

			$numcpu = $job->{slots};

            if($jobName =~ /^$cluster_job_prefix(.*)/) {
				$proj = $1;
				$vital->{$pid}->{PROJ} = $proj;
				$vital->{$pid}->{CPU} = $numcpu;
				$vital->{$pid}->{STAT} = $stat;
				$name2pid->{$proj} = $pid;
			}
		}
	} 
	return ($vital,$name2pid,$error);
}

sub clusterDeleteJob {
    my $job_id = shift;
    my $status  = `$BIN/qdel $job_id 2>&1`;

	# making sure to give enough time for processing deletion request
 	sleep(3);
   	print "$status\n" if $DEBUG;

	# checkign success of deletion
	if($status =~ /has registered (.*) for deletion/ || $status =~ /job (.*) does not exist/) {
	    return 0;#success
	} else {
        return 1;#failure
    }
}

1;
