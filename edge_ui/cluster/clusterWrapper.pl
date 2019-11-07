#!/usr/bin/perl
use XML::Simple;
use strict;

my $DEBUG = 0; 

sub LoadSGEenv {
	my $sys = shift;
	if ($sys->{sge_bin}){
		$ENV{PATH} = "$ENV{PATH}:$sys->{sge_bin}";
	}
	if ($sys->{slurm_bin}){
		$ENV{PATH} = "$ENV{PATH}:$sys->{slurm_bin}";
	}
        $ENV{SGE_ROOT} = "$sys->{sge_root}";
        $ENV{SGE_CELL} = "$sys->{sge_cell}";
        $ENV{SGE_CLUSTER_NAME} = "$sys->{sge_cluster_name}";
}

sub clusterSubmitJob {#submit job
    my ($scheduler, $in, $opt) = @_;
    
    my ($error, $job);
  
    my $submit_cmd = "qsub";
    $submit_cmd = "sbatch" if $scheduler=~ /slurm/;

    if(! -e $in) {
        $error = "Job script file $in does not exist.\n";
	return ('', $error);
    }

    #submit job
    print "$submit_cmd $opt $in\n" if $DEBUG;
    $job = `$submit_cmd $opt $in`;
    print "$job\n" if $DEBUG;
    
    my $job_id;
    my $submit_info_prefix = "Your job ";
    $submit_info_prefix = "Submitted batch job " if $scheduler=~ /slurm/;
    #uge format
    #Your job 15504 ("Sleeper") has been submitted
    if($job =~ /$submit_info_prefix(\d+)/){
        $job_id = $1;
    }

    unless ($job_id) {
        $error = "Queuing system error.\n";
    }
    
    return ($job_id, $error);

} 

sub clusterGetJobs {#get all jobs
    my $scheduler=shift;
    my $error;
    my $stats;
    if ($scheduler=~/slurm/){
        $stats = `squeue -h -a --format="%.18i %.9P %.32j %.10u %.2c %.2t %.10M %.6D %R" | grep EDGE` if $scheduler=~ /slurm/;
    }elsif($scheduler=~/uge|sge/){
        $stats = `qstat -xml`;
    }
    my $joblist;
     
    if ($scheduler=~ /slurm/){
    # JOBID PARTITION                             NAME     USER MI ST       TIME  NODES NODELIST(REASON)
    #    134      defq              EDGE_workshop_21998   andylo  1 PD       0:00      1 (Resources)
        foreach my $line (split/\n/,$stats){
            $line =~ s/^\s+|\s+$//;
            my ($jobid, $partition,$jobName, $user, $numcpu, $state,$time, $numNodes , $reason) = split/\s+/,$line;
            $joblist->{$jobid}->{partition} = $partition;
            $joblist->{$jobid}->{jobName}=$jobName;
            $joblist->{$jobid}->{user}=$user;
            $joblist->{$jobid}->{numcpu}=$numcpu;
            $joblist->{$jobid}->{state}=$state;
            $joblist->{$jobid}->{time}=$time;
            $joblist->{$jobid}->{numNodes}=$numNodes;
            $joblist->{$jobid}->{reason}=$reason;
        }
    }elsif($scheduler=~ /uge|sge/){
        if($stats !~ /<queue_info>/) {
            $error = "Queuing system error.\n";
        } else {
            $joblist = XMLin($stats);
        }
    }
    return ($joblist,$error);
} 

sub checkProjVital_cluster {
        my $scheduler=shift;
	my $cluster_job_prefix = shift;
	my $vital;
	my $name2pid;
	my ($joblist, $error) = clusterGetJobs($scheduler);
	my ($pid, $proj, $stat, $jobName,$numcpu);
	if(!$error) {
		my @jobs;
		if ($scheduler =~ /slurm/){
			@jobs= keys %{$joblist};
			foreach my $pid (@jobs) {
				$jobName=$joblist->{$pid}->{jobName};
				$stat=$joblist->{$pid}->{state};
				$numcpu=$joblist->{$pid}->{numcpu};
				if($stat eq 'F' || $stat eq 'NF') {
					$stat = "Cluster Job Error.";
				}
				
				if($jobName =~ /^$cluster_job_prefix(.*)/) {
					$proj = $1;
					$vital->{$pid}->{PROJ} = $proj;
					$vital->{$pid}->{CPU} = $numcpu;
					$vital->{$pid}->{STAT} = $stat;
					$name2pid->{$proj} = $pid;
				}	
			}
		}elsif($scheduler =~ /uge|sge/){
			foreach my $queue_stat ("queue_info","job_info"){
				# queue_info is running jobs, job_info is for pending jobs
				# if only one running, it is stored in HASH.
				if (ref($joblist->{$queue_stat}->{job_list}) eq "HASH"){
					@jobs = ($joblist->{$queue_stat}->{job_list});
				}elsif(ref($joblist->{$queue_stat}->{job_list}) eq "ARRAY"){
					@jobs = @{$joblist->{$queue_stat}->{job_list}};
				}
				foreach my $job (@jobs) {
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
		}

	}
	return ($vital,$name2pid,$error);
}

sub clusterDeleteJob {
    my $scheduler = shift;
    my $job_id = shift;
    my $job_delete_cmd="qdel";
    $job_delete_cmd = "scancel" if $scheduler=~ /slurm/;
    my $status  = `$job_delete_cmd $job_id 2>&1`;

	# making sure to give enough time for processing deletion request
 	sleep(3);
   	print "$status\n" if $DEBUG;

    if ($scheduler=~ /slurm/){
        return ($status =~ /error/)? 1 : 0;
    }

	# checkign success of deletion
	if($status =~ /has registered (.*) for deletion/ || $status =~ /job (.*) does not exist/) {
	    return 0;#success
	} else {
        return 1;#failure
    }
}
1;
