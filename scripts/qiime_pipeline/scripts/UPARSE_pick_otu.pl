#!/usr/bin/env perl
# chienchi at lanl.gov
# 20150910

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use POSIX qw(strftime);

my $input;
my $outDir="otus";
my $threads=4;
my $usearch_path="usearch8";
my $dryRun;
my $similarity=0.97;
my $min_otu_size=2;
my $suppress_align_and_tree;

GetOptions(
        "input|i=s"	=> \$input,
        "output|o=s"	=> \$outDir,
	"threads|t=i"      => \$threads,
	"usearch_path=s" => \$usearch_path,
	"similarity=f"	=> \$similarity,
	"min_otu_size=i" => \$min_otu_size,
	"dryRun"	=> \$dryRun,
        "help|?"	=> sub{Usage()}
);

sub Usage{
print <<"END";
Usage: perl $0 -i <demultiplexed_fastq>

	-input          demultiplexed fastq from split_libraries_fastq.py --store_demultiplexed_fastq
	-output         Output Directory [default: otus]
  Options:
	-similarity     sequence similarity threshold for mapping reads back to OTUs [default: 0.97]

	-min_otu_size   the minimum otu size (in number of sequences) to retain the otu [default: 2]

	-threads        Number of CPU threads

	-usearch_path   /path/to/usearch if the usearch (version >7) installed in a specific path.

	-suppress_align_and_tree  <bool>  skip the sequence alignment and tree-building steps

	-dryRun		<bool> print commands only

END
exit;
}

$ENV{OMP_NUM_THREADS}=$threads;
$ENV{OMP_THREAD_LIMIT}=$threads;

&Usage unless ($input);

my $outDirTmp="$outDir/tmp";
`mkdir -p $outDirTmp`;

#my $log= "$outDir/UPARSE.log";
my $cmd;
print "Print Command Only\n" if ($dryRun);

# remove low quality reads (trimming not required for paired-end data)
$cmd = "$usearch_path -fastq_filter $input -fastaout $outDirTmp/seqs.filtered.fasta -fastq_maxee 0.5 -threads $threads";
&process_cmd($cmd,"Remove low quality reads") unless (-s "$outDirTmp/seqs.filtered.fasta");

# dereplicate seqs
$cmd = "$usearch_path -derep_fulllength $outDirTmp/seqs.filtered.fasta  -fastaout $outDirTmp/seqs.filtered.derep.fasta -sizeout -threads $threads";
&process_cmd($cmd,"Dereplicate seqs") unless ( -s "$outDirTmp/seqs.filtered.derep.fasta");

# filter singletons
$cmd = "$usearch_path -sortbysize $outDirTmp/seqs.filtered.derep.fasta -minsize 2 -fastaout $outDirTmp/seqs.filtered.derep.mc2.fasta ";
&process_cmd($cmd,"Filter singletons") unless ( -s "$outDirTmp/seqs.filtered.derep.mc2.fasta");

# cluster OTUs
# if vsearch is used: --cluster_smallmem .fasta --centroids .fasta
$cmd = "$usearch_path -cluster_otus $outDirTmp/seqs.filtered.derep.mc2.fasta -otus $outDirTmp/seqs.filtered.derep.mc2.repset.fasta";
&process_cmd($cmd,"Cluster OTUs") unless ( -s "$outDirTmp/seqs.filtered.derep.mc2.repset.fasta");

# reference chimera check
$cmd = "$usearch_path -uchime_ref $outDirTmp/seqs.filtered.derep.mc2.repset.fasta -db $RealBin/../data/RDPClassifier_16S_trainsetNo14/trainset14_032015.fasta -strand plus -nonchimeras $outDirTmp/seqs.filtered.derep.mc2.repset.nochimeras.fasta -threads $threads " ;
&process_cmd($cmd,"Chimera check") unless ( -s "$outDirTmp/seqs.filtered.derep.mc2.repset.nochimeras.fasta");

# label OTUs using UPARSE python script
$cmd = "python $RealBin/drive5_py/fasta_number.py $outDirTmp/seqs.filtered.derep.mc2.repset.nochimeras.fasta OTU_ > $outDir/rep_set.fna ";
&process_cmd($cmd,"Label OTUs") unless ( -s "$outDir/rep_set.fna");

# map the _original_ quality filtered reads back to OTUs
$cmd = "$usearch_path -usearch_global $outDirTmp/seqs.filtered.fasta -db $outDir/rep_set.fna -strand plus -id $similarity -uc $outDirTmp/otu.map.uc -threads $threads ";
&process_cmd($cmd,"Map filtered reads back to OTUs") unless ( -s "$outDirTmp/otu.map.uc");


# make OTU table.
$cmd = "python $RealBin/drive5_py/uc2otutab_mod.py $outDirTmp/otu.map.uc > $outDirTmp/otu-table.txt ";
&process_cmd($cmd,"Make OTU table") unless ( -s "$outDirTmp/otu-table.txt");

# convert to biom
$cmd = "biom convert --table-type=\"OTU table\" -i $outDirTmp/otu-table.txt -o $outDir/otu-table.biom --to-hdf5";
&process_cmd($cmd,"Convert to biom") unless ( -s "$outDir/otu-table.biom");

# assign taxonomy 
$cmd = "parallel_assign_taxonomy_uclust.py -T --jobs_to_start $threads -i $outDir/rep_set.fna -o $outDir/uclust_assigned_taxonomy ";
&process_cmd($cmd,"Assign taxonomy ") unless ( -s "$outDir/uclust_assigned_taxonomy/rep_set_tax_assignments.txt");

# Add taxa to OTU table command 
$cmd = "biom add-metadata -i $outDir/otu-table.biom --observation-metadata-fp $outDir/uclust_assigned_taxonomy/rep_set_tax_assignments.txt -o $outDir/otu-table_w_tax.biom --sc-separated taxonomy --observation-header OTUID,taxonomy ";
&process_cmd($cmd,"Add taxa to OTU table") unless ( -s "$outDir/otu-table_w_tax.biom");

if ($suppress_align_and_tree){

# filter otu by miniumu observation number
$cmd = "filter_otus_from_otu_table.py -i $outDir/otu-table_w_tax.biom -o $outDir/otu_table_mc${min_otu_size}_w_tax.biom -n $min_otu_size "; 
&process_cmd($cmd,"Filter OTU by miniumu observation number") unless ( -s "$outDir/otu_table_mc${min_otu_size}_w_tax.biom");

}else{

# Align sequences command 
$cmd = "parallel_align_seqs_pynast.py -i $outDir/rep_set.fna -o $outDir/pynast_aligned_seqs -T --jobs_to_start $threads ";
&process_cmd($cmd,"Align sequences") unless (-s "$outDir/pynast_aligned_seqs/rep_set_aligned.fasta");

# Filter alignment command
$cmd = "filter_alignment.py -o $outDir/pynast_aligned_seqs -i $outDir/pynast_aligned_seqs/rep_set_aligned.fasta ";
&process_cmd($cmd,"Filter alignment") unless (-s "$outDir/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta");

# Build phylogenetic tree command 
$cmd = "make_phylogeny.py -i $outDir/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta -o $outDir/rep_set.tre";
&process_cmd($cmd,"Build phylogenetic tree") unless (-s "$outDir/rep_set.tre");

# filter otu by miniumu observation number and fail to align otus
$cmd = "filter_otus_from_otu_table.py -i $outDir/otu-table_w_tax.biom -e $outDir/pynast_aligned_seqs/rep_set_failures.fasta -o $outDir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom -n $min_otu_size "; 
&process_cmd($cmd,"Filter OTU by miniumu observation number and fail to align OTUs") unless ( -s "$outDir/otu_table_mc${min_otu_size}_w_tax_no_pynast_failures.biom");

}

sub process_cmd {
    my ($cmd, $msg, $dieCatch) = @_;

    if ($msg) {
        print "\n\n";
        print "###########################################################################\n";
        print &getTimeString."  $msg\n";
        print "###########################################################################\n";
    }
    
    print "CMD: $cmd\n";
	
    my $time_start = time();
    my $ret = system($cmd) if (!$dryRun);
    my $time_end = time();

    if ($ret) {
        if ($dieCatch)
        {
            print $ret ;
        }else
        {
            die "Error, CMD: $cmd died with ret $ret";
        }
    }

    my $number_minutes = sprintf("%.1f", ($time_end - $time_start) / 60);
    
    print "TIME: $number_minutes min.\n";
    
    return $ret;
}

sub getTimeString
{
    my $now_string = strftime "[%Y %b %e %H:%M:%S]", localtime;
    #$now_string =~ s/\s//g;
    return $now_string;
}
