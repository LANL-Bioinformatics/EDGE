#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($RealBin);
use File::Basename;
use lib "$RealBin/../../lib";
use HTML::Template;
use POSIX qw{strftime};
use JSON;

my $out_dir       = $ARGV[0];
my $html_outfile  = $ARGV[1];
my $mode          = $ARGV[2];
$mode         ||= $html_outfile ? "web" : "";
$html_outfile ||= "$out_dir/HTML_Report/report.html";
my @out_dir_parts = split('/', $out_dir);
my $projname = $out_dir_parts[-1];
my $sysconfig    = "$RealBin/../../edge_ui/sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);

## Instantiate the variables
my $getting_paired=0;
my $vars;
my $refname;
my $NUM_READS_FOR_DOWNSTREAM = 0;

&check_analysis();

eval {
	&pull_summary();

	## sample metadata
	&pull_sampleMetadata();
	#END

	&pull_fastqCount();
#	print STDOUT "fastqcount\n";
	&pull_qc();
#	print STDOUT "qc\n";
	&pull_assy();
#	print STDOUT "assay\n";
	&pull_anno();
#	print STDOUT "anno\n";
	&pull_referenceName();
#	print STDOUT "reference name\n";
	&pull_readmapping_contig();
#	print STDOUT "read mapping contig\n";
	&pull_readmapping_ref();
#	print STDOUT "read mapping ref\n";
	&pull_contigmapping();
#	print STDOUT "contig mapping\n";
	&pull_host_rev();
#	print STDOUT "host rev\n";
	&pull_taxa();
#	print STDOUT "taxa\n";
	&pull_contig_taxa();
#	print STDOUT "contig taxa\n";
	&pull_snp_phylogeny();
#	print STDOUT "snp phylogeny\n";
	&pull_specialty_gene_profiling();
#	print STDOUT "specialty genes\n";
	&pull_pcr_contig_valid();
#	print STDOUT "pcr contig valid\n";
	&pull_pcr_ref_valid();
#	print STDOUT "pcr ref valid\n";
	&pull_pcr_design();
#	print STDOUT "pcr design\n";
	&pull_blast();
#	print STDOUT "blast\n";
	&pull_sra_download();
#	print STDOUT "sra download\n";
	&prep_jbrowse_link();
#	print STDOUT "jbrowse link\n";
	&checkImageReady();
#	print STDOUT "check image ready\n";
	&checkProjTarFile();
#	print STDOUT "check proj tar file\n";
};
#print STDOUT "start output html \n";
output_html();
#print STDOUT "end output html\n";
sub pull_referenceName {
	if( -e "$out_dir/Reference/reference.fasta" ){
		#try to parse ref name
		my @fasta_header =`grep "^>" $out_dir/Reference/reference.fasta`;
		foreach (@fasta_header){
			chomp $_;
			$refname->{$1}=$2 if ($_ =~ /^>(\S+)\s+(.+[a-zA-Z0-9])[^a-zA-Z0-9]?$/ );
		}
	}
}

sub checkProjTarFile{
	my $proj_realname = $vars->{PROJNAME};
	my $projID = $vars->{PROJID} || $projname;
	$vars->{OUT_GET_DOWNLOAD} = 0 if $vars->{PROJSTATUS} ne "Complete";
	my $tarFile1 = "$out_dir/$proj_realname.tgz";
	my $tarFile2 = "$out_dir/${proj_realname}_$projID.tgz";
	if (-e "$out_dir/.tarfinished" and (-e $tarFile1 or -e $tarFile2)){
		$vars->{OUT_DOWNLOAD} = 1;
		$vars->{OUT_GET_DOWNLOAD} = 0;
		$vars->{PROJTARFILE} = ( -e $tarFile1)? $tarFile1 : $tarFile2; 
	}
}

sub checkImageReady {
	my $log = "$out_dir/process.log";
	my $tmp = `grep "Converting pdf to png" $log`;
	if( -e "$out_dir/final_report.pdf" && $tmp ){
		$vars->{OUT_IMAGE_READY} = 1;
	}

	#for backfard compatible
	if( -e "$out_dir/HTML_Report/images/Assembly_CovDepth_vs_Len.png" ){
		$vars->{OUT_IMAGE_READY_ASSEMBLY_VALID} = 1;
	}
}

sub _reformat_val {
	my $content = shift;
	return unless defined $content;

	my $outdir_rgx = $out_dir;
	$outdir_rgx =~ s/\//\\\//g;

	if( $content =~ /^\d{4,}$/ ){
		$content =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	}
	elsif( $content =~ /^(\d{4,}) \(([\d\.]+ \%)\)$/ ){
		my ($num,$pct) = ($1,$2);
		$num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
		$content = "$num ($pct)";
	}
	elsif( $content =~ /^(\d+\.\d{2})\d+$/ ){
		$content = $1;
	}
	elsif( $content =~ /$outdir_rgx/ ){
		$content =~ s/$outdir_rgx/\.\./ unless $mode eq "web";
	}
	return $content;
}

sub output_html {
	$vars->{OUTPUTDIR}  = $out_dir;
	#$vars->{PDFREPORT}  = "$out_dir/final_report.pdf";
	#$vars->{PROCESSLOG} = "$out_dir/process.log";
	#$vars->{ERRORLOG}   = "$out_dir/error.log";
	$vars->{PROJNAME} ||= $projname;
	$vars->{PROJID} ||= $projname;
	#reformat number with thousand separator
	foreach my $var ( keys %$vars ){
		next if ($var eq "PROJID");
		#sample metadata
		next if ($var =~ /SMD_/);
		#END sample metadata
		if( ref($vars->{$var}) eq 'ARRAY' ){
			for( my $i=0; $i < scalar(@{$vars->{$var}}); $i++ ){
				if( ref($vars->{$var}[$i]) eq 'HASH' ){
					foreach my $idx ( keys %{$vars->{$var}[$i]} ){
						$vars->{$var}[$i]->{$idx} = &_reformat_val($vars->{$var}[$i]->{$idx});
					}
				}
				else{
					$vars->{$var}[$i] = &_reformat_val($vars->{$var}[$i]);
				}
			}
		}
		else{
			$vars->{$var} = &_reformat_val($vars->{$var});
		}
	}

	my $template = HTML::Template->new(filename => "$RealBin/edge_html_report.tmpl",
		                               strict => 0,
								       die_on_bad_params => 0);
	$template->param(%$vars);

	system("mkdir -p $out_dir/"."HTML_Report");	
	
	system("cp -r $ENV{'EDGE_HOME'}/edge_ui/css $out_dir/HTML_Report/")        if $mode ne "web"; 
	system("cp -r $ENV{'EDGE_HOME'}/edge_ui/images $out_dir/HTML_Report/")     if $mode ne "web";  
	system("cp -r $ENV{'EDGE_HOME'}/edge_ui/javascript $out_dir/HTML_Report/") if $mode ne "web";  

	open(my $htmlfh, ">$html_outfile") or die $!;
	print $htmlfh $template->output();
	close ($htmlfh);
}

sub check_analysis {
	$vars->{OUT_GET_DOWNLOAD} = 1 if $mode eq "web";
	$vars->{OUT_QC_SW}    = 1 if -e "$out_dir/QcReads/runQC.finished";
	$vars->{OUT_HR_SW}    = 1 if -e "$out_dir/HostRemoval/hostclean.stats.txt";
	$vars->{OUT_AS_SW}    = 1 if (glob "$out_dir/AssemblyBasedAnalysis/*.finished");
	$vars->{OUT_AS_PC_SW} = 1 if ( -e "$out_dir/AssemblyBasedAnalysis/processProvideContigs.finished");
	#$vars->{OUT_AS_SW}    = 1 if -e "$out_dir/AssemblyBasedAnalysis/runIdbaAssembly.finished";
	$vars->{OUT_AN_SW}    = 1 if -e "$out_dir/AssemblyBasedAnalysis/Annotation/runAnnotation.finished";
	$vars->{OUT_RA_SW}    = 1 if -e "$out_dir/ReadsBasedAnalysis/readsMappingToRef/runReadsToGenome.finished";
	$vars->{OUT_UM_CP_SW} = 1 if -e "$out_dir/ReadsBasedAnalysis/UnmappedReads/Taxonomy/taxonomyAssignment.finished";
	$vars->{OUT_AL_CP_SW} = 1 if -e "$out_dir/ReadsBasedAnalysis/Taxonomy/taxonomyAssignment.finished";
	$vars->{OUT_RCP_SW}   = 1 if $vars->{OUT_UM_CP_SW} || $vars->{OUT_AL_CP_SW}; 
	$vars->{OUT_CCP_SW}   = 1 if -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/ContigsTaxonomy.finished"; 
	$vars->{OUT_CP_SW}    = 1 if $vars->{OUT_UM_CP_SW} || $vars->{OUT_AL_CP_SW} || $vars->{OUT_CCP_SW}; 
	$vars->{OUT_PA_SW}    = 1 if -e "$out_dir/AssayCheck/";
	$vars->{OUT_SP_SW}    = 1 if -e "$out_dir/SNP_Phylogeny/";
	$vars->{OUT_PA_V_SW}  = 1 if -e "$out_dir/AssayCheck/pcrValidation.pdf";
	$vars->{OUT_PA_D_SW}  = 1 if -e "$out_dir/AssayCheck/PCR.design.primers.txt";
	$vars->{OUT_AB_SW}    = 1 if ( -e "$out_dir/AssemblyBasedAnalysis/Blast/ContigsBlast.finished" and -s "$out_dir/AssemblyBasedAnalysis/Blast/ContigsForBlast");
	$vars->{OUT_qiime_SW}   = 1 if ( -e "$out_dir/QiimeAnalysis/runQiimeAnalysis.finished");
	$vars->{OUT_SG_SW}    = 1 if (-e "$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished" or -e "$out_dir/ReadsBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished" );
	$vars->{OUT_OSG_SW}   = 1 if ( -e "$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished");
	$vars->{OUT_RSG_SW}   = 1 if ( -e "$out_dir/ReadsBasedAnalysis/SpecialtyGenes/runSpecialtyGenesProfiling.finished");
}

sub prep_jbrowse_link {
	$vars->{JB_CTG_ANNO}     = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA";
	$vars->{JB_CTG_ANNO_BAM} = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA%2CBAM"; 
	$vars->{JB_CTG_ANNO_PCR} = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA%2CPCR_V%2CPCR"; 
	$vars->{JB_CTG_ANNO_SG}  = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fctg_tracks&tracks=DNA%2CCDS%2CVFmarker%2CARmarker"; 
	$vars->{JB_REF_ANNO}     = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA"; 
	$vars->{JB_REF_ANNO_BAM} = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA%2CBAM"; 
	$vars->{JB_REF_ANNO_CTG} = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA%2CCTG2REF"; 
	$vars->{JB_REF_ANNO_VCF} = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CRRNA%2CTRNA%2CVCF";
	$vars->{JB_REF_CDS_CTG}  = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CCTG2REF"; 
	$vars->{JB_REF_CDS_BAM}  = "JBrowse/?data=data%2F$projname%2FJBrowse%2Fref_tracks&tracks=DNA%2CCDS%2CBAM"; 
}


## sample metadata
sub pull_sampleMetadata {
	my $metadata = "$out_dir/metadata_sample.txt";
	if($sys->{edge_sample_metadata} && -s $metadata) {
		$vars->{OUT_SAMPLE_METADATA}   = 1;
	}
	if(-e $metadata) {
        	open CONF, $metadata or die "Can't open $metadata $!";
        	while(<CONF>){
      			chomp;
      			next if(/^#/);
           		if ( /(.*)=(.*)/ ){
             			$vars->{SMD_STUDY_TITLE} =$2 if ($1 eq "study_title");
             			$vars->{SMD_STUDY_TYPE} =$2 if ($1 eq "study_type");
             			$vars->{SMD_NAME} =$2 if ($1 eq "sample_name");
             			$vars->{SMD_TYPE} =$2 if ($1 eq "sample_type");
             			$vars->{OUT_SMD_GENDER} =1 if ($1 eq "gender");
             			$vars->{SMD_GENDER} =$2 if ($1 eq "gender");
             			$vars->{OUT_SMD_AGE} =1 if ($1 eq "age");
             			$vars->{SMD_AGE} =$2 if ($1 eq "age");
             			$vars->{OUT_SMD_HOST} =1 if ($1 eq "host");
             			$vars->{SMD_HOST} =$2 if ($1 eq "host");
             			$vars->{OUT_SMD_HOST_CONDITION} =1 if ($1 eq "host_condition");
             			$vars->{SMD_HOST_CONDITION} =$2 if ($1 eq "host_condition");
             			$vars->{SMD_SOURCE} =$2 if ($1 eq "isolation_source");
             			$vars->{SMD_COLLECTION_DATE} =$2 if ($1 eq "collection_date");
             			$vars->{SMD_LOCATION} =$2 if ($1 eq "location");
             			$vars->{SMD_CITY} =$2 if ($1 eq "city");
             			$vars->{SMD_STATE} =$2 if ($1 eq "state");
             			$vars->{SMD_COUNTRY} =$2 if ($1 eq "country");
             			$vars->{SMD_LAT} =$2 if ($1 eq "lat");
             			$vars->{SMD_LNG} =$2 if ($1 eq "lng");
             			$vars->{SMD_EXP_TITLE} =$2 if ($1 eq "experiment_title");
             			$vars->{SMD_SEQ_CENTER} =$2 if ($1 eq "sequencing_center");
             			$vars->{SMD_SEQUENCER} =$2 if ($1 eq "sequencer");
             			$vars->{SMD_SEQ_DATE} =$2 if ($1 eq "sequencing_date");
              		}
		}
        	close CONF;
		pull_other();
        	if($vars->{SMD_TYPE} eq "human") {
				pull_travels();
				pull_symptoms();
		}

	} else {
		$vars->{OUT_SAMPLE_NO_METADATA}   = 1;
	}
}

sub pull_other {
	my $other = "$out_dir/metadata_other.txt";
	return unless -e $other;
	open CONF, $other or die "Can't open $other $!";
        while(<CONF>){
      		chomp;
        	next if(/^#/);
     		if ( /(.*)=(.*)/ ){
			my $sym;
			$sym->{SMD_OTHER_F}=$1;
			$sym->{SMD_OTHER_V}=$2;
			push @{$vars->{LOOP_SMD_OTHER}}, $sym;
		}
	}
	close CONF;
}

sub pull_symptoms {
	my $symptoms = "$out_dir/metadata_symptoms.txt";
	return unless -e $symptoms;
	open CONF, $symptoms or die "Can't open $symptoms $!";
        while(<CONF>){
      		chomp;
        	next if(/^#/);
     		if ( /(.*)\t(.*)/ ){
			$vars->{OUT_SMD_SYMS} =1; 
			my $sym;
			$sym->{SMD_SYM_CAT}=$1;
			$sym->{SMD_SYM}=$2;
			push @{$vars->{LOOP_SMD_SYMS}}, $sym;
		}
	}
	close CONF;
}

sub pull_travels {
	my $travels = "$out_dir/metadata_travels.txt";
	return unless -e $travels;
	open CONF, $travels or die "Can't open $travels $!";
	
	my $from;
	my $to;
        while(<CONF>){
      		chomp;
        	next if(/^#/);
     		if ( /(.*)=(.*)/ ){
			$vars->{OUT_SMD_TVLS} =1; 
			if ($1 eq "travel-date-from") {
				$from = $2;
			} elsif ($1 eq "travel-date-to") {
				$to = $2;
			} elsif ($1 eq "travel-location") {
				my $sym;
				$sym->{SMD_TVL_DATE}="$from ~ $to";
				$sym->{SMD_TVL_LOC}=$2;
				push @{$vars->{LOOP_SMD_TVLS}}, $sym;
			}
		}
	}
	close CONF;
}

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
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
	return $sys;
}

##END sample metadata

sub pull_assy {
	my $err;
	$err = `grep failed $out_dir/AssemblyBasedAnalysis/assembly.log | sed -e 's/\$/<br>/'` if (-e "$out_dir/AssemblyBasedAnalysis/assembly.log");
	if ($err){
		chomp $err;
		$vars->{ASSYERR} = $err;
		return;
	}
	return unless -e "$out_dir/AssemblyBasedAnalysis/contigs_stats.txt";
	open(my $assyfh, "<", "$out_dir/AssemblyBasedAnalysis/contigs_stats.txt") or die $!;
        my $tmp = <$assyfh>;
	while(<$assyfh>){
                my @array = split /\s+/,$_;
                $vars->{ASSYNUMCONTIGS} = $array[1];
                $vars->{ASSYN50} = $array[2];
                $vars->{ASSYN90} = $array[3];
                $vars->{ASSYMAX} = $array[4];
                $vars->{ASSYMIN} = $array[5];
                $vars->{ASSYSIZE} = $array[8];
	#	if($_ =~ /.*contigs.fa\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)/) {
			## File	Contigs_number	N50	N90	Max	Min	Mean	Median	Total_bases	Top10_bases	Top20_bases	Top40_bases	Top100_bases	>100kb_bases	>50kb_bases	>25kb_bases	>10kb_bases	>5kb_bases	>3kb_bases	>2kb_bases	>1kb_bases
	#	}
	}
	close($assyfh);
}

sub pull_contigmapping {
	return unless -e "$out_dir/AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.log";
	open(my $cmfh, "<", "$out_dir/AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.log") or die $!;
	while(<$cmfh>) {
		# v1.1 use Total_reads
		if ($_ =~ /^Total_reads:\s+(\d+)/) 		  { $vars->{CMREADS} = $1; next; }
		# v2.0 correct to Total_contigs
		if ($_ =~ /^Total_contigs:\s+(\d+)/) 		  { $vars->{CMREADS} = $1; next; }
		if ($_ =~ /^Unused Contigs#:\s+(\d+)/)     { $vars->{CMUNUSED} = $1; next; }
		if ($_ =~ /^Avg_coverage_fold:\s(.*)\n/)  { $vars->{CMAVEFOLD} = $1; next; }
		if ($_ =~ /^Reference_Coverage:\s(.*)\n/) { $vars->{CMREFCOV} = $1; next; }
		if ($_ =~ /^Avg_Identity:\s(.*)\n/) 	{ $vars->{CMREFIDY} = $1; next; }
		if ($_ =~ /^Number of SNPs:\s(\d+)\n/) 	  { $vars->{CMSNPS} = $1; next; }
		if ($_ =~ /^Number of INDELs:\s(\d+)\n/)  { $vars->{CMINDELS} = $1; next; } 
	}
	close($cmfh);

	$vars->{CMMAPPED}    = $vars->{CMREADS} - $vars->{CMUNUSED};
	$vars->{CMMAPPEDPCT} = sprintf "%.2f", $vars->{CMMAPPED}/$vars->{CMREADS}*100;

	return unless -e "$out_dir/AssemblyBasedAnalysis/contigMappingToRef/contigsToRef_avg_coverage.table";
	my $ref_display_limit = 30;
	my ($tol_ref_number, $tol_ref_hashit) = (0,0);
	my $cnt=0;
	open(my $reffh, "<", "$out_dir/AssemblyBasedAnalysis/contigMappingToRef/contigsToRef_avg_coverage.table") or die $!;
	while(<$reffh>) {
		next if( /^ID/ );
		my @temp = split /\t/, $_;
		if( scalar @temp >= 5 ){
			my $refinfo;
			if ($cnt++<=$ref_display_limit){
				for my $i (0 .. $#temp) {
					my $idx = $i+1;
					$refinfo->{"CMREF$idx"}=$temp[$i];
				}
				$refinfo->{"CMREFNAME"}=$refname->{$temp[0]};
				push @{$vars->{LOOP_CMREF}}, $refinfo;
			}
			$tol_ref_hashit++ if $temp[3];
			$tol_ref_number++;
		}
	}
	close($reffh);
	$vars->{CMREFTOLREF}      = $tol_ref_number;
	$vars->{CMREFTOLREFHASHIT} = $tol_ref_hashit;
	$vars->{CMREFTABLENOTE} = "Only top $ref_display_limit results in terms of \"Base Recovery %\" are listed in the table." if $tol_ref_number > $ref_display_limit;
	## display unmapped contigs mapping to RefSeq
	$vars->{CMREF_UM} = 1 if -e "$out_dir/ReferenceBasedAnalysis/UnmappedContigs/log.txt";
	if ($vars->{CMREF_UM}){
		open (my $Log_fh, "<", "$out_dir/ReferenceBasedAnalysis/UnmappedContigs/log.txt") or die $!;
		while(<$Log_fh>){
			chomp;
			if ( $_ =~ /Total Contigs: (\d+) \((\d+) bp\); Classified Contigs: (\d+) \((\d+) bp\); Unclassified Contigs: (\d+) \((\d+) bp\);/){
				$vars->{CMREF_UM_NTC}  = $1;
				$vars->{CMREF_UM_NTCB} = $2;
				$vars->{CMREF_UM_NCC}  = $3;
				$vars->{CMREF_UM_NCCB} = $4;
				$vars->{CMREF_UM_NUC}  = $5;
				$vars->{CMREF_UM_NUCB} = $6;
			}
		}
		close $Log_fh;
		# unmapped portion of mapped contigs
		$vars->{CMREF_UM_NUPMC} = $vars->{CMREF_UM_NTC} - $vars->{CMUNUSED};
	}
}	

sub pull_anno {
	my $err;
	$err = `grep "No contigs" $out_dir/AssemblyBasedAnalysis/Annotation/Annotation.log | sed -e 's/\$/<br>/'` if (-e "$out_dir/AssemblyBasedAnalysis/Annotation/Annotation.log");
	if ($err){
		$vars->{ANOERR} = $err;
	}

	return unless -e "$out_dir/AssemblyBasedAnalysis/Annotation/$vars->{PROJNAME}.txt";
	open(my $anno, "<", "$out_dir/AssemblyBasedAnalysis/Annotation/$vars->{PROJNAME}.txt") or die $!;
	while(<$anno>){
		chomp;
		if( /^CDS: (\d+)/ ){
			$vars->{ANOCDS} = $1;
		}
		elsif( /^tRNA: (\d+)/ ){
			$vars->{ANOTRNA} = $1;
		}
		elsif( /^rRNA: (\d+)/ ){
			$vars->{ANORRNA} = $1;
		}
	}
	close $anno;
}

sub pull_qc {
	my $err;
	$err = `grep -i 'ERROR\\|invalid' $out_dir/QcReads/QC.log | sed -e 's/\$/<br>/'` if ( -e "$out_dir/QcReads/QC.log");
	if ($err){
		$out_dir =~ /.*(EDGE_output.*)/;
		$vars->{QCERR} = "$err"."Please check <a target='new_window' data-ajax='false' href='$1/QcReads/QC.log’>QC.log</a>.";
	}
	return unless -e "$out_dir/QcReads/QC.stats.txt";
	open(my $qcfh, "<", "$out_dir/QcReads/QC.stats.txt") or die $!;
	my $after=0;
	foreach(<$qcfh>) {
		chomp;
		if ($_ =~ /^Reads #:\s(\d+ \(\d+\.\d+\s\%\))/)              { $vars->{AFTERREADS} = $1; next; }
		if ($_ =~ /^Total bases:\s(\d+ \(\d+\.\d+\s\%\))/)          { $vars->{AFTERBASES} = $1; next; }
		if ($_ =~ /^Mean Reads Length: (.+)/)                       { $vars->{AFTERMRL} = $1; next; }
		if ($_ =~ /\s+Paired Reads #:\s(\d+ \(\d+\.\d+\s\%\))/)     { $vars->{AFTERPAIRED} = $1; next; }
		if ($_ =~ /\s+Paired total bases:\s(\d+ \(\d+\.\d+\s\%\))/) { $vars->{AFTERPTB} = $1; next; }
		if ($_ =~ /\s+Unpaired Reads #:\s(.*)/)                     { $vars->{AFTERUNPAIRED} = $1; next; }
		if ($_ =~ /\s+Unpaired total bases:\s(.*)/)                 { $vars->{AFTERUTB} = $1; next; }

		$after=1 if /After Trimming/;
		if($after){
			next;
		}
		if ($_ =~ /^Reads #:\s(\d+)/) 	{ $vars->{BEFOREREADS} = $1; next; }
		if ($_ =~ /^Total bases:\s(\d+)/) { $vars->{BEFOREBASES} = $1; next; }
		if ($_ =~ /^Reads Length:\s(.+)/) { $vars->{BEFOREMRL} = $1; next; }
	}
	close ($qcfh);
	$NUM_READS_FOR_DOWNSTREAM = $vars->{AFTERREADS};
}

sub pull_fastqCount {
	return unless -e "$out_dir/QcReads/fastqCount.txt";
	open(my $qcfh, "<", "$out_dir/QcReads/fastqCount.txt") or die $!;
	foreach(<$qcfh>) {
		chomp;
		my @temp = split /\t/, $_;
		if(@temp){
			$vars->{BEFOREREADS} += $temp[1];
			$vars->{BEFOREBASES} += $temp[2];
		}
	}
	$vars->{BEFOREMRL} = sprintf ('%.2f',$vars->{BEFOREBASES}/$vars->{BEFOREREADS});
	$NUM_READS_FOR_DOWNSTREAM = $vars->{BEFOREREADS};
	close ($qcfh);
}

sub pull_host_rev {
	return unless -e "$out_dir/HostRemoval/hostclean.stats.txt";
	open(my $hrfh, "<", "$out_dir/HostRemoval/hostclean.stats.txt") or die $!;
	while(<$hrfh>) {
		chomp;
		if ( /^(.+): (.+)$/   ){
			my $hr;
			$hr->{HRTITLE}=$1;
			$hr->{HRVALUE}=$2;
			$NUM_READS_FOR_DOWNSTREAM = $2 if $1 eq "Total non-host reads";
			push @{$vars->{LOOP_HR}}, $hr;
		}
	}
	close $hrfh;
}

sub sort_by_best_hit {
    my $transRef = shift;
    @$transRef = sort { $a->[8] cmp $b->[8] } @$transRef;
    return $transRef;
}

sub pull_specialty_gene_profiling {
	my $display_limit=40;
	my $proj_realname = $vars->{PROJNAME};
	if ($vars->{OUT_RSG_SW}){
		my $num_input_reads;
		($num_input_reads) = ($NUM_READS_FOR_DOWNSTREAM =~ m/(\d+)/);
		$vars->{SGREADTOTAL} = $num_input_reads;
		# Support old runs
		my ($ar_reads_output, $ar_reads_output_table, $vf_reads_output, $vf_reads_output_json, $vf_reads_output_sbhits, $vf_reads_output_table);
		if ( -e "$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_ShortBRED.txt") 
		{

			$ar_reads_output="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_ShortBRED.txt";
			$ar_reads_output_table="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_ShortBRED_table.txt";
			$vf_reads_output="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED.txt";
			$vf_reads_output_json="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED_table.json";
			$vf_reads_output_sbhits="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_SBhits.txt";
			$vf_reads_output_table="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED_table.txt";
		}
		else 
		{
			$ar_reads_output="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/AR_genes_ShortBRED.txt";
			$ar_reads_output_table="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/AR_genes_ShortBRED_table.txt";
			$vf_reads_output="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/VF_genes_ShortBRED.txt";
			$vf_reads_output_json="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/VF_genes_ShortBRED_table.json";
			$vf_reads_output_sbhits="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/VF_genes_SBhits.txt";
			$vf_reads_output_table="$out_dir/ReadsBasedAnalysis/SpecialtyGenes/VF_genes_ShortBRED_table.txt";
		}
		if ( -e $ar_reads_output) {
			my $ar_read_hit_num=`cat $ar_reads_output_table | wc -l`;
			chomp $ar_read_hit_num;
			$ar_read_hit_num --;
			$vars->{SGRARHIT} = ($ar_read_hit_num)? $ar_read_hit_num : 0;
		}
		if ( -e $ar_reads_output_table) {
			my $total_num=0;
			my $cnt=1;
			open (my $fh, $ar_reads_output_table) or die "Cannot read $ar_reads_output_table\n";
			while(<$fh>){
				next if (/^Family/);
				my $ar_info;
				#Family  Count   Hits    TotMarkerLength ID_Orig Merged.ID       Class   Source  Description
				my @temp = split /\t/, $_;
				if($cnt++<=$display_limit){
					$temp[0] = '<a href=\'http://ardb.cbcb.umd.edu/cgi/search.cgi?db=R&term='."$temp[0]"."'>$temp[0]</a>" if ($temp[-2] eq "ARDB");
					for my $i (0 .. $#temp) {
						my $idx = $i+1;
						$ar_info->{"SGRAR$idx"}=$temp[$i];
						
					}
					push @{$vars->{LOOP_SGRAR}}, $ar_info;
				}
				$total_num++;
			}
			close $fh;
			$vars->{SGRARNOTE} = "Only top $display_limit results in terms of RPKM are listed in the table." if ($total_num>$display_limit);
		}
		if ( -e $vf_reads_output) {
			#my $vf_read_hit_num=`awk 'BEGIN{a=0}{a=a+\$3}END{print a}' $vf_reads_output`;
			my $vf_read_hit_num=`awk 'BEGIN{count=0}\$2>0{count++}END{print count}'  $vf_reads_output`;
			chomp $vf_read_hit_num;
			$vars->{SGRVFHIT} = ($vf_read_hit_num)? $vf_read_hit_num : 0;
		}
		if ( -e $vf_reads_output_sbhits) {
			my $vf_read_count=`cut -f 1 $vf_reads_output_sbhits | sort | uniq | wc -l`;
			chomp $vf_read_count;
			$vars->{SGRVFCOUNT} = ($vf_read_count)? $vf_read_count : 0;
			my $vf_read_percent = ($vf_read_count/$num_input_reads)*100;
			$vars->{SGRVFPERC} = $vf_read_percent;
		}
		if (-e $vf_reads_output_json){
			my $json_text;{ local $/; open (my $jsonfh, "<", $vf_reads_output_json) or die "Cannot read json file\n"; $json_text = <$jsonfh>;close $jsonfh;}
			my $jsondata = decode_json($json_text);
			my ($genusCount, $classCount, $vfCount) = 0;
			for my $class (sort keys %$jsondata){
				my $classInfo;
				my $vfHash = $jsondata->{$class};
				#print STDOUT $vfHash." VFHash Class ".$class."\n";
				for my $vf (sort keys %$vfHash){
					my $vfInfo;
					my $detailsArray = $vfHash->{$vf};
					#print STDOUT $detailsArray." DetailsArray VF ".$vf."\n";
					for my $detailsHash (sort { $a->{vfgene} cmp $b->{vfgene} } @{$detailsArray}){
						my $detailInfo;
						$detailInfo->{VFGENE} = $detailsHash->{vfgene};
						my $vfnumber = $detailsHash->{vfnumber};
						$detailInfo->{VFNUMBER} = '<a href=\'http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID='."$vfnumber"."'>$vfnumber</a>";
						my $vfaccession = $detailsHash->{accession};
						$detailInfo->{VFACCESSION} = '<a href=\'http://www.ncbi.nlm.nih.gov/protein/'."$vfaccession"."'>$vfaccession</a>";
						$detailInfo->{VFTAXID} = $detailsHash->{taxid};
						$detailInfo->{VFSOURCE} = $detailsHash->{source};
						$detailInfo->{VFFOUNDIN} = $detailsHash->{foundin};
						$detailInfo->{VFGI} = $detailsHash->{gi};
						$detailInfo->{VFVFDBINFO} = $detailsHash->{vfdbinfo};
						$genusCount++;
						$classCount++;
						$vfCount++;
						push @{$vfInfo->{LOOP_DETAIL_RES}}, $detailInfo;
						#print STDOUT $detailInfo." DetailInfo\n";
						#print STDOUT $detailsHash." DetailsHash Details ".$details."\n";
						#for my $vfrow (keys $jsondata->{$genus}->{$class}->{$vf}->[$vfrowindex]){
							#print STDOUT $vfrow."\n";
						#}
					} # End Details
					$vfInfo->{VFVF_COUNT} = $vfCount;
					$vfInfo->{VFVF} = $vf;
					$vfCount = 0;
					push @{$classInfo->{LOOP_VF_RES}},$vfInfo;
				} # End VF
				$classInfo->{VFCLASS_COUNT} = $classCount;
				$classInfo->{VFCLASS} = $class;
				$classCount = 0;
				push @{$vars->{LOOP_SGRVF_CLASS_RES}}, $classInfo;					
			} # End Class
		}
		if ( -e $vf_reads_output_table) {
			my $total_num=0;
			my $cnt=1;

			$vars->{SGRVF_KRONA} = "$out_dir/ReadsBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED.krona.html";
			
			open (my $fh, $vf_reads_output_table) or die "Cannot read $vf_reads_output_table\n";
			while(<$fh>){
				next if (/^VFDB_Genus/);
				my $vf_info;
		
		#0Family  1Count   2Hits    3TotMarkerLength Property        4Source  5Gene Name       6Locus Tag       7Gene ID 8GI      9Organism        10Product Function        11Classification 
		#1VFGenus 2VFClass 3VF 4VFGene 5GeneraRepresented 6VFNumber 7GI 8OrganismGenus 9NewGeneAccession 10Score	
				my @temp = split /\t/, $_;
				if($cnt++<=$display_limit){
					$temp[5] = '<a href=\'http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID='."$temp[5]"."'>$temp[5]</a>";
					$temp[8] = '<a href=\'http://www.ncbi.nlm.nih.gov/protein/'."$temp[8]"."'>$temp[8]</a>";
					for my $i (0 .. $#temp) {
						my $idx = $i+1;
						$vf_info->{"SGRVF$idx"}=$temp[$i];
					}
					push @{$vars->{LOOP_SGRVF}}, $vf_info;
				}
				$total_num++;
			}
			close $fh;
			$vars->{SGRVFNOTE} = "Only the first $display_limit results displayed in the table." if ($total_num>$display_limit);
		}
	}
	if ($vars->{OUT_OSG_SW}){
		my ($ar_orf_output, $ar_orf_output_json, $ar_orf_output_table, $vf_orf_output, $vf_orf_output_json, $vf_orf_output_table);
		if ( -e "$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_rgi.json") 
		{
			$ar_orf_output="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_rgi.json";
			$ar_orf_output_json="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_rgi_table.json";
			$ar_orf_output_table="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_AR_genes_rgi.txt";
			$vf_orf_output="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_SBhits.txt";
			$vf_orf_output_json="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED_table.json";
			$vf_orf_output_table="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED_table.txt";
		}	
		else
		{
			$ar_orf_output="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/AR_genes_rgi.json";
			$ar_orf_output_json="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/AR_genes_rgi_table.json";
			$ar_orf_output_table="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/AR_genes_rgi.txt";
			$vf_orf_output="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/VF_genes_SBhits.txt";
			$vf_orf_output_json="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/VF_genes_ShortBRED_table.json";
			$vf_orf_output_table="$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/VF_genes_ShortBRED_table.txt";
		}
		if ( -e $ar_orf_output){
			#my $ar_orf_hit_num = `awk '{print \$2}' $ar_orf_output | sort | uniq | wc -l`;
			my $ar_orf_hit_num = `tail -n +2 $ar_orf_output_table | wc -l`;
			chomp $ar_orf_hit_num;
			$vars->{SGOARHIT} = ($ar_orf_hit_num)? $ar_orf_hit_num : 0;
		}
		if (-e $ar_orf_output_json && $vars->{SGOARHIT} > 0 ){
			my $json_text;{ local $/; open (my $jsonfh, "<", $ar_orf_output_json) or die "Cannot read json file\n"; $json_text = <$jsonfh>;close $jsonfh;}
			my $jsondata = decode_json($json_text);
			my $categoryCount = 0;
			for my $category (sort keys %$jsondata){
				my $categoryInfo;
				my $detailsArray = $jsondata->{$category};
					#print STDOUT $detailsArray." DetailsArray VF ".$vf."\n";
					for my $detailsHash (sort { $a->{bestHitName} cmp $b->{bestHitName} } @{$detailsArray}){
						my $detailInfo;
						#my $detailsHash = $detailsArray->[$details];
						$detailInfo->{ORFID} = $detailsHash->{orfID};
						my $bestHitARO = $detailsHash->{bestHitARO};
						my $bestHitAROName = $detailsHash->{bestHitName};
						$detailInfo->{BESTHIT} = '<a href=\'https://card.mcmaster.ca/aro/'."$bestHitARO"."'>$bestHitAROName</a>";
						$detailInfo->{SNPS} = $detailsHash->{snps};
						$detailInfo->{BITSCORE} = $detailsHash->{bestBitScore};
						$detailInfo->{PASSEVALUE} = $detailsHash->{passEvalue};
						$detailInfo->{EVALUE} = $detailsHash->{bestHitEvalue};
						$detailInfo->{CUTOFF} = $detailsHash->{cutOff};
						$detailInfo->{AROCATEGORIES} = $detailsHash->{bestHitCategories};
						$detailInfo->{OTHERHITS} = $detailsHash->{otherHits};
						$categoryCount++;
						push @{$categoryInfo->{LOOP_DETAIL_RES}}, $detailInfo;
					} # End Details
				$categoryInfo->{ARCATEGORY_COUNT} = $categoryCount;
				$categoryInfo->{ARCATEGORY} = $category;
				$categoryCount = 0;
				push @{$vars->{LOOP_SGOAR_CAT_RES}}, $categoryInfo;					
			} # End Class
		}
		if ( -e $ar_orf_output_table && $vars->{SGOARHIT} > 0 ) {
			my $total_num=0;
			my $cnt=1;
			open (my $fh, $ar_orf_output_table) or die "Cannot read $ar_orf_output_table\n";
			while(<$fh>){
				next if (/^ORF_ID/);
				my $ar_info;
				my @temp = split /\t/, $_;
				if($cnt++<=$display_limit){
					#$temp[0] = '<a href=\'http://ardb.cbcb.umd.edu/cgi/search.cgi?db=R&term='."$temp[0]"."'>$temp[0]</a>" if ($temp[-2] eq "ARDB");
					$temp[9] = '<a href=\'https://card.mcmaster.ca/aro/'."$temp[9]"."'>$temp[8]</a>";
					for my $i (0 .. $#temp) {
						my $idx = $i+1;
						$ar_info->{"SGOAR$idx"}=$temp[$i];
					}
					#my @ar_info = sort_by_best_hit($ar_info);
					push @{$vars->{LOOP_SGOAR}}, $ar_info;
				}
				$total_num++;
			}
			close $fh;
			$vars->{SGOARNOTE} = "Only top $display_limit results in terms of Count are listed in the table." if ($total_num>$display_limit);
		}
		if ( -e $vf_orf_output){
			#my $vf_orf_hit_num=`awk 'BEGIN{count=0}\$2>0{count++}END{print count}'  $vf_orf_output`;
			my $vf_orf_hit_num = `awk '{print \$2}' $vf_orf_output | sort | uniq | wc -l`;
			chomp $vf_orf_hit_num;
			my $vf_orf_perc = ($vf_orf_hit_num/$vars->{ANOCDS})*100;
			$vars->{SGOVFHIT} = ($vf_orf_hit_num)? $vf_orf_hit_num : 0;
			$vars->{SGOVFPERC} = $vf_orf_perc;
		}
		if (-e $vf_orf_output_json && $vars->{SGOVFHIT} > 0 ){
			my $json_text;{ local $/; open (my $jsonfh, "<", $vf_orf_output_json) or die "Cannot read json file\n"; $json_text = <$jsonfh>;close $jsonfh;}
			my $jsondata = decode_json($json_text);
			my ($genusCount, $classCount, $vfCount) = 0;
			for my $class (sort keys %$jsondata){
				my $classInfo;
				my $vfHash = $jsondata->{$class};
				#print STDOUT $vfHash." VFHash Class ".$class."\n";
				for my $vf (sort keys %$vfHash){
					my $vfInfo;
					my $detailsArray = $vfHash->{$vf};
					#print STDOUT $detailsArray." DetailsArray VF ".$vf."\n";
					for my $detailsHash (sort { $a->{vfgene} cmp $b->{vfgene} } @{$detailsArray}){
						my $detailInfo;
						$detailInfo->{VFGENE} = $detailsHash->{vfgene};
						my $vfnumber = $detailsHash->{vfnumber};
						$detailInfo->{VFNUMBER} = '<a href=\'http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID='."$vfnumber"."'>$vfnumber</a>";
						my $vfaccession = $detailsHash->{accession};
						$detailInfo->{VFACCESSION} = '<a href=\'http://www.ncbi.nlm.nih.gov/protein/'."$vfaccession"."'>$vfaccession</a>";
						$detailInfo->{VFTAXID} = $detailsHash->{taxid};
						$detailInfo->{VFSOURCE} = $detailsHash->{source};
						$detailInfo->{VFFOUNDIN} = $detailsHash->{foundin};
						$detailInfo->{VFGI} = $detailsHash->{gi};
						$detailInfo->{VFVFDBINFO} = $detailsHash->{vfdbinfo};
						$genusCount++;
						$classCount++;
						$vfCount++;
						#print $genus." ".$class." ".$vf." ".$detailInfo->{VFGENE}." ".$vfCount."\n";
						push @{$vfInfo->{LOOP_DETAIL_RES}}, $detailInfo;
						#print STDOUT $detailInfo." DetailInfo\n";
						#print STDOUT $detailsHash." DetailsHash Details ".$details."\n";
						#for my $vfrow (keys $jsondata->{$genus}->{$class}->{$vf}->[$vfrowindex]){
							#print STDOUT $vfrow."\n";
						#}
					} # End Details
					$vfInfo->{VFVF_COUNT} = $vfCount;
					$vfInfo->{VFVF} = $vf;
					$vfCount = 0;
					push @{$classInfo->{LOOP_VF_RES}},$vfInfo;
				} # End VF
				$classInfo->{VFCLASS_COUNT} = $classCount;
				$classInfo->{VFCLASS} = $class;
				$classCount = 0;
				push @{$vars->{LOOP_SGOVF_CLASS_RES}}, $classInfo;					
			} # End Class
		}
		if ( -e $vf_orf_output_table) {
			my $total_num=0;
			my $cnt=1;
			

			$vars->{SGOVF_KRONA} = "$out_dir/AssemblyBasedAnalysis/SpecialtyGenes/${proj_realname}_VF_genes_ShortBRED.krona.html";



			open (my $fh, $vf_orf_output_table) or die "Cannot read $vf_orf_output_table\n";
			while(<$fh>){
				next if (/^Family/);
				my $vf_info;
				my @temp = split /\t/, $_;
				if($cnt++<=$display_limit){
					$temp[5] = '<a href=\'http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID='."$temp[5]"."'>$temp[5]</a>";
					$temp[8] = '<a href=\'http://www.ncbi.nlm.nih.gov/protein/'."$temp[8]"."'>$temp[8]</a>";
					for my $i (0 .. $#temp) {
						my $idx = $i+1;
						$vf_info->{"SGOVF$idx"}=$temp[$i];
					}
					push @{$vars->{LOOP_SGOVF}}, $vf_info;
				}
				$total_num++;
			}
			close $fh;
			$vars->{SGOVFNOTE} = "Only top $display_limit results in terms of Count are listed in the table." if ($total_num>$display_limit);
		}
	}
}

sub pull_snp_phylogeny {
	return unless -e "$out_dir/SNP_Phylogeny/SNPtree.finished";
	#return unless -e "$out_dir/SNP_Phylogeny/$vars->{SPDB}";
	my $db = $vars->{SPDB};
	my $proj_realname = $vars->{PROJNAME};
	my (@errs,$errs,@warnings); 
	@errs = `grep ERROR $out_dir/SNP_Phylogeny/results/error.log` if ( -e "$out_dir/SNP_Phylogeny/results/error.log");
	@warnings = `grep "Warnings" $out_dir/SNP_Phylogeny/${proj_realname}_PhaME.log` if ( -e "$out_dir/SNP_Phylogeny/${proj_realname}_PhaME.log");
	$errs = join ("<br/>\n",@errs) if (@errs);
	#$vars->{SPCDS} = 1 if ( -e "$out_dir/SNP_Phylogeny/SNPphyloTree.cds.xml");
	$vars->{SPREFLIST} = "$out_dir/SNP_Phylogeny/annotation.txt";
	$vars->{SPTREEALL} = "$out_dir/SNP_Phylogeny/SNPphyloTree.all.xml";
	$vars->{SPTREECDS} = "$out_dir/SNP_Phylogeny/SNPphyloTree.cds.xml" if ( -e  "$out_dir/SNP_Phylogeny/SNPphyloTree.cds.xml" );
	$vars->{SPSNPALL}  = "$out_dir/SNP_Phylogeny/${proj_realname}_stats.txt";
	$vars->{SPALNALL}  = "$out_dir/SNP_Phylogeny/${proj_realname}_all_snp_alignment.fna";
	$vars->{SPALNCDS}  = "$out_dir/SNP_Phylogeny/${proj_realname}_cds_snp_alignment.fna";
	$vars->{SPERROR}   = "$errs" if ($errs);
	$vars->{SPWARN}    = join ("<br/>\n",@warnings) if (@warnings);
	$vars->{SPDIR}     = "$out_dir/SNP_Phylogeny/";
	my $summary_file = "$out_dir/SNP_Phylogeny/${proj_realname}_summaryStatistics.txt";
	if (-e "$summary_file"){
		open (my $snp_summary_fh, $summary_file) or die "Cannot open $summary_file\n";
		while(<$snp_summary_fh>){
			chomp;
			if($_ =~ /Reference used:\s+(.+)/){ $vars->{SPREFNAME} = $1;}
			if($_ =~ /Reference sequence length:\s+(.+)/){ $vars->{SPREFLEN} = $1;}
			if($_ =~ /Total gap length:\s+(.+)/){ $vars->{SPGAPLEN} = $1;}
			if($_ =~ /Core genome length:\s+(.+)/){ $vars->{SPCORELEN} = $1;}
			if($_ =~ /Total SNPs:\s+(.+)/){ $vars->{SPTOTALSNP} = $1;}
			if($_ =~ /CDS SNPs:\s+(.+)/){ $vars->{SPCDSSNP} = $1;}
		}
		close $snp_summary_fh;
	}
}

sub pull_pcr_contig_valid {
	return unless -e "$out_dir/AssayCheck/pcrContigValidation.txt";
	#$vars->{PV_CTG_RESULT} = "FAILURE!";
	
	open(my $pcfh, "<", "$out_dir/AssayCheck/pcrContigValidation.txt") or die $!;
	my $pcr;
	while(<$pcfh>){
		#$vars->{PV_CTG_RESULT} = "SUCCESS!" if /success/;

		chomp;
		if( /Primer Pair (.*)/){
			my $primer_name = $1;
			$primer_name =~ s/ and/,/;
			$pcr->{PVPRIMER}   = $primer_name;
		}
		if( /The primers amplify (\S+) from (\d+) to (\d+), with size (\d+)/ ){
			my ($ctg,$start,$end,$size) = ($1,$2,$3,$4); 
			$pcr->{PVCONTIG}   = $ctg;
			$pcr->{PVLOCATION} = "$start..$end";
			$pcr->{PVPSIZE}    = $size;
		}
		if( /------/){
			$pcr->{PVCONTIG} = $vars->{PROJNAME}."_Assembly" if (!$pcr->{PVCONTIG});
                        $pcr->{PVLOCATION} = "NA" if (!$pcr->{PVLOCATION});
                        $pcr->{PVPSIZE}    = "NA" if (!$pcr->{PVPSIZE});
			push @{$vars->{LOOP_PV}}, $pcr;
			$pcr={};
		}
	}
	close $pcfh;
}

sub pull_pcr_ref_valid {
	return unless -e "$out_dir/AssayCheck/pcrRefValidation.txt";
	#$vars->{PV_REF_RESULT} = "FAILURE!";
	
	open(my $pcfh, "<", "$out_dir/AssayCheck/pcrRefValidation.txt") or die $!;
	my $pcr;
	$pcr->{PVCONTIG} = "Reference";
	while(<$pcfh>){
		#$vars->{PV_REF_RESULT} = "SUCCESS!" if /success/;

		chomp;
		if( /Primer Pair (.*)/){
			my $primer_name = $1;
			$primer_name =~ s/ and/,/;
			$pcr->{PVPRIMER}   = $primer_name;
		}
		if( /The primers amplify (\S+) from (\d+) to (\d+), with size (\d+)/ ){
			my ($ctg,$start,$end,$size) = ($1,$2,$3,$4); 
			$pcr->{PVCONTIG}   = $ctg;
			$pcr->{PVLOCATION} = "$start..$end";
			$pcr->{PVPSIZE}    = $size;
		}else{
		#	$pcr->{PVCONTIG}   = "";
		##	$pcr->{PVLOCATION} = "";
		#	$pcr->{PVPSIZE}    = "";
		}
		if( /------/){
			$pcr->{PVCONTIG} = "Reference" if (!$pcr->{PVCONTIG});
			$pcr->{PVLOCATION} = "N/A" if (!$pcr->{PVLOCATION});
			$pcr->{PVPSIZE}    = "N/A" if (!$pcr->{PVPSIZE});
			push @{$vars->{LOOP_PV}}, $pcr;
			$pcr={};
		}
	}
	close $pcfh;
}


sub pull_pcr_design {
	return unless -e "$out_dir/AssayCheck/PCR.design.primers.txt";
	my %ctg;
	my $cnt=1;

	open(my $pcfh, "<", "$out_dir/AssayCheck/PCR.design.primers.txt") or die $!;
	while(<$pcfh>){
		chomp;
		next if /-- no qualified primers --/;
		next if /^$/;
		next if /Primer Name/;
		
		my @temp = split /\t/, $_;
		if( scalar @temp > 7 ){
			my $pr;
			my ($contig) = $temp[0] =~ /^(.*)-\d+$/;
			next if $ctg{$contig};
			$ctg{$contig}=1;
			$pr->{PDCTG}   = $contig;
			$pr->{PDNAME}  = $temp[0];
			$pr->{PDFWD}   = $temp[1]; 
			$pr->{PDFWDTM} = $temp[2];
			$pr->{PDREV}   = $temp[3];
			$pr->{PDREVTM} = $temp[4];
			$pr->{PDSIZE}  = $temp[5];
			$pr->{PDBG}    = $temp[7];
			$pr->{PDLOC}   = $temp[8];
			push @{$vars->{LOOP_PA_D}}, $pr;
			last if $cnt++ >= 5;
		}
	}
	close $pcfh;
}

sub pull_blast {
	$vars->{ABRESULT} = "All Contigs Blast to NCBI NT Database";
	$vars->{ABRESULT} = "UnmappedContigs Blast to NCBI NR Database" if ($vars->{OUT_RA_SW});
}

sub pull_sra_download {
	my $err;
	if (! -e "$out_dir/SRA_Download/DownloadSRA.finished"){
		$err = `grep -i Failed $out_dir/SRA_Download/log.txt | sed -e 's/\$/<br>/'` if (-e "$out_dir/SRA_Download/log.txt");
		$err .= `grep "ERROR" $out_dir/SRA_Download/log.txt | sed -e 's/\$/<br>/'` if (-e "$out_dir/SRA_Download/log.txt");
	}
	if ($err){
		$err =~ s/(http\S+)/<a href=\"$1\"  target=\"_blank\">$1<\/a>/;
		$vars->{SRADERR} = $err;
	}
}

sub pull_contig_taxa {
	return unless -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/summary_by_topHitCount.txt";
	open (my $ccpfh, "<", "$out_dir/AssemblyBasedAnalysis/Taxonomy/summary_by_topHitCount.txt") or die $!;
	my $header=<$ccpfh>;
	while(<$ccpfh>){
		chomp;
		my @temp = split /\t/, $_;
		my $row;
		$row->{CCPRANK}=$temp[0];
		if ($temp[0] eq "species"){ 
			for my $i (1..5){
				next if ($temp[$i] =~ /unclassified/ or $temp[$i] eq "NA" );
				$temp[$i]="<a href=\'http://www.ncbi.nlm.nih.gov/genome/?term=\"$temp[$i]\"\' target='_blank' >$temp[$i]</a>"
			}
		}
		$row->{CCPTOP1}=$temp[1];
		$row->{CCPTOP2}=$temp[2];
		$row->{CCPTOP3}=$temp[3];
		$row->{CCPTOP4}=$temp[4];
		$row->{CCPTOP5}=$temp[5];
		push @{$vars->{LOOP_CCP}}, $row;
	}
	close $ccpfh;
	return unless -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/summary_by_hitAccLength.txt";
	open (my $ccplfh, "<", "$out_dir/AssemblyBasedAnalysis/Taxonomy/summary_by_hitAccLength.txt") or die $!;
	$header=<$ccplfh>;
	while(<$ccplfh>){
		chomp;
		my @temp = split /\t/, $_;
		my $row;
		$row->{CCPLRANK}=$temp[0];
		$row->{CCPLTOP1}=$temp[1];
		$row->{CCPLTOP2}=$temp[2];
		$row->{CCPLTOP3}=$temp[3];
		$row->{CCPLTOP4}=$temp[4];
		$row->{CCPLTOP5}=$temp[5];
		push @{$vars->{LOOP_CCPL}}, $row;
	}
	close $ccplfh;
	return unless -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/log.txt";
	open (my $ccLog, "<", "$out_dir/AssemblyBasedAnalysis/Taxonomy/log.txt") or die $!;
	while(<$ccLog>){
		chomp;
		if ( $_ =~ /Total Contigs: (\d+) \((\d+) bp\); Classified Contigs: (\d+) \((\d+) bp\); Unclassified Contigs: (\d+) \((\d+) bp\);/){
			$vars->{CCPNTC}  = $1;
			$vars->{CCPNTCB} = $2;
			$vars->{CCPNCC}  = $3; 
			$vars->{CCPNCCB} = $4;
			$vars->{CCPNUC}  = $5;
			$vars->{CCPNUCB} = $6;
		}
	}
	close $ccLog;
	my $proj_realname = $vars->{PROJNAME};
	return unless -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/$proj_realname.ctg_class.top.csv";
	open (my $ccp_top, "<", "$out_dir/AssemblyBasedAnalysis/Taxonomy/$proj_realname.ctg_class.top.csv") or die $!;
	my %uniq;
	while(<$ccp_top>){
		chomp;
		my @array = split /\t/,$_;
		next if ($array[1] ne "genus");
		$uniq{$array[2]}++;
	}
	close $ccp_top;
	foreach my $spe (sort { $uniq{$b} <=> $uniq{$a} } keys %uniq){
		my $row;
		$row->{CCPSPE}=$spe;
		$row->{CCPSPENUM}=$uniq{$spe};
		push @{$vars->{LOOP_CCPSPE}},$row;
	}
	if ( ! -e "$out_dir/AssemblyBasedAnalysis/Taxonomy/$proj_realname.ctg_class.LCA.json" ){
		my $row_limit = $sys->{edgeui_result_table_rows} || 3000;
		system("perl $RealBin/../tab2Json_for_dataTable.pl -project_dir $out_dir -mode contig -limit $row_limit $out_dir/AssemblyBasedAnalysis/Taxonomy/$proj_realname.ctg_class.LCA.csv > $out_dir/AssemblyBasedAnalysis/Taxonomy/$proj_realname.ctg_class.LCA.json");
	}
}
sub pull_taxa {
	my $um_reads="";  
	my $reads_type="";
	$um_reads      = "."             if $vars->{OUT_AL_CP_SW};
	$reads_type    = "allReads"      if $vars->{OUT_AL_CP_SW};
	$um_reads      = "UnmappedReads" if $vars->{OUT_UM_CP_SW};
	$reads_type    = "UnmappedReads" if $vars->{OUT_UM_CP_SW};

	my $num_input_reads = $NUM_READS_FOR_DOWNSTREAM;
	$num_input_reads = $vars->{RMREFUNMAPPED} if $reads_type eq "UnmappedReads";
	$num_input_reads =~ s/ \(.*\)//;

	my $cnt=0;

	my $temp_num_input_reads = $num_input_reads;
	$temp_num_input_reads = &_reformat_val($temp_num_input_reads);
	$vars->{CPRESULT} = "Profiling results of $temp_num_input_reads unmapped reads to the reference.";
	$vars->{CPRESULT} = "Profiling results of all $temp_num_input_reads reads." if (defined $reads_type and $reads_type eq "allReads");

	return unless -e "$out_dir/ReadsBasedAnalysis/$um_reads/Taxonomy/report/summary.txt";

	$vars->{CPDIR}         = "$out_dir/ReadsBasedAnalysis/$um_reads/Taxonomy/report";
	$vars->{CPRADAR}       = "$vars->{CPDIR}/radarchart_DATASET_$reads_type.species.html";
	$vars->{CPHEATMAP}     = "$vars->{CPDIR}/heatmap_DATASET-$reads_type.species.pdf";
	$vars->{CPHEATMAP_PNG} = "$out_dir/HTML_Report/images/heatmap_DATASET-$reads_type.species.png";
	$vars->{CPSUMMARY}     = "$vars->{CPDIR}/report_summary.xlsx";
	$vars->{CPABUN}        = "$vars->{CPDIR}/report_SEQ1_$reads_type.xlsx";

	$vars->{CPRADAR}   = "" unless -e $vars->{CPRADAR};
	$vars->{CPHEATMAP} = "" unless -e $vars->{CPHEATMAP};

	open(my $cpfh, "<", "$out_dir/ReadsBasedAnalysis/$um_reads/Taxonomy/report/summary.txt") or die $!;
	while(<$cpfh>){
		chomp;
		my @temp = split /\t/, $_;
		my $row;
		next if $temp[1] =~ /gottcha\d?-(\w{3})DB/ && $temp[2] !~ /$1/;
		next if $temp[1] =~ /^TOOL$/;

		if( scalar @temp > 5 ){
			$row->{CPTOOL}=$temp[1];
			$row->{CPRANK}=$temp[2];
			$row->{CPTOP1}=$temp[3];
			$row->{CPTOP2}=$temp[4];
			$row->{CPTOP3}=$temp[5];
			$row->{CPTOP4}=$temp[6];
			$row->{CPTOP5}=$temp[7];
			
			### calculate classified reads
			my $toolname = $row->{CPTOOL};
			my $creads = 0; #classified reads
			my $cur_level ="";
			my $abu_list = "$vars->{CPDIR}/1_$reads_type/$row->{CPTOOL}/$reads_type-$row->{CPTOOL}.list.txt"; 
			if ( -e $abu_list){
				open LIST, $abu_list or die $!;
				while(<LIST>){
					next if /^LEVEL/;
					my @temp = split /\t/, $_;
					$cur_level ||= $temp[0];
					last if $cur_level ne $temp[0];
					my $mapped_reads = $toolname =~ /gottcha-/i ? $temp[8] : $temp[2];
					$creads += $mapped_reads if ($mapped_reads =~ /\d/);
					$creads = "N/A" if $toolname =~ /metaphlan/i;
				}
				close LIST;
			}
			$row->{CPCNUM} = $creads;
			$row->{CPCPCT} = sprintf "%.1f", ($creads/$num_input_reads*100) if $creads =~ /\d+/;
			
			push @{$vars->{LOOP_CP}}, $row;

			### tool result
			next if $temp[2] ne "species";
			my $tool;
			$tool->{CPTOOL_CPABU_PMD} = 1 if $toolname =~ /gottcha-/i;
			$tool->{CPTOOL_GOTTCHA2} = 1 if $toolname =~ /gottcha2-/i;
			$tool->{CPTOOL_PANGIA} = 1 if $toolname =~ /pangia/i;
			
			my $count=0;
			if (-e $abu_list){
				open LIST, $abu_list or die $!;
				while(<LIST>){
					chomp;
					my @t = split /\t/, $_;
					my $res_row;
					
					if( $t[0] eq "species"){
						$res_row->{CPABU_LVL} = $t[0];
						$res_row->{CPABU_TAX} = ($t[1] =~ /unclassied|unassign/i or $t[1] eq "NA" )? $t[1]:
									"<a href=\'http://www.ncbi.nlm.nih.gov/genome/?term=\"$t[1]\"\' target='_blank'>$t[1]</a>";
						$res_row->{CPABU_DOWNLOAD_TAX} = $t[1];
						$res_row->{CPABU_DOWNLOAD_TAX_ID} = $t[1];
					
						if( $toolname =~ /gottcha-/ ){
							$res_row->{CPABU_REA} = _reformat_val($t[8]);
							$res_row->{CPABU_PMD} = sprintf "%.1f", ($t[7]/$t[6]*100);
							$res_row->{CPTOOL_CPABU_PMD} = 1;
							$res_row->{CPABU_ABU} = sprintf "%.1f", ($t[2]*100);
						}
						elsif( $toolname =~ /gottcha2/ ){
							$res_row->{CPABU_REA} = _reformat_val($t[2]);
							$res_row->{GOTTCHA2_LINEAR_LEN} = _reformat_val($t[8]);
							$res_row->{GOTTCHA2_LINEAR_DOC} = sprintf "%.2f", $t[9];
							$res_row->{CPABU_ABU} = sprintf "%.1f", ($t[10]*100);
							$res_row->{CPABU_DOWNLOAD_TAX_ID} = $t[4];
						}
						elsif( $toolname =~ /pangia/ ){
							$res_row->{CPABU_REA} = _reformat_val($t[2]);
							$res_row->{PANGIA_LINEAR_COV}  = sprintf "%.2f", $t[16];
							$res_row->{PANGIA_LINEAR_DOC}  = sprintf "%.2f", $t[11];
							$res_row->{PANGIA_LINEAR_RSNB} = _reformat_val(sprintf "%.2f", $t[7]);
							$res_row->{PANGIA_LINEAR_SCR}  = $t[13];
							$res_row->{CPABU_ABU} = sprintf "%.1f", ($t[14]*100);
							$res_row->{CPABU_DOWNLOAD_TAX_ID} = $t[4];
						}
						elsif( $toolname =~ /metaphlan/ ){
							$res_row->{CPABU_REA} = "N/A";
							$res_row->{CPABU_ABU} = $t[2];
						}
						else{
							$res_row->{CPABU_REA} = _reformat_val($t[2]);
							$res_row->{CPABU_ABU} = sprintf "%.1f", ($t[2]/$creads*100);
						}
						if ($count <5){
							push @{$tool->{LOOP_CPTOOL_RES}}, $res_row;
						}
						push @{$tool->{LOOP_CPTOOL_DOWNLOAD}},$res_row;
						++$count;
					}
				}
				close LIST;
			}
			$tool->{CPTOOL_ID}    = $cnt++;
			$tool->{CPTOOL_LABEL} = $row->{CPTOOL};
			$tool->{CPTOOL_LABEL} =~ s/\b(\w)/\U$1/g;
			$tool->{CPTOOL_LABEL} = "GOTTCHA (bacterial species database)" if $row->{CPTOOL} =~ /gottcha-.*-b/;
			$tool->{CPTOOL_LABEL} = "GOTTCHA (viral species database)"     if $row->{CPTOOL} =~ /gottcha-.*-v/;
			$tool->{CPTOOL_LABEL} = "GOTTCHA2 (bacterial species database)" if $row->{CPTOOL} =~ /gottcha2-.*-b/;
			$tool->{CPTOOL_LABEL} = "GOTTCHA2 (viral species database)"     if $row->{CPTOOL} =~ /gottcha2-.*-v/;
			$tool->{CPTOOL_LABEL} = "Kraken (mini database)"       if $row->{CPTOOL} =~ /kraken_mini/;
			$tool->{CPTOOL_LABEL} = "BWA (reads mapping)"          if $row->{CPTOOL} =~ /bwa/;
			$tool->{CPTOOL_LABEL} = "PanGIA"          if $row->{CPTOOL} =~ /pangia/;

			$tool->{CPTOOL}       = $row->{CPTOOL};
			$tool->{CPTOOL_TREE}  = "$vars->{CPDIR}/1_$reads_type/$row->{CPTOOL}/$reads_type-$row->{CPTOOL}.tree.svg";
			$tool->{CPTOOL_TREE_PNG} = "$out_dir/HTML_Report/images/$reads_type-$row->{CPTOOL}.tree.png";
			$tool->{CPTOOL_KRONA} = "$vars->{CPDIR}/1_$reads_type/$row->{CPTOOL}/$reads_type-$row->{CPTOOL}.krona.html";
			$tool->{CPTOOL_ABU}   = "$vars->{CPDIR}/1_$reads_type/$row->{CPTOOL}/$reads_type-$row->{CPTOOL}.list.txt";
			$tool->{CPTOOL_DIR}   = "$vars->{CPDIR}/1_$reads_type/$row->{CPTOOL}";

			$tool->{CPTOOL_TREE}  = "" unless -e $tool->{CPTOOL_TREE};
			$tool->{CPTOOL_KRONA} = "" unless -e $tool->{CPTOOL_KRONA};
			$tool->{CPABU_DOWNLOAD_LIST} = 1 if ($row->{CPTOOL}=~/gottcha|bwa|pangia/);

			push @{$vars->{LOOP_CPTOOL}}, $tool;
		}
	}
	close $cpfh;
}

sub pull_readmapping_contig {
	return unless -e "$out_dir/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.alnstats.txt";	
	## Reads to Contigs
	open(my $rmfh, "<", "$out_dir/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.alnstats.txt") or die $!;
	while(<$rmfh>) {
		if ($_ =~ /^(\d+) \+ \d+ in total/) 	{ $vars->{RMUSED} = $1; next; }
		if ($_ =~ /^(\d+) \+ \d+ duplicates/) 	{ $vars->{RMDUPS} = $1; next; }
		if ($_ =~ /^(\d+) \+ \d+ mapped/) 		{ $vars->{RMMAPPED} = $1; $vars->{RMUNMAPPED} = $vars->{RMUSED} - $vars->{RMMAPPED}; next; }
		if ($_ =~ /^Avg_coverage_fold:\t(\d+\.\d+)/)	{ $vars->{RMAVECOV} = $1; next; }
		if ($_ =~ /^Coverage:\t(\d+\.\d+\%)/)	{ $vars->{RMTOTALCOV} = $1; next; }
	}
	close($rmfh);
	
	$vars->{RMMAPPEDPCT}   = sprintf "%.2f", $vars->{RMMAPPED}/$vars->{RMUSED}*100;
	$vars->{RMUNMAPPED}    = $vars->{RMUSED} - $vars->{RMMAPPED};
	$vars->{RMUNMAPPEDPCT} = sprintf "%.2f", $vars->{RMUNMAPPED}/$vars->{RMUSED}*100;
	if ( ! -e "$out_dir/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs_coverage.table.json"){
		my $row_limit = $sys->{edgeui_result_table_rows} || 3000;
		system("perl $RealBin/../tab2Json_for_dataTable.pl -project_dir $out_dir -mode contig -limit $row_limit $out_dir/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs_coverage.table > $out_dir/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs_coverage.table.json");
	}
}

sub pull_readmapping_ref {
	return unless -e "$out_dir/ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt";
	
	my ($tol_ref_number, $tol_ref_hashit, $tol_ref_len, $tol_mapped_bases, $tol_non_gap_bases) = (0,0,0,0,0);
	my ($tol_snps, $tol_indels)= (0,0);
	my $ref_display_limit = 30;
	my $ref_display_limit_plot = 4;
	my $ref;

	open(my $reffh, "<", "$out_dir/ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt") or die $!;
	while(<$reffh>) {
		if ($_ =~ /^(\d+) \+ \d+ in total/) 	{ $vars->{RMREFUSED} = $1; next; }
		if ($_ =~ /^(\d+) \+ \d+ duplicates/) 	{ $vars->{RMREFDUPS} = $1; next; }
		if ($_ =~ /^(\d+) \+ \d+ mapped/) 		{ $vars->{RMREFMAPPED} = $1; $vars->{RMREFUNMAPPED} = $vars->{RMREFUSED} - $vars->{RMREFMAPPED}; next; }
	
		next if( /^Ref/ );
		my @temp = split /\t/, $_;
		if( scalar @temp == 11 ){
			my $refinfo;
			$tol_ref_number++;
			$tol_ref_hashit++ if $temp[3];
			$tol_ref_len += $temp[1];
			$tol_mapped_bases += $temp[1] * $temp[5];
			$tol_non_gap_bases += $temp[1] - $temp[8];
			$tol_snps += $temp[9];
			$tol_indels += $temp[10];
			next if $temp[3] == 0;
			for my $i (0 .. $#temp) {
				my $idx = $i+1;
				
				$refinfo->{"RMREFT$idx"}=$temp[$i];
			}
			$refinfo->{"RMREFNAME"}=$refname->{$temp[0]};
			$ref->{$temp[0]}=$refinfo;
		}
	}
	close($reffh);
	$vars->{RMAVEFOLD} = sprintf "%.2f", $tol_mapped_bases/$tol_ref_len;
	$vars->{RMREFCOV} = sprintf "%.2f", $tol_non_gap_bases/$tol_ref_len*100;
	$vars->{RMSNPS} = $tol_snps;
	$vars->{RMINDELS} = $tol_indels;	

	my $cnt=0;
	foreach my $n ( sort {$ref->{$b}->{RMREFT5} <=> $ref->{$a}->{RMREFT5}} keys %$ref )
	{
		push @{$vars->{LOOP_RMREF}}, $ref->{$n};

		my $plot_info;
		my $pid = $ref->{$n}->{RMREFT1};
		$plot_info->{RMREF_PLOT_ID}       = $pid;
		$plot_info->{RMREF_PLOT_COVPATH}  = "$out_dir/ReadsBasedAnalysis/readsMappingToRef/Coverage_plots/readsToRef_${pid}_base_coverage.png";
		$plot_info->{RMREF_PLOT_FOLDPATH} = "$out_dir/ReadsBasedAnalysis/readsMappingToRef/Coverage_plots/readsToRef_${pid}_coverage_histogram.png";
		push @{$vars->{LOOP_RMREF_PLOT}}, $plot_info if $cnt < $ref_display_limit_plot;
		
		last if ++$cnt == 30;
	}

	$vars->{RMREFMAPPEDPCT}   = sprintf "%.2f", $vars->{RMREFMAPPED}/$vars->{RMREFUSED}*100;
	$vars->{RMREFUNMAPPED}    = $vars->{RMREFUSED} - $vars->{RMREFMAPPED};
	$vars->{RMREFUNMAPPEDPCT} = sprintf "%.2f", $vars->{RMREFUNMAPPED}/$vars->{RMREFUSED}*100;
	$vars->{RMREFTOLREF}      = $tol_ref_number;
	$vars->{RMREFTOLREFHASHIT} = $tol_ref_hashit;
	$vars->{RMREFTABLENOTE} = "Only top $ref_display_limit results in terms of \"Base Recovery %\" are listed in the table." if $tol_ref_number > $ref_display_limit;
	$vars->{RMREFSNPFILE}     = 0;
	$vars->{RMREFSNPFILE}     = 1 if -e "$out_dir/ReadsBasedAnalysis/readsMappingToRef/readsToRef.SNPs_report.txt";
	$vars->{RMREFGAPFILE}     = 1 if -e "$out_dir/ReadsBasedAnalysis/readsMappingToRef/GapVSReference.report.json";

	## display unmapped reads mapping to RefSeq
	my $tol_um_org=0;

	$vars->{RMREF_UM} = 1 if -e "$out_dir/ReferenceBasedAnalysis/UnmappedReads/UnmappedReads_coverage.txt";

	if( $vars->{RMREF_UM} )
	{
		open(my $reffh, "<", "$out_dir/ReferenceBasedAnalysis/UnmappedReads/UnmappedReads_coverage.txt") or die $!;
		while(<$reffh>) {
			next if( /^GI/ );
			last if $tol_um_org == 5;
			my @temp = split /\t/, $_;
			my $refinfo;
			($refinfo->{"RMREFNAMESHORT"}) = $temp[8] =~ /^(\S+ \S+)/;
			$refinfo->{"RMREF_UM_FULL"} = $temp[8];
			for my $i (1..7) {
				my $idx = $i+1;
				$refinfo->{"RMREFT$idx"}=$temp[$i];
			}
			$refinfo->{"RMREFNAME"}=$refname->{$temp[0]};
			$ref->{$temp[0]}=$refinfo;
			push @{$vars->{LOOP_RMREF_UMR}}, $refinfo;
			$tol_um_org++;
		}
		close($reffh);

		$vars->{RMREF_UM_NOTE} = "Only top 5 results in terms of \"Mapped Reads\" are listed in the table.";
	}
}

sub pull_summary {
	my @INFILES;
	my $cnt=0;
	my $cpu=0;
	my $tol_running_sec=0;
	my $prog;
	my %map;
	my ($step,$ord,$do,$status,$lastline,$kingdom);

	$vars->{PROJSTATUS} = "Unfinished";

	open(my $sumfh, "cat $out_dir/config.txt $out_dir/process.log|") or die $!;
	while(<$sumfh>) {
		chomp;
		#parse input files
		if( /runPipeline/ ) {
			undef @INFILES;	
		}
		$vars->{LOG_qiime_SW} = 1 if (/qiime_pipeline/);
		if( /runPipeline .*-p (.*) -\w/ || /runPipeline .*-p (.*) >/ || /runPipeline .*-p (.*)$/) {
			push @INFILES, split /\s+/,$1;
		}
		if(/runPipeline .*-u (.*) -\w/ || /runPipeline .*-u (.*) >/ || /runPipeline .*-u (.*)$/){
			push @INFILES, split /\s+/,$1;
		}
		if(/SRA_id=(.*)/)
		{
			push @INFILES, $1;
		}
		#parse reference files
		if( /runPipeline .*-ref (\S+)/){
			$vars->{REFFILE} = $1;
		}
		elsif( /^reference=(\S+)/ ){
			 $vars->{REFFILE} = $1;
		}

		if( /Total Running time: (\d+:\d+:\d+)/){
			$vars->{LASTRUNTIME} = $1;
			next;
		}
		elsif( /^Host=(.*)/ ){
			$vars->{HOSTFILE} = $1;
			next;
		}
		elsif( /^SNPdbName=(.*)/ ){
			$vars->{SPDB} = $1;
			next;
		}
		if (/^kingdom=(.*)/){
			$kingdom = $1;
		}
		if (/^assembler=(.*)/){
			$vars->{ASSEMBLER}=$1;
		}
		if (/^assembledContigs=(.*)/){
			$vars->{ASSEMBLEDCONTIG}=$1;
		}
		if( /^\[(.*)\]/ ){
			my $step = $1;
			if( $step eq "project" or $step eq "system"){
				while(<$sumfh>){
					chomp;
					if ( /^([^=]+)=([^=]+)/ ){
						$vars->{uc($1)}=$2;
					}
					elsif ( /^\[(.*)\]/ ){
						$step = $1;
						last;
					}
				}
			}
			
			if( defined $map{"$step"} ){
				$ord = $map{"$step"};
			}
			else{
				$cnt+=100;
				$prog->{$cnt}->{GNLANALYSIS}=$step;
				$map{"$step"}=$cnt;
			}
		}
		elsif( /Project Start: (.*)/ ){
			$vars->{PROJSUBTIME} = $1;
			$vars->{PROJSTATUS} = "Finished";
		}
		elsif( /^Do.*=(.*)$/ ){
			my $do = $1;
			$do =~ s/ //g;
			$prog->{$cnt}->{GNLRUN}= "Auto" if $do eq 'auto';
			$prog->{$cnt}->{GNLRUN}= "On" if $do eq 1;
			$prog->{$cnt}->{GNLRUN}= "Off" if $do eq 0 && !$prog->{$cnt}->{GNLRUN};
			
			$prog->{$cnt}->{GNLSTATUS}="Skipped";
			$prog->{$cnt}->{GNLSTATUS}="Incomplete" if $prog->{$cnt}->{GNLRUN} eq 'On';
			$prog->{$cnt}->{GNLSTATUS}= "Skipped" if ($prog->{$cnt}->{GNLANALYSIS} =~ /ProPhage/i && ($kingdom =~ /virus/i || /No CDS annotation/));
		}
		elsif( /Produce Final PDF Report/){
			my $pdf_generate_time=<$sumfh>;
			my ($h,$m,$s) = $pdf_generate_time =~ /(\d+):(\d+):(\d+)/;
			$tol_running_sec += $h*3600+$m*60+$s if ( defined $h);
		}
		elsif( /Finished/ ){
			$prog->{$ord}->{GNLSTATUS} = "Skipped (result exists)";
		}
		elsif( /Running time: (\d+:\d+:\d+)/ ){
			$prog->{$ord}->{GNLSTATUS} = "Complete";
			my ($h,$m,$s) = $1 =~ /(\d+):(\d+):(\d+)/;
			$tol_running_sec += $h*3600+$m*60+$s;
			$prog->{$ord}->{GNLTIME} += $h*3600+$m*60+$s;
		}
		elsif( / Running/ ){
			$prog->{$ord}->{GNLSTATUS} = "<span class='edge-fg-orange'>Running</span>";
			$vars->{PROJSTATUS} = "<span class='edge-fg-orange'>Running</span>";
		}
		elsif( /failed/ ){
			$prog->{$ord}->{GNLSTATUS} = "<span class='edge-fg-red'>Failed</span>";
			$vars->{PROJSTATUS} = "<span class='edge-fg-red'>Failure</span>";
			$vars->{PROJSTATUS} = "<span class='edge-fg-red'>Failure. The process has failed unexpectedly. Please contact system administrator. Or try to rerun the job</span>" if (/Unexpected exit/);
		}
		elsif( /All Done/){
			$prog->{$ord}->{GNLSTATUS} = "Complete";
			$vars->{PROJSTATUS} = "Complete";
		}
		$lastline = $_;
	}

	#convert running time to h:m:s
	foreach my $o ( keys %$prog ){
		if (defined $prog->{$o}->{GNLTIME}){
			my $step_sec = $prog->{$o}->{GNLTIME};
			$prog->{$o}->{GNLTIME} = sprintf("%02d:%02d:%02d", int($step_sec / 3600), int(($step_sec % 3600) / 60), int($step_sec % 60));
		}
	}
	#$vars->{RUNTIME} = strftime("\%H:\%M:\%S", gmtime($tol_running_sec));
	$vars->{RUNTIME} = sprintf("%02d:%02d:%02d", int($tol_running_sec / 3600), int(($tol_running_sec % 3600) / 60), int($tol_running_sec % 60));

	$vars->{PROJSTATUS}        = "Unstarted"   if $lastline =~ /EDGE_UI.*unstarted/;
	$vars->{PROJSTATUS}        = "Interrupted" if $lastline =~ /EDGE_UI.*interrupted/;
	$vars->{PROJSTATUS}        = "Archived"	   if $lastline =~ /EDGE_UI.*archived/;
	$prog->{$ord}->{GNLSTATUS} = "Interrupted" if $vars->{PROJSTATUS} eq "Interrupted"; #turn last step to unfinished
	$vars->{PROJSTATUS}        = "Complete"    if $mode ne "web"; # run by end of runPipeline  
	$prog->{$ord}->{GNLSTATUS} = "Complete" if $mode ne "web";
	
	#Reads Taxonomy Classification
	$ord = $map{"Reads Taxonomy Classification"};
	$ord = $map{"Qiime analysis"} if ($vars->{LOG_qiime_SW});
	my %toolmap;
	my $step_runtime;
	open PROC_CUR, "<", "$out_dir/process_current.log" or die $!;
	while(<PROC_CUR>) {
		chomp;
		#parse input files
		if( /^\[RUN_TOOL\] \[(.*)\] (COMMAND|Logfile)/ ||  /^Qiime \[.*\]\s+(.*)/ || /^\[RUN_TOOL\] \[(.*)\] (Result exists)/){
			$step = $1;
			next if defined $toolmap{$step};
			$ord++;
			$toolmap{$step} = $ord;
			$prog->{$ord}->{GNLANALYSIS} = "<span style='margin-left:3em'>$step</span>";
			$prog->{$ord}->{GNLRUN}      = "On";
			$prog->{$ord}->{GNLSTATUS}   = "<span class='edge-fg-orange'>Running</span>";
			$prog->{$ord}->{GNLSTATUS}   = "Skipped (result exists)" if ($2 and $2 =~ /Result exists/);
		}
		elsif( /^\[RUN_TOOL\] \[(.*)\] Error occured/){
			my $ord = $toolmap{$1};
			$prog->{$ord}->{GNLSTATUS}   = "<span class='edge-fg-red'>Error</span>";
		}
		elsif( /^\[RUN_TOOL\] \[(.*)\] Running time: (.*)/){
			my $ord = $toolmap{$1};
			$prog->{$ord}->{GNLTIME}     = $2;
			$prog->{$ord}->{GNLSTATUS}   = "Complete" if ( $prog->{$ord}->{GNLSTATUS} !~ /Error|Skip/i);
		}
		elsif(/ERROR|failed/ and $vars->{LOG_qiime_SW}){
			$prog->{$ord}->{GNLSTATUS}   = "<span class='edge-fg-red'>Error</span>";
			$vars->{ERROR_qiime} = $_ if (/ERROR/i and $step ne "Checking Mapping File");
		}
		elsif( /^Qiime Running time: (.*)/){
			$prog->{$ord}->{GNLSTATUS}   = "Complete";
			my ($h,$m,$s) = $1 =~ /(\d+):(\d+):(\d+)/;
			if ($prog->{$ord}->{GNLTIME}){
				$step_runtime += $h*3600+$m*60+$s;
			}else{
				$step_runtime = $h*3600+$m*60+$s;
			}
			$prog->{$ord}->{GNLTIME}     = sprintf("%02d:%02d:%02d", int($step_runtime / 3600), int(($step_runtime % 3600) / 60), int($step_runtime % 60));
		}
	}
	close PROC_CUR;
	$prog->{$ord}->{GNLSTATUS} = "Interrupted" if $vars->{PROJSTATUS} eq "Interrupted"; #turn last step to unfinished

	$getting_paired = 1 if scalar @INFILES > 1;
	map {  $_=basename($_) if ($_); } @INFILES;
	$vars->{INFILES} = join ", ", @INFILES;
	my @HOST_FILES;
	map {  if ($_){ my $host=basename($_); push @HOST_FILES,$host;}} split /\,/, $vars->{HOSTFILE};
	$vars->{HOSTFILE} = join ", ", @HOST_FILES;
	
	foreach my $ord ( sort {$a<=>$b} keys %$prog ){
		next if $prog->{$ord}->{GNLSTATUS} eq "Skipped";
		push @{$vars->{LOOP_GNL}}, $prog->{$ord};
	}
	
	close ($sumfh);
}
