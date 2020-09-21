#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../../lib";
use Spreadsheet::ParseExcel;
use Spreadsheet::ParseExcel::SaveParser;
use Getopt::Long;
use File::Basename;
use File::Path;
use File::Copy;

my $um; #user management
my $out;
my $project_dir_names;
my $userDir;
my $usage = qq{
Usage: $0
	Required
		-out        	      output file
		-projects        project_dir_names, separated by comma
		-udir		 user dir
};

GetOptions(
		"um=s"        =>  \$um,
		"out=s"        =>  \$out,
		"projects=s"       =>  \$project_dir_names,
		"udir=s"       => \$userDir,
		"help|?"           => sub{print "$usage\n";exit;} 
	);

if (!$project_dir_names && !$out){ print "$usage\n";exit;}

my ($file_prefix, $outputDir, $file_suffix)=fileparse("$out", qr/\.[^.]*/);
mkpath("$outputDir/../NCBI");
my $seqout = "$outputDir/all_sequences.fasta";
my $seqoutNCBI = "$outputDir/../NCBI/all_sequences_ncbi.fasta";
# clean up before start
unlink $seqoutNCBI;
my $source_tsvout = "$outputDir/../NCBI/source.src";
my $comment_tsvout = "$outputDir/../NCBI/comment.cmt";
## https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
my $submission_template = "$outputDir/../NCBI/template.sbt";
## not ready
my $submitter_profile =  "$outputDir/../NCBI/submitter_profile.txt";
copy("$userDir/gisaid_ncbi_submission_profile.txt",$submitter_profile);
my $submission_xml = "$outputDir/../NCBI/submission.xml";
unlink $seqout;

## read template
my $parser   = Spreadsheet::ParseExcel::SaveParser->new();
my $template = $parser->Parse("$FindBin::Bin/20200717_EpiCoV_BulkUpload_Template.xls");

if ( !defined $template ) {
    die $parser->error(), ".\n";
}
my $in_worksheet1 = $template->worksheet('Submissions');
my $in_worksheet0 = $template->worksheet('Instructions');
my $row=2;
my $col=0;
my $seqID="seq0";

my @src_tsv_header = ("Sequence_ID","Organism","isolate","collection-date","country","host","Collected-By");
my @src_tsv_content;
my @cmt_tsv_header = ("SeqID","StructuredCommentPrefix","Assembly Method","Coverage","Sequencing Technology","StructuredCommentSuffix");
my @cmt_tsv_content;

#write metadata to sheets
foreach my $proj_dir (split /,/,$project_dir_names){
	my $vars={};
	$seqID++;
	my $confFile = "$proj_dir/config.txt";
	my $conf = &getParams($confFile);
	#my $otherFile = "$proj_dir/metadata_other.txt";

	my $proj_name = $conf->{'projname'};
	my $owner = $conf->{'projowner'};
	eval {  
		&pull_sampleMetadata($proj_dir,$vars);
		&pull_submissionData($proj_dir,$vars);
		&pull_consensusInfo($proj_dir,$vars,$conf);
		&write_all_sequences($seqout,$vars);
		&write_all_sequences($seqoutNCBI,$vars,$seqID);
	};
	$in_worksheet1->AddCell( $row, $col, $vars->{ID});  ## Submitter (login account id) *
	$in_worksheet1->AddCell( $row, $col+1, "all_sequences.fasta");  ## FASTA filename *
	$in_worksheet1->AddCell( $row, $col+2, $vars->{VIR_NAME});  ## Virus Name * hCoV-19/Country/Identifier/2020, ex:hCoV-19/USA/NM-UNM-00001/2020"  
	$in_worksheet1->AddCell( $row, $col+3, "betacoronavirus");  ## Type *
	$in_worksheet1->AddCell( $row, $col+4, $vars->{VIR_PASSAGE}); ## Passage details/history *
	$in_worksheet1->AddCell( $row, $col+5, $vars->{SM_CDATE}); ## Collection date *
	$in_worksheet1->AddCell( $row, $col+6, $vars->{SM_LOC}); ## Location *
	$in_worksheet1->AddCell( $row, $col+7, ""); ## Additional location information
	$in_worksheet1->AddCell( $row, $col+8, $vars->{SM_HOST});  ## Host *
	$in_worksheet1->AddCell( $row, $col+9, "");  ## Additional host information
	$in_worksheet1->AddCell( $row, $col+10, $vars->{SM_GENDER});  ## Gender *
	$in_worksheet1->AddCell( $row, $col+11, $vars->{SM_AGE});  ## Patient Age *
	$in_worksheet1->AddCell( $row, $col+12, $vars->{SM_STATUS});  ## Patient status *
	$in_worksheet1->AddCell( $row, $col+13, "");  ## Specimen source, Sputum, Alveolar lavage fluid, Oro-pharyngeal swab, Blood, Tracheal swab, Urine, Stool, Cloakal swab, Organ, Feces, Other
	$in_worksheet1->AddCell( $row, $col+14, "");  ## Outbreak
	$in_worksheet1->AddCell( $row, $col+15, "");  ## Last vaccinated
	$in_worksheet1->AddCell( $row, $col+16, "");  ## Treatment
	$in_worksheet1->AddCell( $row, $col+17, $vars->{SM_SEQUENCING_TECH});  ## Sequencing technology *
	$in_worksheet1->AddCell( $row, $col+18, $vars->{ASM_METHOD});  ## Assembly method 
	$in_worksheet1->AddCell( $row, $col+19, $vars->{SM_COV});  ## Coverage
	$in_worksheet1->AddCell( $row, $col+20, $vars->{ORIG_LAB});  ## Originating lab *
	$in_worksheet1->AddCell( $row, $col+21, $vars->{ORIG_ADDRESS});  ## Address *
	$in_worksheet1->AddCell( $row, $col+22, "");  ## Sample ID given by the sample provider
	$in_worksheet1->AddCell( $row, $col+23, $vars->{SUB_LAB});  ## Submitting lab *
	$in_worksheet1->AddCell( $row, $col+24, $vars->{SUB_ADDRESS});  ## Address *
	$in_worksheet1->AddCell( $row, $col+25, "");  ## Sample ID given by the submitting laboratory
	$in_worksheet1->AddCell( $row, $col+26, $vars->{AUTHORS} );  ## Authors, a comma separated list of Authors with complete First followed by Last Name
	$in_worksheet1->AddCell( $row, $col+27, "");  ## Comment
	$in_worksheet1->AddCell( $row, $col+28, "");  ## Comment Icon
	my ($virus,$country,$identifier,$year) = split /\//, $vars->{VIR_NAME};
	my $ncbi_virus_name = "SARS-CoV-2/Homo sapiens/$country/$identifier/$year";
	my $src_tsv_string = join("\t",$seqID,"Severe acute respiratory syndrome coronavirus 2",$ncbi_virus_name,$vars->{SM_CDATE},$country,$vars->{SM_HOST},$vars->{ORIG_LAB});
	push @src_tsv_content, $src_tsv_string;
	my $cmt_tsv_string = join("\t",$seqID,"Assembly-Data",$vars->{ASM_METHOD},$vars->{SM_COV},$vars->{SM_SEQUENCING_TECH},"Assembly-Data");
	push @cmt_tsv_content, $cmt_tsv_string;
	$row++;
}

$template->SaveAs($out);
&write_tsv($source_tsvout, join("\t",@src_tsv_header), join("\n",@src_tsv_content));
&write_tsv($comment_tsvout,join("\t",@cmt_tsv_header), join("\n",@cmt_tsv_content));

sub write_tsv{
	my $outfile = shift;
	my $header = shift;
	my $content = shift;
	open (my $ofh, ">", $outfile) or die "Cannot write to $outfile\n";
	print $ofh join("\n",$header,$content);
	close $ofh;
}

sub write_all_sequences {
	my $outfile = shift;
	my $vars= shift;
	my $seqID = shift;
	my $con_fasta = $vars->{SM_COV_FILE};
	open (my $ofh, ">>", $outfile ) or die "Cannot write to $outfile\n";
	open (my $ifh, "<", $con_fasta ) or die "Cannot read $con_fasta\n";
	while (<$ifh>){
		if(/^>/){
			if ($seqID){
				print $ofh ">".$seqID."\n";
			}else{
				print $ofh ">". $vars->{VIR_NAME}."\n";
			}
		}else{
			print $ofh $_;
		}
	}	
	close $ofh;
	close $ifh;
}

sub pull_consensusFiles {
	my $file_dir = shift;
	my $refname;
	my @consensusFastaFiles = glob("$file_dir/ReadsBasedAnalysis/readsMappingToRef/*consensus.fasta");
	foreach my $file (@consensusFastaFiles){
		open (my $fh, "<", $file) or die "Cannot open $file\n";
		while(my $line=<$fh>){
			chomp $line;
			if ($line =~ /^>(\S+)_consensus_(\S+)/){
				$refname->{$2}=$file;
				last;
			}
		}
		close $fh;
	}
	return ($refname);
}

sub pull_consensusInfo{
	my $file_dir=shift;
	my $vars=shift;
	my $conf=shift;
	my $conFastaFile=&pull_consensusFiles($file_dir);
        # coverage info is at ReadsBasedAnalysis/readsMappingToRef/readsToRef.alnstats.txt
	foreach my $key (keys %$conFastaFile){
		if ($vars->{SM_COV} =~ /$key/ ){
			$vars->{SM_COV_FILE} = $conFastaFile->{$key};
		}
	}
        open (my $log_fh,"<", "$file_dir/ReadsBasedAnalysis/readsMappingToRef/mapping.log");
        my ($tool_version, $align_tool);
        while(<$log_fh>){
                if (/Version:\s+(\S+)/){$tool_version=$1;}
                if (/CMD:\s(\S+)/){$align_tool=$1;}
        }
        close $log_fh;
        $vars->{ASM_METHOD}="EDGE-covid19: $align_tool $tool_version.";

        my $con_min_mapQ = $conf->{r2g_consensus_min_mapQ};
        my $con_min_cov = $conf->{r2g_consensus_min_cov};
        my $con_alt_prop = $conf->{r2g_consensus_alt_prop}*100 . "%";
        my $con_alt_indel = $conf->{r2g_consensus_altIndel_prop}* 100 . "%";

        $vars->{ASM_METHOD} .= " Consensus min coverage: ${con_min_cov}X. min map quality: $con_min_mapQ. Alternate Base > $con_alt_prop. Indel > $con_alt_indel.";
	my $platform = $conf->{FASTQ_SOURCE} || "";
        $vars->{FASTQ_SOURCE}=($platform eq "nanopore")?"Nanopore":"Iumina";
}

sub pull_sampleMetadata {
	my $file_dir=shift;
	my $vars=shift;
        my $metadata = "$file_dir/metadata_gisaid_ncbi.txt";
        if(-e $metadata) {
                open CONF, $metadata or die "Can't open $metadata $!";
                while(<CONF>){
                        chomp;
                        next if(/^#/);
                        if ( /(.*)=(.*)/ ){
                                $vars->{VIR_NAME} =$2 if ($1 eq "virus_name");
                                $vars->{VIR_PASSAGE} =$2 if ($1 eq "virus_passage");
                                $vars->{SM_CDATE} =$2 if ($1 eq "collection_date");
                                $vars->{SM_LOC} =$2 if ($1 eq "location");
                                $vars->{SM_HOST} =$2 if ($1 eq "host");
                                $vars->{SM_GENDER} = $2 if ($1 eq "gender");
                                $vars->{SM_AGE} =$2 if ($1 eq "age");
                                $vars->{SM_STATUS} =$2 if ($1 eq "status");
                                $vars->{SM_SEQUENCING_TECH} =$2 if ($1 eq "sequencing_technology");
                                $vars->{SM_COV} = $2 if ($1 eq "coverage");
                        }
                }
                close CONF;
        }
}
sub pull_submissionData {
	my $file_dir=shift;
	my $vars=shift;
        my $metadata = "$userDir/gisaid_ncbi_submission_profile.txt";
        if(-e $metadata) {
                open CONF, $metadata or die "Can't open $metadata $!";
                while(<CONF>){
                        chomp;
                        next if(/^#/);
                        if ( /(.*)=(.*)/ ){
                                $vars->{ORIG_LAB} =$2 if ($1 eq "originating_lab");
                                $vars->{ORIG_ADDRESS} =$2 if ($1 eq "originating_address");
                                $vars->{SUB_LAB} =$2 if ($1 eq "submitting_lab");
                                $vars->{SUB_ADDRESS} =$2 if ($1 eq "submitting_address");
                                $vars->{AUTHORS} =$2 if ($1 eq "authors");
                                $vars->{SUBMITTER} =$2 if ($1 eq "submitter");
                                $vars->{ID} =$2 if ($1 eq "gisaid_id");
                        }
                  }
                close CONF;
        }
}
sub getParams {
        my $config = shift;
        my $sys;

	if(-e $config) {
		open CONF, $config or die "Can't open $config: $!";
		while(<CONF>){
	      		chomp;
		        next if(/^#/);
		   	if ( /(.*)=(.*)/ ){
		     		$sys->{$1}=$2;
		      	}
		}
		close CONF;
	}
        return $sys;
}

1;
