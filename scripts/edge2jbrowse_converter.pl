#!/usr/bin/env perl
# contig_classifier.pl
# ver 0.1
# 2014/04/08
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

use strict;
use Getopt::Long;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$ENV{EDGE_HOME}/site/JBrowse/bin:/$ENV{PATH}";
$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'proj_outdir=s',             #./
    # generate GFF3 ---------------------------------------------------------------------------------------
    'in-phage-finder=s',         #AssemblyBasedAnalysis/Prophage/phageFinder_summary.txt 
    'in-ref-gff3=s',             
    'in-ctg2ref-coords=s',       #AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.coords
    'in-ctg2ref-snps=s',         #AssemblyBasedAnalysis/Indels/contigsToRef.SNPs_report.txt
    'in-ctg2ref-indels=s',       #AssemblyBasedAnalysis/Indels/contigsToRef.Indels_report.txt
    'gff3out=s',                 #edge_analysis.gff3
    # generate JBrowse tracks -----------------------------------------------------------------------------
    'in-ctg-fa=s',               #AssemblyBasedAnalysis/contigs.fa
    'in-ref-fa=s',               
    'in-ctg-anno-gff3=s',        #AssemblyBasedAnalysis/Annotation/PROKKA.gff
    'in-ctg-adj-primer-gff3=s',  #AssayCheck/PCR.Adjudication.primers.gff3
    'in-read2ctg-bam=s',         #AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.sort.bam
    'in-read2ref-bam=s',         #ReadsBasedAnalysis/readsToRef.sort.bam
    'in-primer2ctg-bam=s',       #AssayCheck/pcrContigValidation.bam
    'in-primer2ref-bam=s',       #AssayCheck/pcrRefValidation.bam
    'in-read2ref-vcf=s',         #ReadsBasedAnalysis/readsToRef.vcf 
    'outdir=s',                  #JBrowse
    'ctg-coord-conf=s',          #\$EDGE_HOME/script/edge2jbrowse_converter.ctg_coord_conf
    'ctg-coord-bam-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'ctg-coord-primer-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'ref-coord-conf=s',          #\$EDGE_HOME/script/edge2jbrowse_converter.ref_coord_conf
    'ref-coord-bam-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'ref-coord-vcf-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'ref-coord-primer-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'ref-coord-bw-conf=s',      #\$EDGE_HOME/script/edge2jbrowse_converter.ctg-coord-ext-conf
    'out-ctg_coord-dir=s',       #[outdir]/jb_ctg_tracks/
    'out-ref-coord-dir=s',       #[outdir]/jb_ref_tracks/
    'help|h|?'
) || &usage();

if ( $opt{help} || scalar keys %opt == 0 ) { &usage(); }

# setting up default values
$opt{'proj_outdir'}            |= ".";
$opt{'in-phage-finder'}        |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/Prophage/phageFinder_summary.txt";
$opt{'in-ctg2ref-coords'}      |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.coords";
$opt{'in-ctg2ref-snps'}        |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/Indels/contigsToRef.SNPs_report.txt";
$opt{'in-ctg2ref-indels'}      |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/Indels/contigsToRef.Indels_report.txt";
$opt{'in-ctg-fa'}              |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/contigs.fa";
$opt{'in-ctg-anno-gff3'}       |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/Annotation/PROKKA.gff";
$opt{'in-ctg-adj-primer-gff3'} |= "$opt{'proj_outdir'}/AssayCheck/PCR.Adjudication.primers.gff3";
$opt{'in-read2ctg-bam'}        |= "$opt{'proj_outdir'}/AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.sort.bam";
$opt{'in-read2ref-bam'}        |= "$opt{'proj_outdir'}/ReadsBasedAnalysis/readsToRef.sort.bam";
$opt{'in-read2ref-vcf'}        |= "$opt{'proj_outdir'}/ReadsBasedAnalysis/readsToRef.vcf";
$opt{'in-primer2ctg-bam'}      |= "$opt{'proj_outdir'}/AssayCheck/pcrContigValidation.bam";
$opt{'in-primer2ref-bam'}      |= "$opt{'proj_outdir'}/AssayCheck/pcrRefValidation.bam";
$opt{'outdir'}                 |= "$opt{'proj_outdir'}/JBrowse";
$opt{'ctg-coord-conf'}         |= "$Bin/edge2jbrowse_converter.ctg_conf.temp";
$opt{'ctg-coord-bam-conf'}     |= "$Bin/edge2jbrowse_converter.ctg_trackList_BAM.temp";
$opt{'ctg-coord-primer-conf'}  |= "$Bin/edge2jbrowse_converter.ctg_trackList_PCR.temp";
$opt{'ref-coord-conf'}         |= "$Bin/edge2jbrowse_converter.ref_conf.temp";
$opt{'ref-coord-bam-conf'}     |= "$Bin/edge2jbrowse_converter.ref_trackList_BAM.temp";
$opt{'ref-coord-vcf-conf'}     |= "$Bin/edge2jbrowse_converter.ref_trackList_VCF.temp";
$opt{'ref-coord-primer-conf'}  |= "$Bin/edge2jbrowse_converter.ref_trackList_PCR.temp";
$opt{'ref-coord-bw-conf'}      |= "$Bin/edge2jbrowse_converter.ref_trackList_BW.temp";
$opt{'gff3out'}                |= "$opt{'outdir'}/edge_analysis.gff3";
$opt{'out-ctg-coord-dir'}      |= "$opt{'outdir'}/ctg_tracks";
$opt{'out-ref-coord-dir'}      |= "$opt{'outdir'}/ref_tracks";

&main();

##########################################################################
#
# MAIN
#
#########################################################################

sub main {
	#generate output directory
	executeCommand("mkdir -p $opt{outdir}");
	die "Can't create output directory: $!\n" unless ( -d $opt{outdir} );

	#converting analysis to GFF3
	print "\n# Converting analysis to GFF3:\n";
	&generateEdgeAnalysisGff3( $opt{gff3out} );
	
	#[Contigs as coordinate]
	#
	#AssayCheck/PCR.Adjudication.primers.gff3
	#AssemblyBasedAnalysis/Annotation/PROKKA.gff
	#AssemblyBasedAnalysis/Prophage/phageFinder_summary.txt (need to convert to gff3)
	#AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.sort.bam

	#preparing files for CONTIG-based tracks
	if( -e $opt{'in-ctg-fa'} ){
		print "\n# Preparing files for JBrowse with CONTIG-based tracks:\n";
		executeCommand("rm -rf $opt{'out-ctg-coord-dir'}");
		executeCommand("mkdir -p $opt{'out-ctg-coord-dir'}");
		executeCommand("cat $opt{gff3out} >> $opt{'out-ctg-coord-dir'}/edge_analysis_merged.gff3");

		if( -e $opt{'in-ctg-anno-gff3'} ){
			print "#  - amending annotation gff3...";
			&fixProkkaGff( $opt{'in-ctg-anno-gff3'}, "$opt{outdir}/contig_annotation.gff3" );
			executeCommand("cat $opt{outdir}/contig_annotation.gff3 >> $opt{'out-ctg-coord-dir'}/edge_analysis_merged.gff3");
			print "Done.\n";
		}

		if( -e $opt{'in-ctg-adj-primer-gff3'} ){
			print "#  - adding PcrAdjudicationPrimers gff3...";
			executeCommand("cat $opt{'in-ctg-adj-primer-gff3'} >> $opt{'out-ctg-coord-dir'}/edge_analysis_merged.gff3");
			print "Done.\n";
		}	

		print "\n# Generating CONTIG-based tracks:\n";
		
		#generate CONTIG-based tracks
		print "#  - adding feature tracks...";
		my $esc_str = $opt{'out-ctg-coord-dir'};
		$esc_str =~ s/\//\\\//g;
		executeCommand("cat $opt{'ctg-coord-conf'} | perl -pi -e \"s/%%GFF3DIR%%/$esc_str/\" > $opt{'out-ctg-coord-dir'}/edge2jbrowse_converter.ctg_conf");
		executeCommand("prepare-refseqs.pl --fasta $opt{'in-ctg-fa'} --out $opt{'out-ctg-coord-dir'}");
		executeCommand("biodb-to-json.pl --quiet --conf $opt{'out-ctg-coord-dir'}/edge2jbrowse_converter.ctg_conf --out $opt{'out-ctg-coord-dir'}");
		print "Done.\n";

		#add read2ctg BAM track
		if( -e $opt{'in-read2ctg-bam'} ){
			print "#  - Adding read2ctg BAM track...";
			executeCommand("samtools view -F4 -h $opt{'in-read2ctg-bam'} | samtools view -bS - > $opt{'out-ctg-coord-dir'}/readsToContigs.mapped.bam");
			executeCommand("samtools sort $opt{'out-ctg-coord-dir'}/readsToContigs.mapped.bam $opt{'out-ctg-coord-dir'}/readsToContigs.mapped.sort");
			executeCommand("samtools index $opt{'out-ctg-coord-dir'}/readsToContigs.mapped.sort.bam");
			executeCommand("add-track-json.pl $opt{'ctg-coord-bam-conf'} $opt{'out-ctg-coord-dir'}/trackList.json");
			print "Done.\n";
		}


		#add pcrValidation track
		if( -e $opt{'in-primer2ctg-bam'} ){
			print "#  - Adding pcrValidation track...";
			executeCommand("samtools sort $opt{'in-primer2ctg-bam'} $opt{'out-ctg-coord-dir'}/pcrContigValidation.sort");
			executeCommand("samtools index $opt{'out-ctg-coord-dir'}/pcrContigValidation.sort.bam");
			executeCommand("add-track-json.pl $opt{'ctg-coord-primer-conf'} $opt{'out-ctg-coord-dir'}/trackList.json");
			print "Done.\n";
		}
		
		print "#  - Indexing features...";
		executeCommand("generate-names.pl --out $opt{'out-ctg-coord-dir'}");
		print "Done.\n";
	}

	#[Reference as coordinate]
	#
	#Reference/NC_000913.fna
	#Reference/NC_000913.gbk (genbank2gff3.pl is in EDGE/scripts/)
	#AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.coords   (need to convert to gff3)
	#ReadsBasedAnalysis/readsToRef.sort.bam
	#ReadsBasedAnalysis/readsToRef.vcf
	#AssemblyBasedAnalysis/Indels/contigsToRef.SNPs_report.txt   (need to convert to gff3)
	#AssemblyBasedAnalysis/Indels/contigsToRef.Indels_report.txt    (need to convert to gff3)

	#generate jb tracks using REF as coordinates
	if( -e $opt{'in-ref-fa'} ){
		print "\n# Preparing files for JBrowse with REF-based tracks:\n";
		executeCommand("rm -rf $opt{'out-ref-coord-dir'}/");
		executeCommand("mkdir -p $opt{'out-ref-coord-dir'}/");
		executeCommand("cat $opt{gff3out} >> $opt{'out-ref-coord-dir'}/edge_analysis_merged.gff3");
		
		if( -e $opt{'in-ref-gff3'} ){
			print "#  - Amending annotation gff3...";
			&fixRefGff3( $opt{'in-ref-gff3'}, $opt{'in-ref-fa'}, "$opt{outdir}/ref_annotation.gff3");
			executeCommand("cat $opt{outdir}/ref_annotation.gff3 >> $opt{'out-ref-coord-dir'}/edge_analysis_merged.gff3");
			print "Done.\n";
		}
		
		#generate CONTIG-based tracks
		print "\n# Generating CONTIG-based tracks:\n";
		
		print "#  - preparing sequence track...";
		executeCommand("prepare-refseqs.pl --fasta $opt{'in-ref-fa'} --out $opt{'out-ref-coord-dir'}");
		print "Done.\n";
		
		print "#  - adding feature tracks...";
		my $esc_str = $opt{'out-ref-coord-dir'};
		$esc_str =~ s/\//\\\//g;
		executeCommand("cat $opt{'ref-coord-conf'} | perl -pi -e \"s/%%GFF3DIR%%/$esc_str/\" > $opt{'out-ref-coord-dir'}/edge2jbrowse_converter.ref_conf");
		executeCommand("biodb-to-json.pl --quiet --conf $opt{'out-ref-coord-dir'}/edge2jbrowse_converter.ref_conf --out $opt{'out-ref-coord-dir'}");
		print "Done.\n";
	
		#add pcrValidation track
		if( -e $opt{'in-primer2ref-bam'} ){
			print "#  - Adding pcrValidation track...";
			executeCommand("samtools sort $opt{'in-primer2ref-bam'} $opt{'out-ref-coord-dir'}/pcrRefValidation.sort");
			executeCommand("samtools index $opt{'out-ref-coord-dir'}/pcrRefValidation.sort.bam");
			executeCommand("add-track-json.pl $opt{'ref-coord-primer-conf'} $opt{'out-ref-coord-dir'}/trackList.json");
			print "Done.\n";
		}

		if( -e $opt{'in-read2ref-bam'} ){
			print "#  - Adding read2ref BAM track...";
			executeCommand("samtools view -F4 -h $opt{'in-read2ref-bam'} | samtools view -bS - > $opt{'out-ref-coord-dir'}/readsToRef.mapped.bam");
	        	executeCommand("samtools sort $opt{'out-ref-coord-dir'}/readsToRef.mapped.bam $opt{'out-ref-coord-dir'}/readsToRef.mapped.sort");
	        	executeCommand("samtools index $opt{'out-ref-coord-dir'}/readsToRef.mapped.sort.bam");
			executeCommand("add-track-json.pl $opt{'ref-coord-bam-conf'} $opt{'out-ref-coord-dir'}/trackList.json");
			print "Done.\n";

			print "#  - Adding BigWig track...";
			executeCommand("convert_bam2bigwig.pl $opt{'out-ref-coord-dir'}/readsToRef.mapped.sort.bam");
			executeCommand("add-track-json.pl $opt{'ref-coord-bw-conf'} $opt{'out-ref-coord-dir'}/trackList.json");
			print "Done.\n";

			print "#  - Adding read2ref VCF track...";
			executeCommand("bgzip -c $opt{'in-read2ref-vcf'} > $opt{'out-ref-coord-dir'}/readsToRef.vcf.gz");
   			executeCommand("tabix -p vcf $opt{'out-ref-coord-dir'}/readsToRef.vcf.gz");
			executeCommand("add-track-json.pl $opt{'ref-coord-vcf-conf'} $opt{'out-ref-coord-dir'}/trackList.json");
			print "Done.\n";
		}
		
		print "#  - Indexing features...";
		executeCommand("generate-names.pl --out $opt{'out-ref-coord-dir'}");	
		print "Done.\n";
	}

	print "\n# ALL DONE.\n";
}

##########################################################################
#
# SUBROUTINES
#
#########################################################################

#GFF3 format
#Column 1: "seqid"
#Column 2: "source"
#Column 3: "type"
#Columns 4 & 5: "start" and "end"
#Column 6: "score"
#Column 7: "strand"
#Column 8: "phase"
#Column 9: "attributes"

sub generateEdgeAnalysisGff3 {
	my $outfile = shift;
	open GFF3OUT, ">$outfile" or die "$!";
	
	my $out_text;
	if( defined $opt{"in-phage-finder"} ){
		print "#  - Parsing PhageFinder result...";
		$out_text = &convertPhageFinder( $opt{"in-phage-finder"} );
		print GFF3OUT $out_text;
		print "DONE.\n";
	}
	if( defined $opt{"in-ctg2ref-coords"} ){
		print "#  - Parsing contig to reference mapping result...";
		$out_text = &convertCtg2refCoords( $opt{"in-ctg2ref-coords"} );
		print GFF3OUT $out_text;
		print "DONE.\n";
	}
	if( defined $opt{"in-ctg2ref-snps"} ){
		print "#  - Parsing SNPs result...";
		$out_text = &convertCtg2refSnps( $opt{"in-ctg2ref-snps"} );
		print GFF3OUT $out_text;
		print "DONE.\n";
	}
	if( defined $opt{"in-ctg2ref-indels"} ){
		print "#  - Parsing InDels result...";
		$out_text = &convertCtg2refIndels( $opt{"in-ctg2ref-indels"} );
		print GFF3OUT $out_text;
		print "DONE.\n";
	}
	
	close GFF3OUT;
}

sub convertPhageFinder {
	my $file = shift;
	my @header;
	my $id=1;
	my $out_text="";
	open IN, $file or die "Can't open PhageFinder output $file. $!\n";
	while(<IN>){
		chomp;
		if(/^#/){
			@header = split /\t/, $_;
			next;
		}

		my @temp = split /\t/, $_;
		my $attr;

		for( my $i=0; $i<=$#temp; $i++){
			push @{$attr->{$header[$i]}}, $temp[$i];
		}
		
		#add ID attribute
		push @{$attr->{ID}}, "phage_finder_".$id++;
		push @{$attr->{Name}}, "$temp[6]_$temp[7]";

		my $attr_str = &gff_attr_string($attr);
		
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t+\t.\t%s\n",
				$temp[0],
				"phage_finder",
				"region",
				$temp[4],
				$temp[5],
				$attr_str;
		
	}
	close IN;

	return $out_text;
}

sub fixRefGff3 {
	my ($gff3,$fa,$gff3out) = @_;
	my $out_text;
	
	if( -r $gff3 && -r $fa ){
		my @fname = `grep '>' $fa | awk -F\\  '{print \$1}'`;
		my @gname = `awk -F\\\\t '{print \$1}' $gff3 | uniq`;
		
		my %name_mapping;
		foreach my $gn ( @gname ){
			chomp $gn;
			foreach my $fn ( @fname ){
				chomp $fn;
				$fn =~ s/^>//;
				if( $fn =~ /$gn/ ){
					$name_mapping{$gn} = $fn;
					last;
				}
			}
		}

		open IN, $gff3;
		while(<IN>){
			next if /^\s*$/;
			
			if( /^##FASTA/ ){
				$/ = ">";
				next;
			}
			
			#parsing FASTA section
			if( $/ eq ">" ){
				next if /^>$/;
				$_ =~ s/\>//g;
				my ($id,@seq) = split /\n/, $_;
				my $seq = join "", @seq;
				if( $seq =~ /^M/i ){ #protein sequence
					$out_text =~ s/ID=$id;([^\n]+)/ID=$id;$1;protein_seq=$seq/;
				}
			}
			#parsing GFF3 section
			else{
				chomp;
				next if /^#/;
				my ($refname) = $_ =~ /^(\S+)\t/;
				s/^$refname/$name_mapping{$refname}/ if defined $name_mapping{$refname};
				s/;product=([^;]+)/;Product=$1;Description=$1/i unless /;Description=/;
				$out_text .= "$_\n";
			}
		}
		$/ = "\n";
		close IN;
	}
	
	open OUT, ">$gff3out" or die "Can't write $gff3out: $!\n";
	print OUT $out_text;
	close OUT;
}

sub fixProkkaGff {
	my ($infile,$outfile) = @_;
	my $out_text;
	
	open IN, $infile or die "Can't run PROKKA GFF3: $!\n";

	while(<IN>){
		next if /^\s*$/;
		
		if( /^##FASTA/ ){
			$/ = ">";
			next;
		}
		
		#parsing FASTA section
		if( $/ eq ">" ){
			next if /^>$/;
			$_ =~ s/\>//g;
			my ($id,@seq) = split /\n/, $_;
			my $seq = join "", @seq;
			if( $seq =~ /^M/i ){ #protein sequence
				$out_text =~ s/ID=$id;([^\n]+)/ID=$id;$1;protein_seq=$seq/;
			}
		}
		#parsing GFF3 section
		else{
			chomp;
			next if /^#/;
			s/;product=([^;]+)/;Product=$1;Description=$1/i unless /;Description=/;
			if( /;gene=/i ){
				s/;gene=([^;]+)/;gene=$1;Name=$1/i;
			}
			else{
				s/\tID=([^;]+)/\tID=$1;Name=$1/i;
			}
			$out_text .= "$_\n";
		}
	}
	$/ = "\n";
	close IN;
	
	open OUT, ">$outfile" or die "Can't write to output file: $!";
	print OUT $out_text;	
	close OUT;
}

sub convertCtg2refCoords {
	my $file = shift;
	my $out_text = "";
	my $id=1;

	open IN, $file;

	#1->ref
	#2->query
	# 0        1       2       3       4       5       6        7        8     9        10     11     12
	#[S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [R_TAGS] [Q_TAGS]
	#refstr  refend  qrystr  qryend  aln_len                 tol_len         aln/tol 

	while(<IN>){
		chomp;
		next if /^\[/;
		my @temp = split /\t/, $_;
		next if scalar @temp != 13;
		my $descStr = &encode("contig_coverage: $temp[10]\%, identity: $temp[6]\%");
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n",
				$temp[11],
				"nucmer",
				"Ctg2Ref_coords",
				$temp[0],
				$temp[1],
				$temp[3]>$temp[2]?"+":"-",
				"ID=Ctg2Ref_$id;Name=$temp[12];Contig_start=$temp[2];Contig_end=$temp[3];Ref_aln_length=$temp[4];Ctg_aln_length=$temp[5];Identity=$temp[6];Contig_len=$temp[8];Ref_coverage=$temp[9];Contig_coverage=$temp[10];Description=$descStr";
		$id++;
	}
	close IN;

	return $out_text;
}

sub convertCtg2refSnps {
	my $file = shift;
	my $out_text = "";
	my $id=1;
	my @header;
	open IN, $file;
	while(<IN>){
		chomp;

		# [0]         [1]             [2]         [3]         [4]     [5]     [6]         [7]     [8]         [9]     [10]
		# Chromosome  SNP_position    Ref_codon   Sub_codon   aa_Ref  aa_Sub  Synonymous  Product CDS_start   CDS_end CDS_strand
	
		if ( $_ =~ /Ref_codon\tSub_codon/ ){
			@header = split /\t/, $_;
			next;
		}
		next if $_ =~ /Merged with SNP /;
	
		my @temp = split /\t/, $_;
		#next if scalar @temp < 11;
	
		# id generater
		$id++;
	
		# fill null
		$temp[4] ||= "-";
		$temp[5] ||= "-";
		$temp[6] ||= "-";
		$temp[8] ||= "-";
		$temp[9] ||= "-";
		$temp[10] = ($temp[10] == -1) ? "-":"+";
		
		my $attr;

		for( my $i=0; $i<=$#temp; $i++){
			push @{$attr->{$header[$i]}}, $temp[$i];
		}
		
		#add ID attribute
		push @{$attr->{ID}}, "SNP_$id";
		push @{$attr->{Name}}, "SNP_$id";
		push @{$attr->{Description}}, "$temp[2]($temp[4])->$temp[3]($temp[5])";

		my $attr_str = &gff_attr_string($attr);
	
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n",
				$temp[0],
				"contigsToRef",
				"SNPs",
				$temp[1],
				$temp[1]+length($temp[2])-1,
				$temp[10],
				$attr_str
		;
	
	}

	close IN;
	return $out_text;
}

sub convertCtg2refIndels {
	my $file = shift;
	my $out_text = "";
	my $id=1;
	my @header;
	open IN, $file;
	while(<IN>){
		chomp;

		# [0] Chromosome [1] INDEL_position [2] Sequence [3] Length [4] Type [5] Product [6] CDS_start [7] CDS_end
	
		if ( $_ =~ /Length\s+Type\s+Product/ ){
			@header = split /\t/, $_;
			next;
		}
	
		my @temp = split /\t/, $_;
		#next if scalar @temp < 8;
	
		# id generater
		$id++;
	
		my $attr;

		for( my $i=0; $i<=$#temp; $i++){
			push @{$attr->{$header[$i]}}, $temp[$i];
		}
		
		#add ID attribute
		push @{$attr->{ID}}, "indels_$id";
		push @{$attr->{Name}}, "indels_$id";
		push @{$attr->{Description}}, "$temp[4]:$temp[2] ($temp[3]bp);";

		my $attr_str = &gff_attr_string($attr);
	
		$out_text .= sprintf "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n",
				$temp[0],
				"contigsToRef",
				$temp[4],
				$temp[1],
				($temp[4] eq "Insertion") ? $temp[1] : ($temp[1]+$temp[3]-1),
				"+",
				$attr_str
		;
	
	}

	close IN;
	return $out_text;
}

sub gff_attr_string {
	my $attr = shift;
	my $has_id = 0;
	my $gff_str = "";

	foreach my $attrn ( keys %$attr ){
		my @attrv = map { &encode($_) } @{$attr->{$attrn}};
		my $temp = join ",", @attrv;

		$attrn =~ s/[<>#]//g;
		if( $attrn =~ /^id$/i ){
			$attrn = uc($attrn); 
			$gff_str = "$attrn=$temp;$gff_str";
		}
		else{
			$gff_str .= "$attrn=$temp;";
		}
	}
	
	return $gff_str;
}

sub encode {
	my $str = shift;
	$str =~ s/,/%2c/g;
	$str =~ s/;/%3b/g;
	return $str;
}

sub executeCommand 
{
    my $command = shift;
    system($command) == 0
         || die "the command $command failed.\n";
}



sub usage {

print <<__END__;

$0 [OPTIONS]
    
[OPTIONS]

    OPTION                    DEFAULT
    -----------------------------------------------------------------------------------------------
    --in-ref-fa               [no default]
    --in-ref-gff3             [no default]
    --proj_outdir             ./

    --in-phage-finder         AssemblyBasedAnalysis/Prophage/phageFinder_summary.txt 
    --in-ctg2ref-coords       AssemblyBasedAnalysis/contigMappingToRef/contigsToRef.coords
    --in-ctg2ref-snps         AssemblyBasedAnalysis/Indels/contigsToRef.SNPs_report.txt
    --in-ctg2ref-indels       AssemblyBasedAnalysis/Indels/contigsToRef.Indels_report.txt
    --in-contig-fa            AssemblyBasedAnalysis/contigs.fa
    --in-ctg-anno-gff3        AssemblyBasedAnalysis/Annotation/PROKKA.gff
    --in-ctg-adj-primer-gff3  AssayCheck/PCR.Adjudication.primers.gff3
    --in-read2ctg-bam         AssemblyBasedAnalysis/readsMappingToContig/readsToContigs.sort.bam
                              (.bai is required in the same directory)
    --in-read2ref-bam         ReadsBasedAnalysis/readsToRef.sort.bam
                              (.bai is required in the same directory)
    --in-read2ref-vcf         ReadsBasedAnalysis/readsToRef.vcf 
    --gff3out                 edge_analysis.gff3
    --outdir                  JBrowse
    --ctg-coord-conf          \$EDGE_HOME/script/edge2jbrowse_converter.ctg_coord_conf
    --ref-coord-conf          \$EDGE_HOME/script/edge2jbrowse_converter.ref_coord_conf
    --out-ctg_coord-dir       [outdir]/ctg_tracks/
    --out-ref-coord-dir       [outdir]/ref_tracks/
    --help|h|?

__END__
exit();
}
