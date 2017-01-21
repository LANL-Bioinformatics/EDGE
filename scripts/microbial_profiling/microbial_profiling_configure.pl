#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;
use JSON;
use FindBin qw($Bin);
use lib "$Bin/../../lib";
use HTML::Template;
#use DataG::Dumper;

my %opt;
my $EDGE_HOME = "$Bin/../..";
GetOptions(\%opt,
           "template=s",
           "tools=s",
           "bwaScoreCut=i",
           "configJson=s",
           "bwa-db=s",
           "metaphlan-db=s",
           "kraken-db=s",
           "gottcha-v-speDB=s",
           "gottcha-b-speDB=s",
           "gottcha-v-strDB=s",
           "gottcha-b-strDB=s",
           "gottcha-v-genDB=s",
           "gottcha-b-genDB=s",
           "gottcha2-v-genDB=s",
           "gottcha2-b-speDB=s",
           "gottcha2-v-speDB=s",
           "gottcha2-e-invDB=s",
           "gottcha2-e-ptzDB=s",
           "gottcha2-e-ptgDB=s",
           "pangia-db=s",
           'help|h|?'
          );
if ( $opt{help} || scalar keys %opt == 0 ) { &usage(); }

# setting up default values
$opt{'template'}        = $opt{'template'}||$ARGV[0]||"$Bin/microbial_profiling.settings.tmpl";
$opt{'tools'}           ||= $ARGV[1];
$opt{'bwaScoreCut'}     ||= 30;
$opt{'bwa-db'}          ||= "$EDGE_HOME/database/bwa_index/NCBI-Bacteria-Virus.fna";
$opt{"metaphlan-db"}    ||= "$EDGE_HOME/database/metaphlan/mpa";
$opt{"kraken-db"}       ||= "$EDGE_HOME/database/minikraken/database.idx";
$opt{"gottcha-v-speDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.species";
$opt{"gottcha-b-speDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species";
$opt{"gottcha-v-strDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain";
$opt{"gottcha-b-strDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.strain";
$opt{"gottcha-v-genDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.genus";
$opt{"gottcha-b-genDB"} ||= "$EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.genus";
$opt{"gottcha2-v-genDB"} ||= "$EDGE_HOME/database/GOTTCHA2/Virus_VIPR_HIVdb_IRD_NCBI_xHuBaAr_noEngStv.genus.fna";
$opt{"gottcha2-b-speDB"} ||= "$EDGE_HOME/database/GOTTCHA2/RefSeq-Release75.Bacteria.species.fna";
$opt{"gottcha2-v-speDB"} ||= "$EDGE_HOME/database/GOTTCHA2/Virus_VIPR_HIVdb_IRD_NCBI_xHuBaAr_noEngStv.species.fna";
$opt{"gottcha2-e-invDB"} ||= "$EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.invertebrate.species.fna";
$opt{"gottcha2-e-ptgDB"} ||= "$EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.pathogen.species.fna";
$opt{"gottcha2-e-ptzDB"} ||= "$EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.protozoa.species.fna";

#PanGIA configs
my $config_json = readListFromJson($opt{"configJson"});
$opt{"pangia-db"} = $config_json->{"edge-taxa-pangia-db"} || "$EDGE_HOME/database/PanGIA/NCBI_genomes_111216_p_GRCh38.fa";
$opt{"pangia-bg"} = "$EDGE_HOME/database/PanGIA/background/$config_json->{'edge-taxa-pangia-bg'}" if $config_json->{'edge-taxa-pangia-bg'};
$opt{"pangia-ra"} = $config_json->{"edge-taxa-pangia-ra"} || "READ_COUNT";
$opt{"pangia-ms"} = $config_json->{"edge-taxa-pangia-ms"} || "0" ;
$opt{"pangia-mr"} = $config_json->{"edge-taxa-pangia-mr"} || "3" ;
$opt{"pangia-mb"} = $config_json->{"edge-taxa-pangia-mb"} || "1" ;
$opt{"pangia-ml"} = $config_json->{"edge-taxa-pangia-ml"} || "50";
$opt{"pangia-rc"} = $config_json->{"edge-taxa-pangia-rc"} || "R_MAT";
$opt{"pangia-opts"} = "";
$opt{"pangia-opts"} = "-ps" if $config_json->{"edge-taxa-pangia-ps-sw"};

my @tools = split /,/, $opt{'tools'};
print STDERR &usage() if !-e $opt{'template'};
print STDERR &usage() if scalar @tools < 1;

# remove suffix if any to meet the tool db format
foreach my $db ( keys %opt) { 
	$opt{"$db"} =~ s/\.?(amb|ann|bwt|fai|pac|sa|parsedGOTTCHA\.dmp)$// if ($db =~ /bwa-db|gottcha.*DB|pangia/ );
	$opt{"$db"} =~ s/(\.rev)?.\d\.bt2$// if ($db =~ /metaphlan-db/);
	(my $tmp_filenaem,$opt{"kraken-db"},my $tmp_suffix)=&fileparse($opt{"kraken-db"}) if ($db =~ /kraken-db/);
	
}

foreach my $tool ( @tools ){
	next unless $tool;
	$opt{"$tool"} = 1;
}

#print Dumper (\%opt);
my $template = HTML::Template->new(filename => $opt{'template'},  die_on_bad_params => 0  );
#my $template = HTML::Template->new(filename => $opt{'template'}, die_on_bad_params => 1 );
$template->param( %opt );
print $template->output();


sub readListFromJson {
    my $json = shift;
    my $list = {};
    if( -r $json ){
         open JSON, $json;
         flock(JSON, 1);
         local $/ = undef;
         $list = decode_json(<JSON>);
         close JSON;
    }
    return $list;
}

sub usage {
print <<USAGE;
$0 [template.tmpl] [tools] > microbial_profiling_configure.settings.ini 
	
    [Options]           [Default]
    -------------------------------------------------------------------
    -template           $Bin/microbial_profiling.settings.tmpl
    -tools              comma separated tools: bwa,kraken-mini, ...
    -bwaScoreCut        minimum score to output for BWA [30]
    -bwa-db             $EDGE_HOME/database/bwa_index/NCBI-Bacteria-Virus.fna
    -metaphlan-db       $EDGE_HOME/database/metaphlan/mpa
    -kraken-db          $EDGE_HOME/database/minikraken
    -gottcha-v-speDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.species
    -gottcha-b-speDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species
    -gottcha-v-strDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain
    -gottcha-b-strDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.strain
    -gottcha-v-genDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.genus
    -gottcha-b-genDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.genus
    -gottcha2-v-genDB   $EDGE_HOME/database/GOTTCHA2/Virus_VIPR_HIVdb_IRD_NCBI_xHuBaAr_noEngStv.genus.fna
    -gottcha2-b-speDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release75.Bacteria.species.fna
    -gottcha2-v-speDB   $EDGE_HOME/database/GOTTCHA2/Virus_VIPR_HIVdb_IRD_NCBI_xHuBaAr_noEngStv.species.fna
    -gottcha2-e-invDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.invertebrate.species.fna
    -gottcha2-e-ptzDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.protozoa.species.fna
    -gottcha2-e-ptgDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-release75.Euk_only.pathogen.species.fna

USAGE
    #
exit;
} 

exit 0;

