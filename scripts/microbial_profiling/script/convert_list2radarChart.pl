#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../../../lib";
use HTML::Template;

$|=1;
$ENV{PATH} = "$Bin:$Bin/bin/:$Bin/scripts/:$ENV{PATH}";

my $template;
my %opt;
my $res=GetOptions(\%opt,
                   'level|l=s',
                   'outdir|o=s',
                   'outprefix|p=s',
                   'filter_taxa|f=s',
                   'title=s',
                   'template|t=s',
                   'top=i',
                   'help|h|?') || &usage();

&usage() if ( defined $opt{help} || scalar @ARGV < 1 );

$opt{template}    ||= "$Bin/convert_list2radarChart.html.tmpl";
$opt{outdir}      ||= ".";
$opt{outprefix}   ||= "";
$opt{top}         ||= 5;
$opt{filter_taxa} ||= "";
$opt{level}       ||= "genus";

if( -r $opt{template}){
	$template = HTML::Template->new(filename => $opt{template});
}
else{
	die "ERROR: Can't find template: $opt{template}.";
}

my %filter_taxa = map { $_ => 1 } split /,/,$opt{filter_taxa};

my @listfile = @ARGV;
my $matrix;
my $count;
my $order;
my $data;
my $tools;
my %tool_order=(
	'blastn'          => 10,
	'bwa'             => 20,
	'bwa-target'      => 21,
	'metaphlan'       => 22,
	'metaphlan2'      => 23,
	'motus'           => 24,
	'kraken_mini'     => 25,
	'kraken'          => 26,
	'kraken2'         => 27,
	'centrifuge'      => 28,
	'gottcha-genDB-b' => 30,
	'gottcha-speDB-b' => 40,
	'gottcha2-speDB-b'=> 41,
	'gottcha-strDB-b' => 50,
	'gottcha-genDB-v' => 60,
	'gottcha-speDB-v' => 70,
	'gottcha2-speDB-v'=> 71,
	'gottcha-strDB-v' => 80,
	'sequedex-opg'    => 90,
	'sequedex-tol'    => 100,
	'metacv'          => 120,
	'metaphyler-bn'   => 130,
	'metaphyler-bx'   => 140,
	'metaphyler-srv'  => 150,
	'pangia'          => 160,
	'diamond'         => 170
);

############################
# parsing list files
############################

foreach my $file ( @listfile ){
	my ($dataset,$tool) = $file =~ /\d+_([^\/]+)\/([^\/]+)\/([^\/]+)\.[^\.]+$/;

	if( $opt{level} =~ /species/ ){
		next if $tool =~ /gottcha-genDB/;
	}	
	elsif( $opt{level} =~ /strain/ ){
		next if $tool =~ /gottcha-genDB/;
		next if $tool =~ /gottcha-speDB/;
	}

	$tools->{$tool}=1;
	my $filtered_taxa=0;

	open LIST, $file or die "Can't open $file\n";
	while(<LIST>){
		chomp;
		my @fields = split /\t/, $_;
		next if $fields[0] ne $opt{level};
        next if defined $filter_taxa{$fields[1]};

		$matrix->{$dataset}->{$tool}->{$fields[1]} = $fields[2];

		# counting total reads/abundance
        $count->{$dataset}->{$tool} ||= 0;
		$count->{$dataset}->{$tool} += $fields[2];
    }
    close LIST;
}

#############################
# Create orders
#############################

my $top_taxa;

foreach my $dataset ( keys %$matrix ){
	foreach my $tool ( keys %{$matrix->{$dataset}} ){
        my $cnt=0;
        my $top_sum=0;
	foreach my $taxa ( sort {
                        $matrix->{$dataset}->{$tool}->{$b} <=> $matrix->{$dataset}->{$tool}->{$a} 
                    } keys %{$matrix->{$dataset}->{$tool}} ){
        	last if $cnt == $opt{top};
        	$top_taxa->{$dataset}->{$taxa}=1;
	        $top_sum += $matrix->{$dataset}->{$tool}->{$taxa};
        	$cnt++;
	}
        $matrix->{$dataset}->{$tool}->{others} = $count->{$dataset}->{$tool} - $top_sum;
    }
}

#############################
# Generating matrix
#############################

my ($legend_text, $data_text);

foreach my $dataset ( keys %$top_taxa ) {
    my @legend = sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools;
    my @data_tool;

    print STDERR "Generate chart: $dataset\n";

    foreach my $tool ( sort { $tool_order{$a} <=> $tool_order{$b} } keys %$tools ){
        my @data;
        my @top_taxas = sort { $matrix->{$dataset}->{bwa}->{$a} <=> $matrix->{$dataset}->{bwa}->{$b}
                                || $matrix->{$dataset}->{'gottcha-genDB-b'}->{$a} <=> $matrix->{$dataset}->{'gottcha-genDB-b'}->{$b} }
                            keys %{$top_taxa->{$dataset}};
	    foreach my $taxa ( ("others",@top_taxas) ){
            #next unless defined $matrix->{$dataset}->{$tool}->{$taxa};
            my $val = $matrix->{$dataset}->{$tool}->{$taxa};
            my $tol = $count->{$dataset}->{$tool};
            $tol = 1 if $val==0 && $tol==0;
            my $pct = sprintf "%.4f", $val/$tol;
            push @data, "{axis:'$taxa',value: $pct}";
            print "\t$tool\t$taxa\t$val\t$pct\n";
        }
        my $temp = join ",\n", @data;
        push @data_tool, "[ $temp ]";
    }
    
    $legend_text = join "','", @legend;
    $legend_text = "'$legend_text'";

    $data_text = join ",\n", @data_tool;

    $template->param(TITLE  => $opt{title} ? $opt{title} : $dataset);
    $template->param(LEGEND => $legend_text);
    $template->param(DATA   => $data_text);

    $opt{outprefix} = $opt{outprefix}."_" if $opt{outprefix};
    open OUTPUT, ">$opt{outdir}/$opt{outprefix}$dataset.$opt{level}.html";
    print OUTPUT $template->output;
    close OUTPUT;
}

sub usage {
	print STDERR <<__USAGE__;
USAGE: $0 <DATASET>/<TOOL>/<DATASET>.out.list [OPTIONS]

EXAMPLE: 
    $0 1_CDC_10_clean/bwa/*.out.list 1_CDC_10_clean/metaphlan/*.out.list

    $0 \\
		--level genus \\
		--top 5 \\
		--template convert_list2radarChart.html.tmpl \\
		1_454_EVEN/*/*.out.list

OPTIONS:

  --level <RANK>

      Default is "species". Any ranks from "kingdom" to "strain" are allowable.

  --top (INT)

      Only display a certain number of classifications.

  --exclude <TAXA1(,TAXA2,...)>
  --outdir <DIR>
  --outprefix <STRING>
  --template <HTML TEMPLATE>
      
  --help

__USAGE__
	exit;
}
