#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin/../../../lib";
use Excel::Writer::XLSX;
use Getopt::Long;

$|=1;
my %opt;
my $res=GetOptions(\%opt,
                   'list|l=s',
                   'setting|s=s',
                   'highlight_list|h=s',
                   'output|o=s',
                   'help|?') || &usage();

if ( $opt{help} ) { &usage(); }

#read settings
my $ini_file = $opt{setting};
my $tools = &restore_settings( $ini_file );

#default settings
my $extract  = $tools->{system}->{EXTRACT_NUM};
my $p_outdir = $tools->{system}->{OUTDIR};
$p_outdir    = $opt{output} if defined $opt{output};
my $p_repdir = $p_outdir."/".$tools->{system}->{REPDIR};
my $fileinfo_out  = "$p_repdir/".$tools->{system}->{FILEINFO_OUT};
my $res_usage_out = "$p_repdir/".$tools->{system}->{RESUSAGE_OUT};
my $summary_out   = "$p_repdir/".$tools->{system}->{SUMMARY_OUT};

# retrieve input files
my $file_info;
my $cmd;

$file_info = &restore_filelist($opt{list});

#create excel file for each samples
my $count=1;
foreach my $idx ( sort {$a<=>$b} keys %$file_info )
{
	my $fnb = $file_info->{$idx}->{PREFIX};
	my $outfile = "$p_repdir/report_SEQ$idx\_$fnb.xlsx";
	print STDERR "Generating $outfile...\n";
	my $workbook = Excel::Writer::XLSX->new( $outfile );

	foreach my $tool ( sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools )
	{
		next if $tool eq 'system';
		next unless -e "$p_repdir/$idx\_$fnb/$tool/$fnb-$tool.list.txt";	
		my $worksheet = $workbook->add_worksheet( $tool );
		$worksheet->set_default_row( 15 );

		my $fmt_h = $workbook->add_format( align => 'center', bg_color => '#D9D9D9', size => 12 );
		my $fmt_b = $workbook->add_format( size => 12 );
	
		my @data = &csv2array( "$p_repdir/$idx\_$fnb/$tool/$fnb-$tool.list.txt" );
		&worksheet_set_column_wid( $worksheet, @data );

		my @header = shift @data;
		$worksheet->write_col( 'A1', \@header, $fmt_h );
		$worksheet->write_col( 'A2', \@data, $fmt_b );
		&worksheet_conditional_formatting( $workbook, $worksheet, "B:B", $opt{highlight_list}, "genus" );

		$worksheet->freeze_panes( 1, 0 );
	}
	$count++;
	$workbook->close();
}

#create excel file for each tools
foreach my $tool ( sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools )
{
	next if $tool eq 'system';
	my $order=$tools->{$tool}->{ORDER};
	my $outfile = "$p_repdir/report_TOOL$order\_$tool.xlsx";
	print STDERR "Generating $outfile...\n";

	$count=1;
	my $workbook = Excel::Writer::XLSX->new( $outfile );
	foreach my $idx ( sort {$a<=>$b} keys %$file_info )
	{
		my $fnb = $file_info->{$idx}->{PREFIX};
		my ($sheetname) = $fnb =~ /^(.{1,30})/;
		my $worksheet = $workbook->add_worksheet( $sheetname );
		$worksheet->set_default_row( 15 );

		my $fmt_h = $workbook->add_format( align => 'center', bg_color => '#D9D9D9', size => 12 );
		my $fmt_b = $workbook->add_format( size => 12 );
		
		next unless -e "$p_repdir/$idx\_$fnb/$tool/$fnb-$tool.list.txt";
		my @data = &csv2array( "$p_repdir/$idx\_$fnb/$tool/$fnb-$tool.list.txt" );
		&worksheet_set_column_wid( $worksheet, @data );

		my @header = shift @data;
		$worksheet->write_col( 'A1', \@header, $fmt_h );
		$worksheet->write_col( 'A2', \@data, $fmt_b );
		&worksheet_conditional_formatting( $workbook, $worksheet, "B:B", $opt{highlight_list}, "genus" );

		# freeze header
		$worksheet->freeze_panes( 1, 0 );
	
		$count++;
	}
	$workbook->close();
	unlink $outfile if ! -e "$p_repdir/1_allReads//$tool/allReads-$tool.list.txt";
}

$count=1;
print STDERR "Generating $p_repdir/report_summary.xlsx...\n";
my $workbook = Excel::Writer::XLSX->new( "$p_repdir/report_summary.xlsx" );

my @files = (
	"$fileinfo_out",
	"$summary_out",
	"$res_usage_out"
);

foreach my $file ( @files ){
	my ($name) = $file =~ /([^\/]+)\.[^\.]+$/;
	my $worksheet = $workbook->add_worksheet( $name );

	my $fmt_h = $workbook->add_format( align => 'center', bg_color => '#D9D9D9', size => 12 );
	my $fmt_b = $workbook->add_format( size => 12 );

	my @data = &csv2array( $file );
	&worksheet_set_column_wid( $worksheet, @data );

	my @header = shift @data;
	$worksheet->write_col( 'A1', \@header, $fmt_h );
	$worksheet->write_col( 'A2', \@data, $fmt_b );

	if( $name =~ /summary/ ){
		$worksheet->set_column( "D:H", 31);
		&worksheet_conditional_formatting( $workbook, $worksheet, "D:H" );

		# merge dataset cells
		my $tolrow = scalar @data;
		my $fmt = $workbook->add_format( valign => 'vcenter', align  => 'left' );
		# merging tool cells
		for( my $i=2; $i<$tolrow+2; $i+=3 ){
			$worksheet->merge_range( "B$i:B".($i+2), $data[$i-2][1], $fmt );
		}
		# merging dataset cells
		$fmt = $workbook->add_format( valign => 'vcenter', align  => 'left' );
		my ($dataset,$start) = ($data[0][0],2);
		for( my $i=2; $i<$tolrow+3; $i++ ){
			my $temp = defined $data[$i-2][0] ? $data[$i-2][0] : "";
			if( $dataset ne $temp ){
				$worksheet->merge_range( "A$start:A".($i-1), $dataset, $fmt );
				$dataset = $temp;
				$start = $i;
			}
		}
	}
	
	$worksheet->freeze_panes( 1, 0 );
}

$workbook->close();

###############################################################

sub restore_settings {
	my ( $file ) = @_;
	open FILE, $file or die "Can't open settings $file: $!\n";
	my $set;
	my $section;
	my $count=0;
	while(<FILE>){
		chomp;
		next if /^$/;
		next if /^#/;
		next if /^;;/;
		next if (! /\S/);
		if ( /^\[(.+)\]$/ ){
			$section = $1;
			$set->{$section}->{ORDER} = $count++;
			next;
		}
		my ($key, $val) = $_ =~ /^(.*)=(.*)$/;
		$set->{$section}->{$key} = $val;
	}
	close FILE;
	return $set;
}

sub restore_filelist {
	my $file = shift;
	open LIST, $file or die "Can't open $file: $!\n";
	#<FASTQ> <FASTQp1,FASTQp2> <NAME> <FASTA> <FASTA_EXTRACT>
	my $count=1;
	my $file_info;
	my @filelist_header;
	while(<LIST>){
		chomp;
		next if /^--/;
		next if /^\s*$/;
		next if /^#/;
			
		if (/^PREFIX/){
		    @filelist_header = split /\t/, $_; 
		    next;
		}   
		
		my @fields = split /\t/, $_; 
		for (my $i=0; $i<=$#filelist_header; $i++){
		    $file_info->{$count}->{$filelist_header[$i]} = $fields[$i];
		}   
		$count++;
	}
	close LIST;
	return $file_info;
}

sub csv2array {
    my $file = shift;
    my @all_data;
    open FILE, $file or die "Can't open file $file: $!\n";
    while (<FILE>) {
        chomp;
		next if /^$/;
		next if /^-+$/;
		my $line = $_;
		$line =~ s/ {2,}/\t/g;
        push @all_data, [ split /\t/, $line ];
    }
    close FILE;
    return @all_data;
}

sub worksheet_set_column_wid {
	my ($worksheet, @data) = @_;
	
	my @width;
	for( my $i=0; defined $data[$i]; $i++ ){
		for( my $j=0; defined $data[$i][$j]; $j++ ){
			$width[$j] = length $data[$i][$j] unless defined $width[$j];
			if( length $data[$i][$j] > $width[$j] ){
				$width[$j] = length $data[$i][$j];
			}
		}
	}

	for( my $i=0; defined $width[$i]; $i++ ){
		#my $index = &convertToAlpha($i+1);
		$worksheet->set_column( $i, $i, $width[$i]+2 > 70 ? 70 : $width[$i]+2 );
	}
}

sub worksheet_conditional_formatting {
	my ($workbook, $worksheet, $range, $file, $lvl) = @_;

	my $fmt_na = $workbook->add_format( color => '#D9D9D9', size => 12 );	
	my $fmt_red = $workbook->add_format( color => '#9C0006', size => 12 );

	$worksheet->conditional_formatting( $range,{
			type     => 'text',
			criteria => 'begins with',
			value    => 'N/A',
			format   => $fmt_na,
		}
	);

	if( defined $file && -e $file ){
		open LIST, $file or die "FAILED: Can't open $file: $!\n";
		while(<LIST>){
			chomp;
			my ($g,$s) = $_ =~ /^(\w+) (\w+)/;
			($g) = $_ =~ /^(\w+)/ unless defined $g;
			my $patten = $g;
			$patten = "$g $s" if ( $lvl eq 'species' || $lvl eq 'strain' || $lvl eq 'replicon' );
			die "Unknow patten" if $patten eq " ";
	
			$worksheet->conditional_formatting( $range,{
					type     => 'text',
					criteria => 'begins with',
					value    => $patten,
					format   => $fmt_red,
				}
			);
		}
		close LIST;
	}

	#$worksheet->conditional_formatting( $range,{
	#		type     => 'text',
	#		criteria => 'begins with',
	#		value    => 'Bacillus',
	#		format   => $fmt_red,
	#	}
	#);
	#$worksheet->conditional_formatting( $range,{
	#		type     => 'text',
	#		criteria => 'begins with',
	#		value    => 'Yersinia',
	#		format   => $fmt_red,
	#	}
	#);
	#$worksheet->conditional_formatting( $range,{
	#		type     => 'text',
	#		criteria => 'begins with',
	#		value    => 'Francisella',
	#		format   => $fmt_red,
	#	}
	#);
}

