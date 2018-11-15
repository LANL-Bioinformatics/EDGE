#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opt;
my $res=GetOptions(\%opt,
                   'list=s',
                   'level|l=s',
                   'setting|s=s',
                   'top|t=i',
                   'output|o=s',
                   'help|h|?') || &usage();

$opt{level} = "genus,species,strain" unless defined $opt{level};
$opt{top} = 5 unless ( defined $opt{top} && $opt{top} );

my %levels_h = map { $_ => 1 } split /\s*,\s*/, $opt{level};
my @levels_a = split /\s*,\s*/, $opt{level};

my $result;

#read settings
my $ini_file = $opt{setting};
my $tools = &restore_settings( $ini_file );

my $file_info;
$file_info = &restore_filelist($opt{list});

#default settings
my $p_outdir = $tools->{system}->{OUTDIR};
$p_outdir    = $opt{output} if defined $opt{output};

my $p_repdir = $p_outdir."/".$tools->{system}->{REPDIR};
my $count=1;

#submit tools
foreach my $idx ( sort {$a<=>$b} keys %$file_info )
{
	my $fnb = $file_info->{$idx}->{PREFIX};

	foreach my $tool ( sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools )
	{
		next if $tool eq 'system';
		my $file = "$p_repdir/$count\_$fnb/$tool/$fnb-$tool.list.txt";
		print STDERR "Processing $file...\n";

		open LIST, $file || die "Can't open $file: $!\n";
		while(<LIST>){
			chomp;
			my @fields = split /\t/, $_;
			my $lvl = $fields[0];
			next unless $levels_h{$lvl};
			next unless $fields[2];
			$result->{$fnb}->{$tool}->{$lvl}->{$fields[1]} = $fields[2];
		}
		close LIST;
	}
	$count++;
}

print STDERR "Creating TOP$opt{top} list...\n";

#print header
print "DATASET\tTOOL\tLEVEL";
for(my $i=1; $i<=$opt{top}; $i++){
	print "\tTOP$i";
}
print "\n";

#print content
foreach my $idx ( sort {$a<=>$b} keys %$file_info )
{
	my $fnb = $file_info->{$idx}->{PREFIX};
	foreach my $tool ( sort {$tools->{$a}->{ORDER}<=>$tools->{$b}->{ORDER}} keys %$tools )
	{	
		next if $tool eq 'system';
		foreach my $lvl ( @levels_a )
		{
			print "$fnb\t$tool\t$lvl";
			
			#no results at certain level
			unless ( defined $result->{$fnb}->{$tool}->{$lvl} ){
				print "\tN/A"x$opt{top},"\n";
				next;
			}
	
			#print top # results
			my $count=0;
			foreach my $taxa ( sort {$result->{$fnb}->{$tool}->{$lvl}->{$b} <=> $result->{$fnb}->{$tool}->{$lvl}->{$a} or $a cmp $b } keys %{$result->{$fnb}->{$tool}->{$lvl}} ){
				print "\t$taxa";
				last if $count++ == $opt{top}-1;
			}

			#fill "N/A" if the results do not fulfill top number
			print "\tN/A"x($opt{top}-$count);
			print "\n";
		}
	}
}
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
		next if /^;/;
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
