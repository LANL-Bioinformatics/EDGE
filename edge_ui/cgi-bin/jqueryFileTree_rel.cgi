#!/usr/bin/env perl
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use HTML::Entities ();

#-----------------------------------------------------------
#  jQuery File Tree Perl Connector
#
#  Version 1.0
#
#  Oleg Burlaca
#  http://www.burlaca.com/2009/02/jquery-file-tree-connector/
#  12 February 2009
#-----------------------------------------------------------

# for security reasons,  specify a root folder 
# to prevent the whole filesystem to be shown
# for ex: the root folder of your webbrowser

# read system params from sys.properties
my $config_tmpl = "$RealBin/../sys.properties";
my $sys         = &getSysParamFromConfig($config_tmpl);
my $root        = $sys->{edgeui_wwwroot};

#----------------------------------------------------------

my $params = &getCGIParams();
print "Content-type: text/html\n\n";

my $dir = $params->{dir};
if ($dir =~ /[^0-9a-zA-Z_\/\-\.]/){print "Error\n"; exit;}
if ($dir =~ /\.\.\//){print "Error\n"; exit;}
if ($dir !~ /EDGE_output\/+\w+/){print "Error\n"; exit;}
my $fullDir = $root . $dir;

exit if ! -e $fullDir;

opendir(BIN, $fullDir) or die "Can't open $dir: $!";
my (@folders, @files);
my $total = 0;
while( defined (my $file = readdir BIN) ) {
	next if $file =~ /^\./;
    $total++;
    if (-d "$fullDir/$file") {
	push (@folders, $file);
    } else {
	push (@files, $file);
    }
}
closedir(BIN);

return if $total == 0;
print "<ul class=\"jqueryFileTree\" style=\"display: none;\">";

# print Folders
foreach my $file (sort @folders) {
    next if ! -e  $fullDir . $file;
    
    print '<li class="directory collapsed"><a href="#" rel="' . 
          &HTML::Entities::encode($dir . $file) .'/">' . 
          &HTML::Entities::encode($file) . '</a></li>';
}

# print Files
foreach my $file (sort @files) {
    next if ! -e  $fullDir . $file;

    $file =~ /\.([^\.]+)(\.gz)*$/;
    my $ext = $1;
    print '<li class="file ext_' . $ext . '"><a href="#" rel="' . 
    &HTML::Entities::encode( $dir . $file) . '">' .
    &HTML::Entities::encode($file) . '</a></li>';
}

print "</ul>\n";




#--------------------------------------------------------------------------------------------------
sub getCGIParams {
    my $line;
    
    if ($ENV{'REQUEST_METHOD'} eq "POST") {
        read(STDIN, $line, $ENV{'CONTENT_LENGTH'});
    } else {
        $line = $ENV{'QUERY_STRING'};
    }

    my (@pairs) = split(/&/, $line);
    my ($name, $value, %F);
        
    foreach (@pairs) {
        ($name, $value) = split(/=/);
        $value =~ tr/+/ /;
        $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
        
        if (! exists $F{$name}) {
            $F{$name} = $value;
        } elsif (exists $F{$name} and ref($F{$name}) ne 'ARRAY') {
            my $prev_value = $F{$name};
            delete $F{$name};
            $F{$name} = [ $prev_value, $value ];
	} else { push @{ $F{$name} }, $value }
    }
    return \%F;
}

sub getSysParamFromConfig {
	my $config = shift;
	my $sys;
	my $flag=0;
	open CONF, $config or die "Can't open $config: $!";
	while(<CONF>){
		if( /^\[system\]/ ){
			$flag=1;
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
	die "Incorrect system file\n" if (!$flag);
	return $sys;
}

#--------------------------------------------------------------------------------------------------
