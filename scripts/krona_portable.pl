#!/usr/bin/env perl
# krana_portable.pl
# ver 0.1
# 2014/02/05
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

use strict;
use Getopt::Long;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$ENV{PATH}";
my $ktpath;

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'inhtml|i=s',
    'outhtml|o=s',
    'jspath|p=s',
    'help|?') || &usage();

if ( $opt{help} || !-e $opt{inhtml} || !defined $opt{outhtml} ) { &usage(); }

# find full path of krona-2.0.js
$ktpath = $opt{jspath};
unless( defined $ktpath && -e $ktpath ){
    $ktpath ||= `ktGetLibPath`;
    die "Can't find Krona path of krona-2.0.js" unless $ktpath;
    $ktpath .= "/../src/krona-2.0.js";
}

# read krona-2.0.js
my $js="";
open JS, $ktpath or die "Can't open $ktpath\n";
while(<JS>){
    my $line = $_;
    if( /if \( image\.complete \)/ ){
	while(<JS>){
	    $line = $_;
	    last if /^}/;
	}
    }
    $js .= $line;
}
close JS;

# replace and output
open INHTML, $opt{inhtml} or die "Can't open $opt{inhtml}.\n";
open OUTHTML, ">$opt{outhtml}" or die "Can't open $opt{outhtml}.\n";
foreach(<INHTML>){
    if( /<script src="http:\/\/krona\.sourceforge\.net\/src\/krona-2\.0\.js"><\/script>/){
        print OUTHTML "<script type=\"text/javascript\">\n$js\n</script>\n";
    }
    else{
        print OUTHTML $_;
    }
}
close INHTML;
close OUTHTML;

sub usage {
print <<__END__;
Convert internet-required KRONA html to green and portable version. Note that images still display
improperly without internet connection, but all functions are working.

$0 [OPTIONS] --inhtml <KRONA_HTML> --outhtml <HTML>

    --inhtml  | -i <STRING>  internet-required Krona html
    --outhtml | -o <STRING>  output green Krona html

[OPTIONS]

    --jspath  | -p <STRING>  full path of krona-2.0.js

__END__
exit();
}
