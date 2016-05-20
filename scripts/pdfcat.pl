#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PDF::API2;


#   Process command line arguments and populate corresponding variables

my %pages = ();
my @input = ();

GetOptions(
    'h|help'            =>  sub {usage()},
    'i|input=s'         =>  \@input,
    'f|footer=s'		=>	\( my $footer = '' ),
    'o|output=s'        =>  \( my $output = '' ),
    'p|page|pages=s'    =>  sub {

        if (scalar @input > 0) {

        #   If an input file name has previously been defined, associa+te the given page 
        #   ranges to be extracted with the last input file name supplied.

            my @files = split /,/, $input[-1];
            push @{ $pages{ $_ } }, $_[1] foreach @files;
        }
    }
);

exit 1 unless scalar @input > 0 and length $output > 0;


#   Split the input files specified on any comma characters present - This 
#   allows for multiple input files to be specified either by multiple--input 
#   arguments or by a single argument in a comma delimited fashion.

@input = map { split /,/ } @input;


#   Open the PDF file for output (via the PDF::API2 object constructor)

my $pdf = PDF::API2->new( -file => $output );
my $root = $pdf->outlines;
my $font = $pdf->corefont('Helvetica');

#   Step through each of the input files specified and extract the document 
#   pages with the options specified.

my $import_page = 0;

foreach my $file ( @input ) {

    my $input = PDF::API2->open( $file );


    #   Expand the page list and range definitions passed with the --page argument 
    #   associated with the given input file.  By default, all pages of the input
    #   file are included in the output.

    my @pages = ();
    if ( exists $pages{ $file } ) {

    @pages = map { split /,/ } @{ $pages{ $file } };
    @pages = map { /^(\d+)-(\d+)$/ ? $1 .. $2 : $_ } @pages;
    }
    else {

        @pages = 1 .. $input->pages;
    }


    #   Import the pages from the input file input the output PDF file being 
    #   constructed

    if (scalar @pages > 0) {


        #   Extract the filename of the input file without the file extension for 
        #   incorporation into the document outline.

        my ($name, undef, undef) = fileparse($file, '\.[^\.]*');

        my $outline = $root->outline;
        $outline->title( $name );


        #   Step through each of the pages to be imported, import the page and add an 
        #   entry to the document outline.

        my $document_page = 0;
        foreach (@pages) {

            ++$import_page;
            ++$document_page;

            my $page = $pdf->importpage($input, $_, $import_page);
            if ($footer){
               my $text = $page->text();
               my ($llx, $lly, $urx, $ury) = $page->get_mediabox;
               $text->font($font, 9);
               $text->translate($urx-70-length($footer), 5);
               $text->text("$footer");
            }
    
            my $bookmark = $outline->outline;
            $bookmark->title("Page $document_page");
            $bookmark->dest($page);
            $outline->dest($page) if $document_page == 1;
        }
    }
}

$pdf->preferences( -outlines => 1 );
$pdf->update;
$pdf->end;

sub usage {

print <<USAGE;
pdfcat.pl [input files ...] [options] [output file]

    -i|--input [filename]

        Specify an input file for concatenation into the output file. If a single file is specified with the --page parameter, this script can also be used for extracting specific page ranges.

    -o|--output [filename]

        Specify the output file for concatenated PDF output.

    -p|--page|--pages

        This argument, which follows an input file argument, defines the pages to be extracted for concatenation from a given input file. If this argument is not defined, all pages from the input file are concatenated. The pages specified for extraction may be separated by commas or designed by ranges.
        
    -f|--footer
    
    	Add footer text in the right bottom for every pdf page.

USAGE
exit;
} 

exit 0;

