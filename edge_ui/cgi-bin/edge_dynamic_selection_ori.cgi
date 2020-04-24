#!/usr/bin/env perl

use strict;
use FindBin qw($RealBin);
use CGI qw(:standard);
use lib "$RealBin/../../lib";
use JSON;

my $info; # info to return
my $cgi     = CGI->new;
my %opt     = $cgi->Vars();
my $queries = $opt{'query'} || $ARGV[0];
&stringSanitization(\%opt);
# read system params from sys.properties
my $sysconfig    = "$RealBin/../sys.properties";
my $sys          = &getSysParamFromConfig($sysconfig);
my $edgeui_wwwroot = $sys->{edgeui_wwwroot};

my $ref_json_file = "$edgeui_wwwroot/data/Ref_list.json";

my $list = &readListFromJson($ref_json_file);
my @ref_list = keys %$list;


my @names = grep {  length($_)>=3 } split / /,$queries;
my $grep_term = join("|",@names);
if (@names){
	map { $info->{$_}=$list->{$_};} grep { /$grep_term/i } @ref_list;
}

&returnStatus();

sub returnStatus {
        my $json;
        $json = to_json( $info ) if $info;
        $json = to_json( $info, { ascii => 1, pretty => 1 } ) if $info && $ARGV[0];
        print $cgi->header('application/json'), $json;
        exit;
}

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
sub stringSanitization{
	my $opt=shift;
	foreach my $key (keys %opt){
		my $str = $opt->{$key};
		if($str =~ /[^0-9a-zA-Z\,\-\_\^\@\=\:\\\.\/\+ ]/){
			$info->{INFO} = "Invalid characters detected \'$str\'.";
			&returnStatus();
		}
	}
}
