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
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/../..";
$ENV{PATH} = "$EDGE_HOME/bin:$ENV{PATH}";

my $ref_json_file = "$edgeui_wwwroot/data/Ref_list.json";

my @names = grep {  length($_)>=3 } split / /,$queries;
my $search_str = join " or ", map { "contains(".'"'."$_".'")' } @names;
# case incenstive search is 2x slower. 
#my $search_str = join " or ", map { "test(".'"'."$_".'";"i")' } @names;
my $jq_cmd = "jq -c '.| with_entries( select(.key|$search_str))' $ref_json_file";
#print STDERR $jq_cmd,"\n";

my $return_json = `$jq_cmd`;

print $cgi->header('application/json');
print $return_json;
exit;
#&returnStatus();

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
