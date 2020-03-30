#!/usr/bin/env perl
# For import project directories into EDGE development version  with a user mangement system. 
use strict;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use LWP::UserAgent;
use HTTP::Request::Common;
use Data::Dumper;
use Digest::MD5 qw(md5_hex);
use JSON;

my $old_project_dir = $ARGV[0];
my $username = $ARGV[1];
my $password = $ARGV[2];
my $configFile = "$old_project_dir/config.txt";
my $EDGE_HOME = $ENV{EDGE_HOME};
$EDGE_HOME ||= "$RealBin/..";

if (!@ARGV){print "$0 import_project_dir Username Password\n";  exit;}

my $configuration=&readConfig($configFile);
my $sys = &getSysParams("$EDGE_HOME/edge_ui/sys.properties");

my $owner=$configuration->{projowner};
my $projname = $configuration->{projname};
my $projdesc = $configuration->{projdesc};
my $projcode = $configuration->{projcode};
my $projid = $configuration->{projid};
my $edge_output= $sys->{edgeui_output};
my $edge_input=$sys->{edgeui_input};
my $user_myproject_dir="$edge_input/". md5_hex(lc($username))."/MyProjects";
my $um_url      = $sys->{edge_user_management_url};

print join("\n","Old/Import project", "project_owner: $owner","project_name: $projname","project_description: $projdesc","\n");
print join("\n","New system", "edge_output_dir: $edge_output","um_system: $um_url","\n");

if (!$username && !$password){ print "No username/password\n"; exit;}
## Add project to database and assign/update new id/code.
my %data = (
	email => $username,
	password => $password,
	project_name => $projname,
	description => $projdesc,
);

my $data = to_json(\%data);
my $url = $um_url ."WS/project/add";
my $browser = LWP::UserAgent->new;
my $req = POST $url;
$req->header('Content-Type' => 'application/json');
$req->header('Accept' => 'application/json');
$req->header( "Content-Length" => length($data) );
$req->content($data);
my $response = $browser->request($req);
my $result_json = $response->decoded_content;
#print $result_json,"\n";
my $info =  from_json($result_json);
my $new_projid =  $info->{"id"};
my $new_projcode = $info->{"code"};

# make sure we got new id/code from rest api(databasae) before moving project
if ($new_projid && $new_projcode){
	print "Importing $projname\n";
	# update JBrowse symlink
	system("rm -f $EDGE_HOME/edge_ui/JBrowse/data/$new_projcode");
	system("ln -s $edge_output/$new_projcode $EDGE_HOME/edge_ui/JBrowse/data/$new_projcode");
        system("ln -s $edge_output/$new_projcode $user_myproject_dir/${projname}_$new_projid");
	# move data
	system("mv -f $old_project_dir $edge_output/$new_projcode");
	# update owner projid projcode;
	my $new_config = "$edge_output/$new_projcode/config.txt";
	my $new_processLog =  "$edge_output/$new_projcode/process.log";
	unlink "$edge_output/$new_projcode/HTML_Report/.complete_report_web";
	system("sed -i.bak 's/projowner=[[:graph:]]*/projowner=$username/; s/projid=[[:graph:]]*/projid=$new_projid/; s/projcode=[[:graph:]]*/projcode=$new_projcode/' $new_config");
	system("sed -i.bak 's/projowner=[[:graph:]]*/projowner=$username/g; s/projid=[[:graph:]]*/projid=$new_projid/g; s/projcode=[[:graph:]]*/projcode=$new_projcode/g; s/EDGE_output\\\/$projid/EDGE_output\\\/$new_projid/g' $new_processLog");
	system("sed -i.bak 's/$projcode/$new_projcode/g' $new_processLog") if ($projcode);
	my $pangia_vis_dir = "$edge_output/$new_projcode/ReadsBasedAnalysis/Taxonomy/report/1_allReads/pangia/pangia-vis";
	if ( -d "$pangia_vis_dir/$projcode"){
		system("mv $pangia_vis_dir/$projcode $pangia_vis_dir/$new_projcode");
		system("mv $pangia_vis_dir/$projcode.pangia.log $pangia_vis_dir/$new_projcode.pangia.log");
		system("mv $pangia_vis_dir/$projcode.tsv $pangia_vis_dir/$new_projcode.tsv");
	}
	print "Done importing \"$projname\". Please check EDGE to confirm it worked.\n";
}


######################## SUBS ########################
sub getSysParams {
        my $config = shift;
        my $sys;
        open CONF, $config or die "Can't open $config: $!";
        while(<CONF>){
                if( /^\[system\]/ ){
                        while(<CONF>){
                                chomp;
                                last if /^\[/;
                                next if(/^#/);
                                if ( /^([^=]+)=(.*)/ ){
                                        $sys->{$1}=$2;
                                }
                        }
                }
                last;
        }
        close CONF;
        return $sys;
}

sub readConfig
{
    my $file=shift;
    my %hash;
    open (my $fh , $file) or die "No config file $!\n";
    while (<$fh>)
    {
        chomp;
        next if (/^#/);
        if (/=/)
        {
            my ($key,$value)=split /=/,$_;
            if ( defined $value)
            {
               $value =~ s/\"//g;
               if ($key eq "Host")
               {
                 foreach my $each_host(split(/,/,$value))
                 {
                   push @{$hash{$key}} , $each_host;
                 }
               }
               elsif($key eq "reference")
               {
                 foreach my $each_ref(split(/,/,$value))
                 {
                   push @{$hash{$key}} , $each_ref;
                 }
               }
               else
               {
                   $hash{$key}=$value;
               }
            }
            else
            {
               $hash{$key}="";
            }
        }
    }
    close $fh;
    return \%hash;
}


