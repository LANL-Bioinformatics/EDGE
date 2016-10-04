#!/usr/bin/perl

use strict;
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use FindBin qw($RealBin);
use lib "$RealBin/edge_ui/metadata_scripts/lib";
use lib "$RealBin/../metadata_scripts/lib";
use SampleMetadata;
use SamplePathogen;

sub pushSampleMetadata {
	my $action = shift;
	my $proj_dir = shift;
	my $sys = shift;

	## Global hash reference for configurations.
	my $configFile = "$proj_dir/config.txt";
	my $configuration=&readConfig($configFile);
	my $logFile = "$proj_dir/metadata_bsve.log";

	my $success = 1;
	open MLG, ">>$logFile";
	my $timestamp = localtime(time);
	print MLG $timestamp."\n";

	my $metadataFile = "$proj_dir/sample_metadata.txt";
	my $pathogensFile = "$proj_dir/pathogens.txt";
	if(-e $metadataFile) {
		my $metadata = &getMetadataParams($metadataFile);

		if($action eq "delete") {
			my $obj = new SampleMetadata($metadata->{'bsve_id'});
			my $data = $obj->toJson($sys->{'sample_metadata_api_key'}, $sys->{'sample_metadata_api_token'});

			print MLG "Delete sample metadata\n";
			print MLG $data."\n" if($sys->{'sample_metadata_api_debug'});

			#Set the request parameters
			my $url = $sys->{'sample_metadata_api_url'}."/metadata/delete";
			print MLG "$url\n" if($sys->{'sample_metadata_api_debug'}) ;
			my $wsinfo =  apiWS_service($url, $data,"PUT");

			if($wsinfo->{error}) {
				print MLG "Failed to delete sample metadata to EDGE API server: ".$wsinfo->{error}."\n";
				$success = 0;
			} else {
				#remove bsve_id from sample_metadata.txt
				`perl -ni -e 'if(!/bsve_id=/){print;}' $metadataFile`;
			}
		} else {
			#for multiple projects share-metadata-with-bsve action
			if($metadata->{'bsve_id'}) {
				$action = "update";
			} else {
				$action = "add";
			}
			#end

			my $obj = new SampleMetadata($metadata->{'bsve_id'}, ,$metadata->{'study_title'},$metadata->{'sample_name'},$metadata->{'type'},$metadata->{'experiment_title'}, $metadata->{'host'}, $metadata->{'host_condition'}, $metadata->{'gender'}, $metadata->{'age'}, $metadata->{'source'}, $metadata->{'source_detail'}, $metadata->{'collection_date'}, $metadata->{'location'},$metadata->{'city'}, $metadata->{'state'}, $metadata->{'country'}, $metadata->{'lat'}, $metadata->{'lng'}, $metadata->{'seq_date'}, $metadata->{'seq_platform'}, $metadata->{'sequencer'},$metadata->{'instrument_model'},$metadata->{'center_name'});

			if($metadata->{'run_host'}) {
				$configuration->{'projrunhost'} = $metadata->{'run_host'};
			}
			if($metadata->{'run_id'}) {
				$configuration->{'projid'} = $metadata->{'run_id'};
			}
			my $data = $obj->toJson($sys->{'sample_metadata_api_key'}, $sys->{'sample_metadata_api_token'},  $configuration->{'projrunhost'}, $configuration->{'projid'});

			print MLG "$action sample metadata\n";
			print MLG $data."\n" if($sys->{'sample_metadata_api_debug'}) ;

			#Set the request parameters
			my $url = $sys->{'sample_metadata_api_url'}."/metadata/add";
			print MLG "$url\n" if($sys->{'sample_metadata_api_debug'}) ;
			my $wsinfo =  apiWS_service($url, $data,"POST");

			if($wsinfo->{error}) {
				print MLG "Failed to $action sample metadata to EDGE API server: ".$wsinfo->{error}."\n";
				$success = 0;
			} else {
				my $id =  $wsinfo->{"id"};

				if($action eq "add") {
					#append id to sample_metadata.txt
					open META, ">>$metadataFile";
					print META "bsve_id=$id\n";
					close META;
				}

				print MLG "bsve_id = $id\n";
				#push pathogens to  bsve WS
				if(-e $pathogensFile) {
					my $purl = $sys->{'sample_metadata_api_url'}."/metadata/pathogen/add";
					print MLG "$purl\n" if($sys->{'sample_metadata_api_debug'}) ;

					#parse file
					open PG, "$pathogensFile";
					my $top = 0;
					while(<PG>) {
						chomp;
						if(/pathogen\thost\(s\)\tdisease/) {
							next;
						}
						$top ++;
						my @parts = split(/\t/);
						my $pathogen = new SamplePathogen($id, $parts[0], $parts[1], $parts[2], $top);
						my $pdata = $pathogen->toJson($sys->{'sample_metadata_api_key'}, $sys->{'sample_metadata_api_token'});
						print MLG "Add pathogen: ".$pdata."\n" if($sys->{'sample_metadata_api_debug'}) ;

						#submit to WS
						my $wsinfo = apiWS_service($purl, $pdata,"POST");
						if($wsinfo->{error}) {
							print MLG "Failed to add pathogen to WS: ".$wsinfo->{error}."\n";
							$success = 0;
		                              	}
					}
					close PG;
				} else {
					print MLG "File $pathogensFile not found.\n";
				}
			}
		}
	} else {
		print MLG "File $metadataFile not found.\n";
	}

	close MLG;
	return $success;
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

sub getMetadataParams {
        my $config = shift;
        my $sys;
        open CONF, $config or die "Can't open $config: $!";
        while(<CONF>){
      		chomp;
                next if(/^#/);
           	if ( /(.*)=(.*)/ ){
             		$sys->{$1}=$2;
              	}
        }
        close CONF;
        return $sys;
}

sub apiWS_service {
	my $url = shift;
	my $data = shift;
	my $method = shift;
	my $browser = LWP::UserAgent->new;
	my $req;
	if ($method eq "POST") {
		$req = POST $url;
	} else {
		$req = PUT $url;
	}
	$req->header('Content-Type' => 'application/json');
	$req->header('Accept' => 'application/json');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;

	my $wsinfo =  from_json($result_json);
	if ($result_json =~ /\"error_msg\":"(.*)"/) {
		$wsinfo->{error}=$1;
   	}

	return $wsinfo;
}

1;
