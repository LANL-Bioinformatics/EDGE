#!/usr/bin/perl

use strict;
use utf8;
use Encode qw(encode_utf8);
use JSON;
use LWP::UserAgent;
use HTTP::Request::Common;
use FindBin qw($RealBin);
use lib "$RealBin/edge_ui/metadata_scripts/lib";
use lib "$RealBin/../metadata_scripts/lib";
use SampleMetadata;
use Pathogen;
use EdgeSite;
use Travel;
use Symptom;

sub pushEdgeSite {
	my $metadata = shift;
	my $obj = new EdgeSite($metadata->{'edgesite-organization'}, $metadata->{'edgesite-acronym'}, $metadata->{'edgesite-location'},$metadata->{'edgesite-city'}, $metadata->{'edgesite-state'}, $metadata->{'edgesite-country'}, $metadata->{'edgesite-lat'}, $metadata->{'edgesite-lng'});
	my $data = $obj->toJson();

	#Set the request parameters
	my $url = $metadata->{'bsve_api_url'}."/user/register";
	my $wsinfo =  apiWS_service($url, $data,"POST");

	if($wsinfo->{error}) {
		return ("error", $wsinfo->{error});
	} else {
		my $apiKey =  $wsinfo->{"api_key"};
		my $apiToken =  $wsinfo->{"api_token"};
		return ($apiKey, $apiToken);
	}
}

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

	my $metadataFile = "$proj_dir/metadata_sample.txt";
	my $runFile = "$proj_dir/metadata_run.txt";
	my $travelFile = "$proj_dir/metadata_travels.txt";
	my $symptomFile = "$proj_dir/metadata_symptoms.txt";
	my $pathogensFile = "$proj_dir/pathogens.txt";
	my $otherFile = "$proj_dir/metadata_other.txt";
	if(hasPathogens($pathogensFile)) {
		my $metadata = &getMetadataParams($metadataFile);
		my $run = &getMetadataParams($runFile);
		my $other = &getMetadataParams($otherFile);

			#for multiple projects share-metadata-with-bsve action
			if($run->{'bsve_id'}) {
				$action = "update";
			} else {
				$action = "add";
			}
			#end

			my $study_id=$metadata->{'study_id'};
			if($study_id =~ /^\d+$/) {
				$study_id = 'S'.$sys->{'bsve_api_key'}.$metadata->{'study_id'};
			}
	
			#use site location for sample location 
			if($metadata->{'study_type'} ne "SRA") {
				$metadata->{'location'} = $sys->{'edgesite-location'} unless $metadata->{'location'};
				$metadata->{'city'} = $sys->{'edgesite-city'} unless $metadata->{'city'};
				$metadata->{'state'} = $sys->{'edgesite-state'} unless $metadata->{'state'};
				$metadata->{'country'} = $sys->{'edgesite-country'} unless $metadata->{'country'};
				$metadata->{'lat'} = $sys->{'edgesite-lat'} unless $metadata->{'lat'} ;
				$metadata->{'lng'} = $sys->{'edgesite-lng'} unless $metadata->{'lng'};
			}

			my $obj = new SampleMetadata($run->{'bsve_id'},$metadata->{'sra_run_accession'},$study_id,$metadata->{'study_title'},$metadata->{'study_type'},$metadata->{'sample_name'},$metadata->{'sample_type'},$metadata->{'host'}, $metadata->{'host_condition'}, $metadata->{'gender'}, $metadata->{'age'}, $metadata->{'isolation_source'}, $metadata->{'collection_date'}, $metadata->{'location'},$metadata->{'city'}, $metadata->{'state'}, $metadata->{'country'}, $metadata->{'lat'}, $metadata->{'lng'}, $metadata->{'experiment_title'}, $metadata->{'sequencing_center'}, $metadata->{'sequencer'}, $metadata->{'sequencing_date'});

			my $run_id = $run->{'edge-run-id'};
			$run_id = $configuration->{'projid'} unless $run_id;
			my $run_host = $configuration->{'projrunhost'};
			if($run_id =~ /^\d+$/) {
				$run_id = 'R'.$sys->{'bsve_api_key'}.$run_id;
			} else {
				#SRA run
				#$run_host = "SRA";
			}
			my $data = $obj->toJson($sys->{'bsve_api_key'}, $sys->{'bsve_api_token'}, $run_host, $run_id, $other->{'lanl_only'});

			print MLG "$action sample metadata\n";
			print MLG $data."\n" if($sys->{'bsve_api_debug'}) ;

			#Set the request parameters
			my $url = $sys->{'bsve_api_url'}."/metadata/add";
			print MLG "$url\n" if($sys->{'bsve_api_debug'}) ;
			my $wsinfo =  apiWS_service($url, $data,"POST");

			if($wsinfo->{error}) {
				print MLG "Failed to $action sample metadata to BSVE API server: ".$wsinfo->{error}."\n";
				$success = 0;
			} else {
				my $id =  $wsinfo->{"id"};

				if($action eq "add") {
					#append id to sample_metadata.txt
					open META, ">>$runFile";
					print META "bsve_id=$id\n";
					close META;
				}

				print MLG "bsve_id = $id\n";

				#push travels to  bsve WS
				if(-e $travelFile) {
					my $turl = $sys->{'bsve_api_url'}."/metadata/travel/add";
					print MLG "$turl\n" if($sys->{'bsve_api_debug'}) ;

					#parse file
					open TVL, "$travelFile";
					my ($from, $to, $location, $city, $state, $country, $lat, $lng);
				 	while(<TVL>){
				      		chomp;
						next if(/^#/);
				     		if ( /(.*)=(.*)/ ){
							if ($1 eq "travel-date-from") {
								$from = $2;
							} elsif ($1 eq "travel-date-to") {
								$to = $2;
							} elsif ($1 eq "travel-location") {
								$location = $2;
							} elsif ($1 eq "city") {
								$city = $2;
							} elsif ($1 eq "state") {
								$state = $2;
							} elsif ($1 eq "country") {
								$country = $2;
							} elsif ($1 eq "lat") {
								$lat = $2;
							} elsif ($1 eq "lng") {
								$lng = $2;
								my $travel = new Travel($id, "$from ~ $to", $location, $city, $state, $country, $lat, $lng);
								my $tdata = $travel->toJson($sys->{'bsve_api_key'}, $sys->{'bsve_api_token'}, $other->{'lanl_only'});
								print MLG "Add travel: ".$tdata."\n" if($sys->{'bsve_api_debug'}) ;

								#submit to WS
								my $wsinfo = apiWS_service($turl, $tdata,"POST");
								if($wsinfo->{error}) {
									print MLG "Failed to add travel to WS: ".$wsinfo->{error}."\n";
									$success = 0;
						              	}
							}
						}
					}
					close TVL;
				} else {
					print MLG "File $travelFile not found.\n";
				}

				#push symptoms to  bsve WS
				if(-e $symptomFile) {
					my $surl = $sys->{'bsve_api_url'}."/metadata/symptom/add";
					print MLG "$surl\n" if($sys->{'bsve_api_debug'}) ;

					#parse file
					open SM, "$symptomFile";
					while(<SM>) {
						chomp;
						next if(/^#/);
				     		if ( /(.*)\t(.*)/ ){
							my $symptom = new Symptom($id, $1, $2);
							my $sdata = $symptom->toJson($sys->{'bsve_api_key'}, $sys->{'bsve_api_token'}, $other->{'lanl_only'});
							print MLG "Add symptom: ".$sdata."\n" if($sys->{'bsve_api_debug'}) ;

							#submit to WS
							my $wsinfo = apiWS_service($surl, $sdata,"POST");
							if($wsinfo->{error}) {
								print MLG "Failed to add symptom to WS: ".$wsinfo->{error}."\n";
								$success = 0;
				                      	}
						}
					}
					close SM;
				} else {
					print MLG "File $symptomFile not found.\n";
				}

				#push pathogens to  bsve WS
				if(-e $pathogensFile) {
					my $purl = $sys->{'bsve_api_url'}."/metadata/pathogen/add";
					print MLG "$purl\n" if($sys->{'bsve_api_debug'}) ;

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
						my $pathogen = new Pathogen($id, $parts[0], $parts[1], $parts[2], $top);
						my $pdata = $pathogen->toJson($sys->{'bsve_api_key'}, $sys->{'bsve_api_token'}, $other->{'lanl_only'});
						print MLG "Add pathogen: ".$pdata."\n" if($sys->{'bsve_api_debug'}) ;

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
	} else {
		print MLG "File $metadataFile not found.\n";
	}

	close MLG;
	return $success;
}

sub hasPathogens {
	my $file = shift;
	my $top = 0;

	if(-e $file) {
	    	open (my $fh , $file) or die "No config file $!\n";
	   	 while (<$fh>) {
	       	 	chomp;
			if(/pathogen\thost\(s\)\tdisease/) {
				next;
			}
			$top ++;
	    	}
	    	close $fh;
	}

 	return $top;
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

	if(-e $config) {
		open CONF, $config or die "Can't open $config: $!";
		while(<CONF>){
	      		chomp;
		        next if(/^#/);
		   	if ( /(.*)=(.*)/ ){
		     		$sys->{$1}=$2;
		      	}
		}
		close CONF;
	}
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
	$data = encode_utf8($data);
	$req->header('Content-Type' => 'application/json;charset=UTF-8');
	$req->header('Accept' => 'application/json;charset=UTF-8');
	#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
	$req->header( "Content-Length" => length($data) );
	$req->content($data);

	my $response = $browser->request($req);
	my $result_json = $response->decoded_content;
#print STDERR "$result_json\n";
	my $wsinfo =  from_json($result_json);
	if ($result_json =~ /\"error_msg\":"(.*)"/) {
		$wsinfo->{error}=$1;
   	}

	return $wsinfo;
}

1;
