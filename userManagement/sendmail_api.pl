#!/usr/bin/env perl

use Getopt::Long;
use strict;
use HTTP::Request::Common;
use LWP::UserAgent;
use JSON;
$ENV{'PERL_LWP_SSL_VERIFY_HOSTNAME'} = 0;

my $apiUrl = "http://edge_ncov_edge_mail_1:5009/mail/sendmail";
my $token = "McQ%NHzMxrxb+&Wj";

my ($from, $to, $cc, $bcc, $subject, $message, $text);

GetOptions(
           "from=s"           => \$from,
           "to=s"             => \$to,
           "cc=s"             => \$cc,
           "bcc=s"             => \$bcc,
           "subject=s"        => \$subject,
           "message=s"        => \$message,
           "text"             => \$text,
           "help|?"           => sub{Usage()} );
          
Usage() unless $from && $to && $subject && $message;

# call sendmail api
my $format = "html";
if($text) {
     $format = "text";
}
my @recipients = split(/,/, $to);
my @ccList = split(/,/, $cc);
my @bccList = split(/,/, $bcc);
my $data = {
     from => $from,
     to => \@recipients,
     cc => \@ccList,
     bcc => \@bccList,
     subject => $subject,
     $format => $message,
     token => $token
};

my $data_json = encode_json($data);

my $browser = LWP::UserAgent->new;
my $req = POST $apiUrl;
$req->header('Content-Type' => 'application/json');
$req->header('Accept' => 'application/json');
#must set this, otherwise, will get 'Content-Length header value was wrong, fixed at...' warning
$req->header( "Content-Length" => length($data_json) );
$req->content($data_json);
print "$data_json\n";
my $response = $browser->request($req);
my $result_json = $response->decoded_content;
print $result_json,"\n";

sub Usage {
print <<"END";
     usage: perl $0 [options]
          -from     required, email sender, must be a valid gmail address if a recipient is a gmail address
          -to       required, email recipient(s), Comma separated list of recipients email addresses
          -cc       email recipient(s), Comma separated list of recipients email addresses
          -bcc      email recipient(s), Comma separated list of recipients email addresses
          -subject  required, email subject
          -message  required, email body
          -text     message format is plain text, default is html
END
exit;
} 

