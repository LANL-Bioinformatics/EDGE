#!/usr/bin/perl

use Getopt::Long;
use strict;

my ($from, $to, $subject, $message);

GetOptions(
           "from=s"        => \$from,
           "to=s"        => \$to,
           "subject=s"        => \$subject,
           "message=s"        => \$message,
           "help|?"       => sub{Usage()} );
          
Usage() unless $from && $to;

open(MAIL, "|/usr/sbin/sendmail -t");

print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject\n";
print MAIL "Content-type: text/html\n\n";

# Email Body
print MAIL $message;
#
close(MAIL);
 
sub Usage {
print <<"END";
     usage: perl $0 [options] 
     		-from 		email sender
            -to         email recipient
            -subject	email subject
            -message	email body
END
exit;
} 

1;          