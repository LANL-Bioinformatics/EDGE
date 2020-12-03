import os
from requests import Request, Session
from twilio.rest import Client
from twilio.http import HttpClient
from twilio.http.response import Response
import logging
import logging.config

# LOGGING_CONF = os.path.join(os.path.dirname(__file__), "docker_deploy", "logging_config.ini")
#
# logging.config.fileConfig(LOGGING_CONF)
# logger = logging.getLogger(__name__)

class ProxiedTwilioHttpClient(HttpClient):
    """
    General purpose HTTP Client for interacting with the Twilio API
    """
    def request(self, method, url, params=None, data=None, headers=None, auth=None, timeout=None,
                allow_redirects=False):

        session = Session()
        session.proxies = {
                              "https" : "http://proxyout.lanl.gov:8080"
                          }

        request = Request(method.upper(), url, params=params, data=data, headers=headers, auth=auth)

        prepped_request = session.prepare_request(request)
        response = session.send(
            prepped_request,
            allow_redirects=allow_redirects,
            timeout=timeout,
        )

        return Response(int(response.status_code), response.content.decode('utf-8'))


def send_sms_message(message):
    # the following line needs your Twilio Account SID and Auth Token
    client = Client("AC17c6a64e13a57006df7415a2abeddbce", "a3d357b0f257695697d95f641d63ced5", http_client=ProxiedTwilioHttpClient())

    # change the "from_" number to your Twilio number and the "to" number
    # to the phone number you signed up for Twilio with, or upgrade your
    # account to send SMS to any phone number
    client.messages.create(to="+15054129358",
                           from_="+15053932147",
                           body=message)

def send_email_message(message):
    SENDMAIL = "/usr/sbin/sendmail"  # sendmail location

    p = os.popen("%s -t" % SENDMAIL, "w")
    p.write("To: mcflynn617@gmail.com\n")
    p.write("Subject: test\n")
    p.write("\n")# blank line separating headers from body
    p.write("{0}\n".format(message))
    sts = p.close()
    # if sts != 0:
    #     logger.info("Sendmail exit status", sts)