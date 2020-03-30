import tenacity
from apscheduler.schedulers.blocking import BlockingScheduler
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.chrome.options import Options
import time
import logging
import os
import sys

import logging
print(sys.path)
from send_notifications import send_sms_message, send_email_message

sched = BlockingScheduler()
logger = logging.getLogger('watch_bokeh')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('/Devel/edge/bokeh_watcher.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

@tenacity.retry(wait=tenacity.wait_fixed(10),
                    retry=tenacity.retry_if_exception_type(IOError))
@sched.scheduled_job('interval', seconds=900)
def run_check():
    CHROME_PATH = '/usr/bin/google-chrome'
    CHROMEDRIVER_PATH = '/usr/local/bin/chromedriver'
    WINDOW_SIZE = "1920,1080"

    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--window-size=%s" % WINDOW_SIZE)

    driver = webdriver.Chrome(CHROMEDRIVER_PATH, chrome_options=chrome_options)
    try:

        pangia_url = 'http://edge-prod.lanl.gov/pangia-vis?r=pangia-vis/data/b456d8c4464a0ab4f45037d786be21ed.tsv'

        driver.get(pangia_url)
        if driver.title == '503 Service Unavailable':
            print(driver.title)
            msg = driver.title
            send_sms_message(msg)
            logger.debug(msg)
            exit(0)
        element = WebDriverWait(driver, 10).until(
            lambda driver: driver.find_element_by_tag_name('h1'))
        title = driver.find_element_by_tag_name('h1')
        print(title.text)
        logger.debug('Pinged bokeh server')
    finally:
        driver.close()
@sched.scheduled_job('interval', seconds=288000)
def email_log():
    send_sms_message('Still pinging...')

sched.start()
