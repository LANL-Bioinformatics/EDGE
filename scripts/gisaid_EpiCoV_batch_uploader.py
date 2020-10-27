#!/usr/bin/env python

import os
import time
import sys
import argparse as ap
import json
import atexit
from selenium import webdriver
#from selenium.common.exceptions import InvalidSessionIdException
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def parse_params():
        p = ap.ArgumentParser(prog='gisaid_EpiCoV_batch_uploader.py',
                                                  description="""SARS-CoV2 sequcnes submit to GISAID""")

        p.add_argument('-u', '--username',
                                   metavar='[STR]', nargs=1, type=str, required=True,
                                   help="GISAID username")

        p.add_argument('-p', '--password',
                                   metavar='[STR]', nargs=1, type=str, required=True,
                                   help="GISAID password")

        p.add_argument('-f', '--fasta',
                                   metavar='[FILE]', type=ap.FileType(), required=True,
                                   help="sequence file in FASTA format")

        p.add_argument('-m', '--metadata',
                                        metavar='[FILE]', type=ap.FileType(), required=True,
                                        help='metadata file for sample')

        p.add_argument('-t', '--timeout',
                                   metavar='[INT]', type=int, required=False, default=90,
                                   help="set action timeout seconds. Default is 90 secs.")

        p.add_argument('-r', '--retry',
                                   metavar='[INT]', type=int, required=False, default=5,
                                   help="retry how many times when the action fails. Default is 5 times.")

        p.add_argument('-i', '--interval',
                                   metavar='[INT]', type=int, required=False, default=3,
                                   help="time interval between retries in second(s). Default is 3 seconds.")

        p.add_argument('--headless',
                                   action='store_true', help='turn on headless mode')

        p.add_argument('--debug',
                                   action='store_true', help='turn on debug mode')

        p.add_argument('--version', action='version', version='%(prog)s 1.0.0')

        args_parsed = p.parse_args()
        return args_parsed

def quit_driver(driver,outdir):
        screenshot= driver.get_screenshot_as_file(outdir+"/exit_gisaid_screenshot.png")
        driver.quit()

def waiting_sys_timer(wait, sec=1):
        """wait for system timer"""
        wait.until(EC.invisibility_of_element_located(
                (By.XPATH,  "//div[@id='sys_timer']")))
        time.sleep(sec)

def fill_EpiCoV_upload(uname, upass, seq, metadatafile, to, rt, iv, headless, debug):
        
        outdir= os.path.dirname(metadatafile.name)
        # MIME types
        mime_types = "application/octet-stream"
        mime_types += ",application/excel,application/vnd.ms-excel"
        mime_types += ",application/pdf,application/x-pdf"

        print("Opening browser...")
        profile = webdriver.FirefoxProfile()
        profile.set_preference("browser.download.folderList", 2)
        profile.set_preference("browser.download.manager.showWhenStarting", False)
        profile.set_preference(
                "browser.helperApps.neverAsk.saveToDisk", mime_types)
        profile.set_preference(
                "plugin.disable_full_page_plugin_for_types", mime_types)
        profile.set_preference("pdfjs.disabled", True)

        options = Options()
        if headless:
                options.add_argument("--headless")

        driver = webdriver.Firefox(firefox_profile=profile, options=options)

        ## quit the browser if there is any raised exception
        atexit.register(quit_driver,driver,outdir)
        
        # driverwait
        driver.implicitly_wait(20)
        wait = WebDriverWait(driver, to)

        # open GISAID
        print("Opening website GISAID...")
        driver.get('https://www.epicov.org/epi3/frontend')
        waiting_sys_timer(wait)
        print(driver.title)
        assert 'GISAID' in driver.title

        # login
        print("Logining to GISAID...")
        username = driver.find_element_by_name('login')
        username.send_keys(uname)
        password = driver.find_element_by_name('password')
        password.send_keys(upass)
        driver.execute_script("return doLogin();")

        waiting_sys_timer(wait)

        # navigate to EpiFlu
        print("Navigating to EpiCoV...")
        epicov_tab = driver.find_element_by_xpath("//div[@id='main_nav']//li[3]/a")
        epicov_tab.click()

        waiting_sys_timer(wait)

        # access uploading page
        # WARNING: different users might have different uploading options
        print("Accessing batch uploading page...")
        try:
                batch_upload_tab = driver.find_element_by_xpath('//div[@class=sys-actionbar-action][contains(text(), "Batch Upload")]')
                batch_upload_tab.click()
                waiting_sys_timer(wait)
        except:
                upload_tab = wait.until(EC.element_to_be_clickable(
                        (By.CSS_SELECTOR, 'div.sys-actionbar-action:nth-child(4)')))
                upload_tab.click()
                waiting_sys_timer(wait)
                
        
        try:
                time.sleep(iv)      
                wait.until(EC.presence_of_element_located((By.XPATH, "//iframe")))
                iframe = driver.find_element_by_xpath("//iframe")
                if iframe.is_displayed() and iframe.get_attribute('id').startswith('sysoverlay'):
                        print("Popup window detected...")
                        driver.switch_to.frame(iframe)
                        button = wait.until(
                                EC.presence_of_element_located(
                                        (By.XPATH, "//td[2]"))
                        )
                
                        print("Choosing batch upload option...")
                        #button = driver.find_element_by_xpath('//td[1]')
                        button.click()

                        driver.switch_to.default_content()
                        waiting_sys_timer(wait)
        except:
                pass

        print("Send Excel metadata file...")
        iframe = iframe=driver.find_elements_by_tag_name('iframe')[0]
        driver.switch_to.frame(iframe)

        button = wait.until(
                                        EC.presence_of_element_located(
                                                (By.XPATH, "//input[@name='data'][@type='file']"))
                                )
        button.send_keys(metadatafile.name)
        time.sleep(iv)
        driver.switch_to.default_content()
        waiting_sys_timer(wait)

        print("Send Sequence fasta file...")
        iframe2 = iframe=driver.find_elements_by_tag_name('iframe')[1]
        driver.switch_to.frame(iframe2)

        button = wait.until(
                                        EC.presence_of_element_located(
                                                (By.XPATH, "//input[@name='data'][@type='file']"))
                                )
        button.send_keys(seq.name)
        time.sleep(iv)
        driver.switch_to.default_content()
        waiting_sys_timer(wait)

        if debug:
            screenshot= driver.get_screenshot_as_file(outdir+"/submit_gisaid_screenshot.png")
            return

        if not headless:
                # wait until the user to close browser
                print("Please review the form and submit for review...")
                while True:
                        try:
                                _ = driver.window_handles
                        except:
                                print("Browser closed by user.")
                                break
                        time.sleep(1)
        else:
                button = driver.find_element_by_xpath('//button[contains(text(), "Check and Submit")]')
                button.click()
                waiting_sys_timer(wait)

                warnings = driver.find_elements_by_xpath( "//div[@class='bd']")
                for msg in warnings:
                        if msg.is_displayed():
                                print(msg.text)
                
                reportMSG = driver.find_elements_by_xpath( "//div[@class='sys-form-fi-multiline-ro']")
                for msg in reportMSG:
                        print(msg.text)

                screenshot= driver.get_screenshot_as_file(outdir+"/submit_screenshot.png")
                        
        # close driver  No need close driver here since quit_driver function done when program exit
        #driver.quit()

def main():
        argvs = parse_params()
        print(f"--- Ingest at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")

        fill_EpiCoV_upload(
                argvs.username,
                argvs.password,
                argvs.fasta,
                argvs.metadata,
                argvs.timeout,
                argvs.retry,
                argvs.interval,
                argvs.headless,
                argvs.debug
        )
        print("GISAID submit Completed.")


if __name__ == "__main__":
        main()
