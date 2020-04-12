#!/usr/bin/env python3

import os
import time
import sys
import argparse as ap
import json
from selenium import webdriver
#from selenium.common.exceptions import InvalidSessionIdException
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def parse_params():
    p = ap.ArgumentParser(prog='gisaid_EpiCoV_uploader.py',
                          description="""Download all EpiCoV sequcnes from GISAID""")

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

    args_parsed = p.parse_args()
    return args_parsed


def fill_EpiCoV_upload(uname, upass, seq, metadata, to, rt, iv, headless):
    """Download sequences and metadata from EpiCoV GISAID"""

    # add sequence to metadata
    metadata["sequence"] = seq

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

    # driverwait
    driver.implicitly_wait(20)
    wait = WebDriverWait(driver, to)

    # open GISAID
    print("Opening website GISAID...")
    driver.get('https://platform.gisaid.org/epi3/frontend')
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
    print("Access uploading page...")
    upload_tab = wait.until(EC.element_to_be_clickable(
        (By.CSS_SELECTOR, 'div.sys-actionbar-action:nth-child(4)')))
    upload_tab.click()
    waiting_sys_timer(wait)

    iframe = None
    retry = 1
    while retry <= rt:
        try:
            iframe = driver.find_element_by_xpath("//iframe")
            if iframe:
                break
            else:
                raise
        except:
            print(f"retrying...#{retry} in {iv} sec(s)")
            if retry == rt:
                print("Failed")
                sys.exit(1)
            else:
                time.sleep(iv)
                retry += 1
    
    driver.switch_to.frame(iframe)

    button = wait.until(
        EC.presence_of_element_located(
            (By.XPATH, "//td[1]"))
    )
    #button = driver.find_element_by_xpath('//td[1]')
    button.click()

    driver.switch_to.default_content()
    waiting_sys_timer(wait)

    # keyword mapping
    entry_keys_mapping = {
        # text
        0  : "virus_name", #Virus name*: hCoV-19/Country/Identifier/2020
        1  : "virus_passage", #Passage details/history*: Example: Original, Vero
        2  : "collection_date", #Collection date* Example: 2020-04-01
        3  : "location", #location*: Continent / Country / Region
        4  : "", #Additional location information: Example: Cave, Live animal market
        5  : "host", #Host*
        6  : "", #Additional host information: Example: Cruise Ship, Convention, Live animal market
        7  : "gender", #Gender*
        8  : "age", #Patient age*
        9  : "", #Patient status: Example: Hospitalized, Released, Live, Deceased, unknown
        10 : "isolation_source", #Specimen source: Example: Nasal
        11 : "", #Outbreak Detail: Example: Date, Place, Family cluster
        12 : "", #Last vaccinated
        13 : "", #Treatment: Example: Include drug name, dosage
        14 : "", #Sequencing technology: Nanopore MinION
        15 : "assembly_method", #Assembly method
        16 : "coverage", #Coverage
        17 : "", #Sample ID given by the provider
        18 : "", #Sample ID given by the Submitting lab
        # textarea
        19 : "originating_lab", #Originating lab*
        20 : "originating_address", #Originating lab address*
        21 : "submitting_lab", #Submitting lab*: Los Alamos National Lab
        22 : "submitting_address", #Submitting lab address*
        23 : "authors", #Authors*
        24 : "", #Submitter information: address
        25 : "sequence" #custom
    }

    # fill the webform
    text_inputs = driver.find_elements_by_xpath("//input[@type='text']")
    textareas = driver.find_elements_by_xpath("//textarea")

    num = 0
    for inputs in text_inputs, textareas:
        for text_input in inputs:
            meta_key = entry_keys_mapping[num]
            if meta_key and meta_key in metadata:
                text_input.send_keys(metadata[meta_key])
            num += 1
    
    waiting_sys_timer(wait)

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
        #print("Get Button")
        button = driver.find_element_by_xpath("//button[contains(text(), 'Submit for Review')]")
        button.click()
        waiting_sys_timer(wait)
        
        warnings = driver.find_elements_by_xpath( "//div[@class='sys-form-fi-message']")
        for msg in warnings:
            if msg.is_displayed():
                print(warnings.text)

    # close driver
    driver.quit()


def parseMetadata(metadata):
    """parse metadata from EDGE """
    meta = {}
    for line in metadata:
        (key, value) = line.strip().split("=")
        meta[key] = value

    return meta

def parseFasta(fasta):
    """parse metadata from EDGE """
    seq = ""
    header_found = False
    for line in fasta:
        if line.startswith(">"):
            if header_found:
                print("ERROR: Only allow 1 sequence.")
                sys.exit(1)
            else:
                header_found = True
        else:
            line = line.strip()
            seq += line

    return seq


def waiting_sys_timer(wait, sec=1):
    """wait for system timer"""
    wait.until(EC.invisibility_of_element_located(
        (By.XPATH,  "//div[@id='sys_timer']")))
    time.sleep(sec)


def waiting_table_to_get_ready(wait, sec=1):
    """wait for the table to be loaded"""
    wait.until(EC.invisibility_of_element_located(
        (By.XPATH,  "//tbody[@class='yui-dt-message']")))
    time.sleep(sec)


def download_finished(file, timeout=60):
    sec = 0
    while sec < timeout:
        if os.path.exists(file):
            return True
        else:
            sec += 1
    return False


def main():
    argvs = parse_params()
    print(f"--- Ingest at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")

    seq = parseFasta(argvs.fasta)
    metadata = parseMetadata(argvs.metadata)

    fill_EpiCoV_upload(
        argvs.username,
        argvs.password,
        seq,
        metadata,
        argvs.timeout,
        argvs.retry,
        argvs.interval,
        argvs.headless,
    )
    print("Completed.")


if __name__ == "__main__":
    main()
