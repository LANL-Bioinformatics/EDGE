#!/usr/bin/env python

import os
import time
import sys
import argparse as ap
import json
import atexit
from collections import OrderedDict
from selenium import webdriver
#from selenium.common.exceptions import InvalidSessionIdException
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC


def parse_params():
        p = ap.ArgumentParser(prog='NCBI_SARS-CoV2_batch_submitter.py',
                              description="""Batch Upload SARS-CoV2 sequcnes to NCBI""")

        p.add_argument('-u', '--username',
                                   metavar='[STR]', nargs=1, type=str, required=True,
                                   help="NCBI username")

        p.add_argument('-p', '--password',
                                   metavar='[STR]', nargs=1, type=str, required=True,
                                   help="NCBI password")

        p.add_argument('-f', '--fasta',
                                   metavar='[FILE]', type=ap.FileType(), required=True,
                                   help="sequence file in FASTA format")

        p.add_argument('-s', '--source',
                                        metavar='[FILE]', type=ap.FileType(), required=True,
                                        help='source modifier metadata file for sample')
                                        
        p.add_argument('-c', '--comment',
                                        metavar='[FILE]', type=ap.FileType(), required=True,
                                        help='comment metadata file for sample')
                                        
        p.add_argument('-a', '--authors',
                                        metavar='[FILE]', type=ap.FileType(), required=True,
                                        help='authors metadata file for sample')

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

def parseTableToDict(file,extDict=OrderedDict()):
        header=file.readline().strip().split("\t")
        for line in file:
                rows = line.strip().split("\t")
                if rows[0] not in extDict.keys():
                    extDict[rows[0]]=dict()
                for i in range(1,len(header)):
                        if 'human' in rows[i].lower() and header[i] == 'host':
                                extDict[rows[0]][header[i]]='Home sapiens'
                        else:
                                extDict[rows[0]][header[i]]=rows[i]    
        return extDict
        
def parseMetadata(metadata):
        """metadata FileType input"""
        """parse metadata from EDGE """
        meta = {}
        for line in metadata:
                (key, value) = line.strip().split("=")
                meta[key] = value
        return meta

def parseFasta(fasta):
        """fasta FileType input"""
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

def waiting_table_to_get_ready(wait, sec=1):
        """wait for the table to be loaded"""
        wait.until(EC.invisibility_of_element_located(
                (By.XPATH,  "//tbody[@class='yui-dt-message']")))
        time.sleep(sec)
        
def set_viewport_size(driver, width, height):
    window_size = driver.execute_script("""
        return [window.outerWidth - window.innerWidth + arguments[0],
          window.outerHeight - window.innerHeight + arguments[1]];
        """, width, height)
    driver.set_window_size(*window_size)

def quit_driver(driver,outdir):
        screenshot= driver.get_screenshot_as_file(outdir+"/exit_ncbi_screenshot.png")
        driver.quit()

def fill_NCBI_upload(uname, upass, seqfile, source, comment , outdir, authorsMetadata, to, rt, iv, headless,debug):

        seqfile_abs = os.path.abspath(seqfile.name)
        s_metadata = parseTableToDict(source)
        # merge two dictionary
        metadata = parseTableToDict(comment,s_metadata)
        
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
        set_viewport_size(driver, 1280, 1280)
        ## quit the browser if there is any raised exception
        atexit.register(quit_driver,driver,outdir)

        # driverwait
        driver.implicitly_wait(20)
        wait = WebDriverWait(driver, to)

        # open NCBI
        print("Opening website NCBI...")
        driver.get('https://www.ncbi.nlm.nih.gov/account/?back_url=https%3A//submit.ncbi.nlm.nih.gov/sarscov2/')
        time.sleep(iv)
        print(driver.title)
        assert 'NCBI' in driver.title

        # login
        print("Logining to NCBI...")
        wait.until(EC.frame_to_be_available_and_switch_to_it((By.ID,'loginframe')))
        username = driver.find_element_by_name('uname')
        username.send_keys(uname)
        password = driver.find_element_by_name('upasswd')
        password.send_keys(upass)
        singInButton = driver.find_element_by_name('signinBtn')
        singInButton.click()

        driver.switch_to.default_content()
        time.sleep(iv)


        # navigate to  GenBank
        print("Navigating to GenBank Submission Page...")
        genbank_tab = wait.until(EC.element_to_be_clickable((By.XPATH,"//a[@href='/subs/genbank/']")))
        genbank_tab.click()
        time.sleep(iv)

        new_submission = wait.until(EC.element_to_be_clickable((By.ID,'id_sub_new')))
        new_submission.click()
        time.sleep(iv)


        try:
                createAnyway = driver.find_element_by_id('create-new-submission-anyway').click()
                time.sleep(iv)
        except:
                pass

        # Step 1: Submission Type
        # click on SARS-CoV-2 type
        print("1. Submission Type")
        wait.until(EC.text_to_be_present_in_element((By.XPATH,"//h2"),"Submission Type"))
        driver.find_elements_by_xpath("//ul[@id='id_subtype-subgentype']/li")[2].click()
        driver.find_elements_by_xpath("//div[@id='field_subgroup_viruses-choice']//li")[0].click()
        driver.find_element_by_name('subtitle-title').send_keys(metadata[list(metadata)[0]]['isolate'])
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(iv)

        # Step 2: Submitter
        # this page should pull from user profile, no fill action required, press continue
        print("2. Submitter")
        wait.until(EC.text_to_be_present_in_element((By.XPATH,"//h2"),"Submitter"))
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(iv)

        # Step 3: Sequencing Technology
        print("3. Sequencing Technology")
        wait.until(EC.text_to_be_present_in_element((By.XPATH,"//h2"),"Sequencing Technology"))
        seqPlatform=['minion','nanopore','gridion']
        
        if 'illumina' in metadata[list(metadata)[0]]['Sequencing Technology'].lower():
                driver.find_elements_by_xpath("//ul[@id='id_seqtech_method-methods']//label")[3].click()

        if any(x in metadata[list(metadata)[0]]['Sequencing Technology'].lower() for x in seqPlatform):
                driver.find_elements_by_xpath("//ul[@id='id_seqtech_method-methods']//label")[7].click()
                method_detail = wait.until(
                                                EC.presence_of_element_located(
                                                        (By.NAME, 'seqtech_method-other_detail'))
                                        )
                method_detail.send_keys(metadata[list(metadata)[0]]['Sequencing Technology'])
                
        driver.find_elements_by_xpath("//ul[@id='id_seqtech_assembly_state-assembly_state']//label")[1].click()
        i=0
        for key, value in metadata.items(): 
                asm_method=metadata[key]['Assembly Method'].strip().split(" ")
                program_field = 'seqtech_assembly_program-' + str(i) + '-name'
                version_field = 'seqtech_assembly_program-' + str(i) + '-version'
                driver.find_element_by_name(program_field).send_keys(asm_method[1])
                driver.find_element_by_name(version_field).send_keys(asm_method[2])
                driver.find_elements_by_xpath( "//p[@class='add-row']/a")[0].click()
                i += 1
                time.sleep(0.5)
        
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(iv)

        # Step 4: Sequences
        print("4. Sequences")
        wait.until(EC.text_to_be_present_in_element((By.XPATH,"//h2"),"Sequences"))
        driver.find_elements_by_xpath("//ul[@id='id_release-choice']//label")[0].click()
        if "release" in authorsMetadata and authorsMetadata['release']:
            driver.find_elements_by_xpath("//ul[@id='id_release-choice']//label")[1].click()
            driver.find_element_by_name('release-release_date').send_keys(authorsMetadata['release'])
        driver.find_element_by_name('sequences_file-sequences_file').send_keys(seqfile_abs)
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "b.filesize")))
        time.sleep(iv)
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(10)
        try:
            warning = driver.find_elements_by_xpath("//div[contains(@class,'warning')]")[0]
            if warning.is_displayed:
                screenshot= driver.get_screenshot_as_file(outdir+"/seq_warning_or_error.png")
        except:
            pass
            
        try:
            driver.find_elements_by_xpath("//ul[@id='id_internal_ns-gap_type']//label")[0].click()
        except:
            pass
                
        retry = 0
        while retry <= rt and 'Source Modifiers' != driver.find_elements_by_xpath("//h2")[0].text:
                try:
                        driver.find_elements_by_xpath("//ul[@id='id_internal_ns-gap_type']//label")[0].click()
                except:
                        pass    
                driver.find_element_by_id('id_sub_continue').click()
                time.sleep(10)
                retry += 1

        # Step 5: Source Modifiers
        print("5. Source Modifiers")
        ### batch
        try:
            # > 1 sequences
            driver.find_elements_by_xpath("//ul[@id='id_1toall_style-srcmods_style']//label")[1].click()
        except:
            pass
    
        try:
            # only one sequence
            driver.find_element_by_id('id_switch_to_single_seq_file_upload').click()
        except:
	        pass 
        
        srcmodsFile=os.path.abspath(source.name)
        driver.find_element_by_name('srcmods_file-srcmods_file').send_keys(srcmodsFile)
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "b.filesize")))
        time.sleep(iv)
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(10)
        try:
            warning = driver.find_elements_by_xpath("//div[contains(@class,'warning')]")[0]
            if warning.is_displayed:
                screenshot= driver.get_screenshot_as_file(outdir+"/srcmods_warning_or_error.png")
        except:
            pass
        
        retry = 0
        while retry <= rt and 'References' != driver.find_elements_by_xpath("//h2")[0].text:
                driver.find_element_by_id('id_sub_continue').click()
                time.sleep(10)
                retry += 1

        # Step 6: References
        print("6. References")
        authors = authorsMetadata['authors'].strip().split(",")
        for i, name in enumerate(authors):
                firstname_field = 'sequence_author-' + str(i) + '-first_name'
                lastname_field = 'sequence_author-' + str(i) + '-last_name'
                fullname = name.strip().split(" ")
                driver.find_element_by_name(firstname_field).send_keys(fullname[0])
                driver.find_element_by_name(lastname_field).send_keys(fullname[1])
                driver.find_elements_by_xpath( "//p[@class='add-row']/a")[0].click()
                time.sleep(0.5)

        driver.find_elements_by_xpath("//ul[@id='id_reference-publication_status']//label")[0].click()
        driver.find_element_by_name('reference-reference_title').send_keys(metadata[list(metadata)[0]]['isolate'])
        driver.find_elements_by_xpath("//ul[@id='id_select_reference_authors-same_reference_and_sequence_authors']//label")[0].click()
        driver.find_element_by_id('id_sub_continue').click()
        time.sleep(10)

        if debug:
            screenshot= driver.get_screenshot_as_file(outdir+"/submit_ncbi_screenshot.png")
            driver.find_element_by_id('delete_submission_wizard').click()
            obj = driver.switch_to.alert
            obj.accept()
            return

        # Step 7
        print("7. Review and Submit")
        wait.until(EC.text_to_be_present_in_element((By.XPATH,"//h2"),"Submit"))
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
                driver.find_element_by_id('id_sub_submit').click()
                time.sleep(iv)
                screenshot= driver.get_screenshot_as_file(outdir+"/submit_ncbi_screenshot.png")

        # delete submission button
        #driver.find_element_by_id('delete_submission_wizard').click() 

        # close driver  No need close driver here since quit_driver function done when program exit
        # driver.quit()


def main():
        argvs = parse_params()
        print(f"--- Ingest at {time.strftime('%Y-%m-%d %H:%M:%S')} ---")

        authors_metadata = parseMetadata(argvs.authors)
        outdir = os.path.dirname(argvs.source.name)
        
        fill_NCBI_upload(
                argvs.username,
                argvs.password,
                argvs.fasta,
                argvs.source,
                argvs.comment,
                outdir,
                authors_metadata,
                argvs.timeout,
                argvs.retry,
                argvs.interval,
                argvs.headless,
                argvs.debug,
        )
        print("NCBI submit Completed.")


if __name__ == "__main__":
        main()
