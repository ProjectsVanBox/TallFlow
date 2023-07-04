#!/usr/bin/python

import time
import urllib
import zipfile36 as zipfile
import glob
import os
import sys

import psutil
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.firefox_profile import FirefoxProfile
from selenium.webdriver.remote.command import Command
from selenium.webdriver.support.ui import Select

def communicate( fasta, output, species = "human", rss_type = "both" ):
        ''' Uploads the generated .fasta file to RSSsite and analyses the sequences for potential RSSs '''

        driver = webdriver.Firefox()
        driver.get('http://www.itb.cnr.it/rss/analyze.html')

        select_species_option = Select(driver.find_element("name", "species"))
        select_species_option.select_by_visible_text( species )

        select_rss_option = Select(driver.find_element("name",'spacer'))
        select_rss_option.select_by_visible_text( rss_type )

        upload_file = driver.find_element("name", 'upfile')
        upload_file.send_keys( fasta )

        analyse_sequence_button = driver.find_element("xpath",'/html/body/table[2]/tbody/tr/td[3]/table/tbody/tr[1]/td/table/tbody/tr[8]/td[2]/input')
        analyse_sequence_button.click()

        time.sleep(5)
        cpu = psutil.cpu_percent()
        while cpu > 20:
            cpu = psutil.cpu_percent()
        
        download_link = driver.find_element("link text",'Click here to get the zipped txt tab separated version of this table')
        download_address = download_link.get_attribute('href')

        urllib.request.urlretrieve(download_address, output)
        driver.close()

        file_to_extract = zipfile.ZipFile( output )
        file_to_extract.extractall()
        files = glob.glob('file*.txt')
        rss_file = files.pop(0)
        os.rename(rss_file, output)

        for f in files:
            os.remove(f)
        
        return rss_file

def count_rss(rss_file):
    '''Counts RSSs with a RIC score above the threshold '''
    def keep_pass(rss_file):
        with open(rss_file) as reader:
            rss = [x for x in reader if 'PASS' in x.split('\t')[-1]]
        return rss

    rss_passes = keep_pass(rss_file)
    rss_unique = {x.split('\t')[0] for x in rss_passes}

    return len(rss_unique)

fasta_file = os.path.realpath( sys.argv[1] )
out_file = fasta_file.replace(".fasta","_rss.txt")
rss_file = communicate( fasta_file, out_file )
#rss_count = count_rss( out_file )
#print( rss_count )