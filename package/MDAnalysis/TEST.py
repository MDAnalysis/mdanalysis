#!/usr/bin/env python3
## File for testing changes 
## Will get deleted when PR is DONE

import shutil
import os

from pathlib import Path
from MDAnalysis.fetch.fetchers import StaticFetcher


DEFAULT_CACHE_FOLDER='/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs'
downloader = StaticFetcher()

## Inital Download Case (create db)
print("Inital Download Case (create db)\n")
shutil.rmtree(DEFAULT_CACHE_FOLDER, ignore_errors=True)
path = downloader.fetch(base_url='https://files.wwpdb.org/download/', file_name='1AKE.pdb', db_name='test.txt')


## Cache Case (has DB)
print('Cache Case (has DB)\n')
downloader = StaticFetcher()
path = downloader.fetch(base_url='https://files.wwpdb.org/download/', file_name='1AKE.pdb', db_name='test.txt')


## No Database (just download)
print('No Database (just download)\n')
shutil.rmtree(DEFAULT_CACHE_FOLDER, ignore_errors=True)
path = downloader.fetch(base_url='https://files.wwpdb.org/download/', file_name='1AKE.pdb', db_name=None)
