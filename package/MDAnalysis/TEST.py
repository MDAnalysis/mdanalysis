#!/usr/bin/env python3
## File for testing changes 
## Will get deleted when PR is DONE

import shutil
import os

from pathlib import Path
from MDAnalysis.fetch.fetchers import StaticFetcher



shutil.rmtree('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs', ignore_errors=True)
#test_path = Path('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs')
#test_db = test_path / 'test.txt'

downloader = StaticFetcher()
path = downloader.fetch(base_url='https://files.wwpdb.org/download/', file_name='1AKE.pdb', db_name='test.txt')

print(path)
print('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs')
print(os.listdir('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs'))


print('trying cache')
path = downloader.fetch(base_url='https://files.wwpdb.org/download/', file_name='1AKE.pdb', db_name='test.txt')
