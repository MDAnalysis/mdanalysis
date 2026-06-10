#!/usr/bin/env python3
## File for testing changes 
## Will get deleted when PR is DONE

import shutil

from pathlib import Path
from MDAnalysis.fetch.fetchers import StaticFetcher



#shutil.rmtree('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs')
test_path = Path('/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs')
test_db = test_path / 'test.txt'

s = StaticFetcher()
t = s.fetch(base_url='https://files.wwpdb.org/download/')

print(t)

