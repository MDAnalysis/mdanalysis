#!/usr/bin/env python3

from MDAnalysis.fetch.pdb import from_PDB

import shutil

shutil.rmtree("/nfs/homes3/jauy1/.cache/MDAnalysis_pdbs", ignore_errors=True)

print(from_PDB(['1AKE']))
print(from_PDB(['1AKE', '4AKE']))
