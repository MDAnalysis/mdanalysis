#!/usr/bin/env python
"""Simple script that processes input files to be usable with CHARMM

usage: pdb2charmm file [file ...]

The input file can be other structure files but PDB works best because
it typically provides chain/segment identifiers. CHARMM CRD files with
segments also work well.
"""

import sys
import MDAnalysis
from MDAnalysis.builder.charmm import Preprocessor

import logging

logger = logging.getLogger('MDAnalysis.app')

MDAnalysis.start_logging()

for filename in sys.argv[1:]:
    logger.info("Processing %r", filename)
    P = Preprocessor(filename)
    P.fix_names()
    P.split_chains(filename, format="crd")

MDAnalysis.stop_logging()
