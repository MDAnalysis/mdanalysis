#!/usr/bin/env python
# Convolve spin label rotamers with MD trajectories.
# Copyright (c) 2011-2013 Philip Fowler
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSSL at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission
"""\
%prog [options] TOPOLOGY TRAJECTORY [TRAJECTORY ...]
       %prog [options] STRUCTURE

Calculate a distribution of spin label distances from an MD trajectory
or an arbitrary ensemble of conformations by fitting Gunnar Jeschke's
rotamer library of MTSS (at 298 K).

Provide either a topology and trajectory data (multiple trajectories
are concatenated on the fly) or just a structure file such as a PDB or
GRO file.

For details on the method see Stelzl et al, J Mol Biol 426 (2014),
735-751, doi:10.1016/j.jmb.2013.10.024.
"""

import os.path

import MDAnalysis
import numpy as np
import time

import rotcon.convolve

import logging
logger = logging.getLogger("MDAnalysis.app")


if __name__ == "__main__":
        from optparse import OptionParser
        parser = OptionParser(usage=__doc__)
        parser.add_option("--resid", type=int, nargs=2, dest="residues", default=None,
                          help ="REQUIRED: the pair of residues to compute DEER distances for")
        parser.add_option("--discard", dest="discardFrames", type=int, default=0,
                          help="discard the first N frames [%default]")
        parser.add_option("--histogramBins", dest="histogramBins", type=float, nargs=3, default=(0,80,1),
                          help ="the start, end and incr distance for the histogram of the distances in Angstroms [%default]")
        parser.add_option("--clashDistance", dest="clashDistance", type=float, default=2.2,
                          help ="the distance between heavy-atoms of the label and protein within which "
                          "they are assumed to be clashing [%default Angstroms]")
        parser.add_option("--output", dest="outputFile", default="output.dat",
                          help ="the path and name of the output histogram file; the filename will "
                          "have resid 1 and resid 2 inserted before the extension [%default]")
        parser.add_option("--dcdfilename", dest="dcdFilename", metavar="FILENAME",
                          help ="the path and stem of the DCD files of the fitted MTSS rotamers")
        parser.add_option("--libname", dest="libname", metavar="NAME", default="MTSSL 298K",
                          help="name of the rotamer library [%default]")
        parser.add_option("--plotname", dest="plotname", metavar="FILENAME", default=None,
                          help="plot the histogram to FILENAME (the extensions determines the format) "
                          "By default <outputFile>.pdf.")


        options, args = parser.parse_args()

        MDAnalysis.start_logging()
        logger.info("Rotamer Convolve MD --- Copyright (c) Philip W Fowler, Oliver Beckstein 2011-2013")
        logger.info("Released under the GNU Public Licence, version 2 (or higher)")
        logger.info("Please cite: LS Stelzl, PW Fowler, MSP Sansom, O Beckstein. J Mol Biol 426 (2014), 735-751, doi:10.1016/j.jmb.2013.10.024")

        # load the reference protein structure
        try:
            proteinStructure = MDAnalysis.Universe(*args)
        except:
            logger.critical("protein structure and/or trajectory not correctly specified")
            raise
        if options.residues is None or len(options.residues) != 2:
            raise ValueError("Provide residue ids in --residues R1 R2")

        logger.info("Loading trajectory data as Universe({0})".format(*args))

        if not options.dcdFilename:
                options.dcdFilename = options.outputFile + "-tmp"

        startTime = time.time()
        R = rotcon.convolve.RotamerDistances(proteinStructure, options.residues,
                                               outputFile=options.outputFile, dcdFilename=options.dcdFilename,
                                               libname=options.libname, discardFrames=options.discardFrames,
                                               clashDistance=options.clashDistance,
                                               histogramBins=options.histogramBins)
        logger.info("DONE with analysis, elapsed time %6i s" % (int(time.time() - startTime)))

        if options.plotname is None:
            root, ext = os.path.splitext(R.outputFile)
            options.plotname = root + ".pdf"
        R.plot(filename=options.plotname, linewidth=2)

        MDAnalysis.stop_logging()


