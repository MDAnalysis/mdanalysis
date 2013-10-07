# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; encoding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# Convolve MTSS rotamers with MD trajectory.
# Copyright (c) 2011-2013 Philip Fowler, Oliver Beckstein
# Published under the GNU Public Licence, version 2 (or higher)
#
# Includes a rotamer library for MTSS at 298 K by Gunnar Jeschke,
# which is published under the same licence by permission.

import MDAnalysis
import MDAnalysis.analysis.align
import MDAnalysis.KDTree.NeighborSearch as KDNS
from MDAnalysis.core.distances import distance_array
import MDAnalysis.analysis.distances

import numpy as np
import os.path

import logging
logger = logging.getLogger("MDAnalysis.app")


import library


def rms_fit_trj(*args, **kwargs):
    """Silenced :func:`MDAnalysis.analysis.align.rms_fit_trj`"""
    kwargs['quiet'] = True
    return MDAnalysis.analysis.align.rms_fit_trj(*args, **kwargs)



class RotamerDistances(object):
    """Calculation of distance distributions between two spin labels."""
    def __init__(self, *args, **kwargs):
        """RotamerDistances(universe, residue_list, **kwargs)

        :Arguments:
           *universe*
              :class:`MDAnalysis.Universe`
           *residue_list*
              list of two residue numbers ``(r1, r2)`` that indicate
              the labelled sites

        :Keywords:
           *dcdFilename*
              name of the temporary files with rotamers fitted [``'trj'``]
           *outputFile*
              stem of the name of the file containing the distance histogram
              (the final name will be ``<outputFile><resid_1>-<resid_2>.dat``
              [``'distances'``]
           *libname*
              library name; the library is loaded with
              :class:`rotamers.library.RotamerLibrary` [``'MTSSL 298K'``]
           *discardFrames*
              skip initial frames < *discardFrames* [``0``]
           *clashDistance*
              discard rotamer if any distance between rotamer atoms
              and protein atoms is < *clashDistance*. Values down to
              1.5 Å are reasonable. The default is conservative. [``2.2`` Å]
           *histogramBins*
             tuple ``(dmin, dmax, delta)`` in Ångström [``(0.0, 100.0, 1.0)``]
        """
        proteinStructure = args[0]
        residues = args[1]
        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 residue numbers: current value {0}.".format(residues))

        outputFile,ext = os.path.splitext(kwargs.pop('outputFile', 'distances'))
        ext = ext or ".dat"
        self.outputFile = "{0}-{1[0]}-{1[1]}{2}".format(outputFile, residues, ext)

        dcdFilename,ext = os.path.splitext(kwargs.pop('dcdFilename', 'trj'))
        ext = ext or ".dcd"
        tmptrj = ["{0}-1{1}".format(dcdFilename, ext),  # or make this temp files?
                  "{0}-2{1}".format(dcdFilename, ext),  # or make this temp files?
                  ]

        kwargs.setdefault('discardFrames', 0)
        self.clashDistance = kwargs.pop('clashDistance', 2.2)  # Ångström
        histogramBins = kwargs.pop('histogramBins', (0.0, 100.0, 1.0))

        self.lib = library.RotamerLibrary(kwargs.get('libname', 'MTSSL 298K'))

        # setup the main lists
        distances = []
        weights = []

        logger.info("Starting rotamer distance analysis of trajectory {0}...".format(proteinStructure.trajectory.filename))
        logger.info("clashDistance = {0} A; rotamer library = '{1}'".format(self.clashDistance, self.lib.name))
        logger.debug("Temporary trajectories for rotamers 1 and 2 (only last frame of MD trajectory): {0[0]} and {0[1]}".format(tmptrj))
        logger.debug("Results will be written to {0}.".format(self.outputFile))

        progressmeter = MDAnalysis.core.log.ProgressMeter(proteinStructure.trajectory.numframes, interval=1)
        for protein in proteinStructure.trajectory:
            progressmeter.echo(protein.frame)
            if protein.frame < kwargs['discardFrames']:
                continue
            # define the atoms used to fit the rotamers. Note that an
            # ordered list has to be created as the ordering of C CA N is
            # different in both. Fit the rotamers onto the protein:
            self.fit_rotamers(self.lib.rotamers,proteinStructure,residues[0],tmptrj[0])
            rotamersSite1 = MDAnalysis.Universe(self.lib.rotamers.filename,tmptrj[0])
            (rotamer1_clash,rotamer1_clash_total) = self.find_clashing_rotamers(rotamersSite1,proteinStructure,residues[0])

            self.fit_rotamers(self.lib.rotamers,proteinStructure,residues[1],tmptrj[1])
            rotamersSite2 = MDAnalysis.Universe(self.lib.rotamers.filename,tmptrj[1])
            (rotamer2_clash,rotamer2_clash_total) = self.find_clashing_rotamers(rotamersSite2,proteinStructure,residues[1])

            # define the atoms to measure the distances between
            rotamer1nitrogen = rotamersSite1.selectAtoms("name N1")
            rotamer2nitrogen = rotamersSite2.selectAtoms("name N1")

            # loop over all the rotamers on the first site
            for rotamer1 in rotamersSite1.trajectory:
                # only proceed if it isn't the PDB frame and there is no clash
                if rotamer1.frame > 1 and not rotamer1_clash[rotamer1.frame]:
                    # loop over all the rotamers on the second site
                    for rotamer2 in rotamersSite2.trajectory:
                        # only proceed if it isn't the PDB frame and there is no clash
                        if rotamer2.frame > 1 and not rotamer2_clash[rotamer2.frame]:
                            # measure and record the distance
                            (a,b,distance) = MDAnalysis.analysis.distances.dist(rotamer1nitrogen,rotamer2nitrogen)
                            distances.append(distance[0])
                            # create the weights list
                            try:
                                # subtract 1 because frame is 1-based
                                weight = self.lib.weights[rotamer1.frame-1] * self.lib.weights[rotamer2.frame-1]
                            except IndexError:
                                logger.error("oppps: no weights for rotamer 1 #{0} - rotamer 2 #{1}".format(rotamer1.frame, rotamer2.frame))
                            weights.append(weight)

        # check that at least two distances have been measured
        if len(distances) < 2:
            logger.critical("no distances found between the spin pair!")
            raise RuntimeError("no distances found between the spin pair!") # should this really be an exception?

        # calculate Nbins and min and max so that we cover at least
        # the requested lower and upper bounds with the given fixed
        # bin width
        bins=MDAnalysis.core.util.fixedwidth_bins(histogramBins[2], histogramBins[0], histogramBins[1])
        # use numpy to histogram the distance data, weighted appropriately
        (a,b) = np.histogram(distances,weights=weights,density=True,bins=bins['Nbins'],range=(bins['min'],bins['max']))

        with open(self.outputFile,'w') as OUTPUT:
            for (i,j) in enumerate(a):
                OUTPUT.write("%6.2f %8.3e\n" %  ((0.5*(b[i] + b[i+1])), j))
        logger.info("Distance distribution for residues {0[0]} - {0[1]} was written to {1}".format(residues, self.outputFile))

    def plot(self, **kwargs):
        """Load data file and plot"""
        import matplotlib.pyplot as plt
        filename = kwargs.pop('filename', None)
        fig = kwargs.pop('fig', None)
        if fig is None:
            fig = plt.figure(figsize=(5,5))
        ax = kwargs.pop('ax', None)
        if ax is None:
            ax = fig.add_subplot(111)

        dist, prob = np.loadtxt(self.outputFile, unpack=True)
        ax.plot(dist, prob, **kwargs)
        ax.set_xlabel(r"spin-label distance $d$ ($\AA$)")
        ax.set_ylabel("probability density")

        if filename:
            ax.figure.savefig(filename)
            logger.info("Plotted distance distribution to {0}".format(filename))

        return ax

    def fit_rotamers(self,rotamers,protein,site_resid,dcdfile):
        """Produce a temporary trajectory of the rotamers.

        The backbone of the rotamers is fitted to the backbone of the
        spin labelled residues.
        """
        # create an ordered list allowing the rotamer to be fitted onto the backbone of the protein
        fittingSelection = (["name C","name CA","name N"],
                            ["protein and name C and resid {0}".format(site_resid),
                             "protein and name CA and resid {0}".format(site_resid),
                             "protein and name N and resid {0}".format(site_resid)
                             ])
        # fit the rotamer library onto the protein
        rms_fit_trj(rotamers,protein,select=fittingSelection,mass_weighted=True,filename=dcdfile)
        return dcdfile

    def find_clashing_rotamers(self,fitted_rotamers,protein,site_resid):
        """Detect any rotamer that clashes with the protein."""
        # make a KD tree of the protein neighbouring atoms
        proteinNotSite = protein.selectAtoms("protein and not name H* and not (resid " + str(site_resid) + " or (resid " + str(site_resid-1) + " and (name C or name O)) or (resid " + str(site_resid+1) + " and name N))")
        proteinNotSiteLookup = KDNS.AtomNeighborSearch(proteinNotSite)

        rotamerSel = fitted_rotamers.selectAtoms("not name H*")

        rotamer_clash = [True]
        rotamer_clash_counter=0
        for rotamer in fitted_rotamers.trajectory:
            bumps = proteinNotSiteLookup.search_list(rotamerSel,self.clashDistance)
            if bumps:
                rotamer_clash.append(True)
                rotamer_clash_counter+=1
            else:
                rotamer_clash.append(False)
        return (rotamer_clash, rotamer_clash_counter)

