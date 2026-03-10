# bench_align.py
# Benchmarks for MDAnalysis.analysis.align
# No benchmarks existed for this module before this file

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small


class AlignBenchmark:
    """
    Benchmarks for MDAnalysis.analysis.align module.
    Tests alignto() and AlignTraj() which are the two
    most commonly used functions in this module.
    """

    def setup(self):
        # setup() runs ONCE before any timing starts
        # load data here — never in time_ methods
        # PSF = protein structure file (topology)
        # DCD = trajectory file (coordinates over time)
        # PDB_small = single reference structure
        self.mobile = mda.Universe(PSF, DCD)
        self.reference = mda.Universe(PSF, PDB_small)

    def time_alignto_all_atoms(self):
        """Time aligning all atoms of a single frame"""
        align.alignto(
            self.mobile,
            self.reference,
            select="all"
        )

    def time_alignto_backbone(self):
        """Time aligning backbone atoms only (common in research)"""
        align.alignto(
            self.mobile,
            self.reference,
            select="backbone"
        )

    def time_alignto_mass_weighted(self):
        """Time aligning with mass weighting"""
        align.alignto(
            self.mobile,
            self.reference,
            select="backbone",
            weights="mass"
        )

    def time_aligntraj_all_frames(self):
        """Time aligning entire trajectory — most expensive operation"""
        alignment = align.AlignTraj(
            self.mobile,
            self.reference,
            select="backbone",
            in_memory=True  # avoids writing to disk during benchmark
        )
        alignment.run()

    def time_average_structure(self):
        """Time computing average structure across trajectory"""
        avg = align.AverageStructure(
            self.mobile,
            self.reference,
            select="backbone"
        )
        avg.run()