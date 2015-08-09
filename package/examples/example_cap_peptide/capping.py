#!/usr/bin/env python
# coding=utf-8
"""
:Author: Jan Domanski
:Year: 2014
:Copyright: GNU Public License v3

MDAnalysis example: Cap peptide termini
========================================================

This example uses the Merge functionality to add ACE and NMA/NME cap to ends 
of a simple AAQAA3 peptide. This functionality is present in commonly used
packages such as Schrodinger's maestro. 

.. Note::
    Cap structure files are included in the script. 
    The script assumes that aaqaa.gro contains only protein. 
    The caps are added un-minimized by a simple procedure, steric clashes and 
    entanglements may occur.  

"""

from MDAnalysis import Universe, Merge
from MDAnalysis.analysis.align import alignto
import numpy as np


def capping(ref, ace, nma, output):
    resids = ref.select_atoms("all").resids
    resid_min, resid_max = min(resids), max(resids)

    for a in ace.atoms:
        a.resid += resid_min - max(ace.atoms.resids)
    for a in nma.atoms:
        a.resid = resid_max

    # TODO pick the first residue in the protein (how should we cap the chains?)
    # TODO consider a case when the protein resid is 1 and all peptide has to be shifted by +1, put that in docs as a
    #  post-processing step
    alignto(
        ace,
        ref,
        select={
            "mobile": "resid {} and backbone".format(resid_min),
            "reference": "resid {} and backbone".format(resid_min)}
    )
    alignto(
        nma,
        ref,
        select={
            "mobile": "resid {} and backbone and not (resname NMA or resname NME)".format(resid_max),
            "reference": "resid {} and (backbone or name OT2)".format(resid_max)}
    )

    # TODO remove the Hydrogen closest to ACE's oxygen
    u = Merge(
        ace.select_atoms("resname ACE"),
        ref.select_atoms("not (resid {} and name HT*) and not (resid {} and (name HT* or name OT1))".format(resid_min,
                                                                                                           resid_max)),
        nma.select_atoms("resname NME or resname NMA"))
    u.trajectory.ts.dimensions = ref.trajectory.ts.dimensions
    u.atoms.write(output)
    return u


ref = Universe("aaqaa.gro")
ace = Universe("ace.pdb")
nma = Universe("nma.pdb")

capping(ref, ace, nma, "merge.gro")
