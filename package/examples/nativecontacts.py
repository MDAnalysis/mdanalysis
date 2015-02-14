#!/usr/bin/env python
# -*- coding: utf-8 -*-
# MDAnalysis example: native contact analysis (q1-q2)

"""
Example: Calculating the percentage of native contacts
======================================================

The percentag of native contacts, q, is the number of observed contacts
relative to the number of contacts in a reference state (typically, the
native state). The number of  *contacts* is the total number of unique CA atoms
within in a fixed radius of all CA atoms.

This example script takes a trajectory of AdK in which the enzyme transitions
from the closed to the open state. Two sets of native contacts are computed:

*q1*
   contacts relative to the initial (closed) state, pdb 1AKE
*q2*
   contacts relative to the final (open) state, pdb 4AKE

The trajectory was generated with the dynamic importance sampling method (DIMS)
and taken from [Beckstein2009].

The final result is stored in a bzipped data file `adk_dims_q1q2.dat.bz2` and
the q1-q2 graph is also written to `figures/nativecontacts.pdf`
(requires :mod:`matplotlib`).

[Beckstein2009] O. Beckstein, E.J. Denning, J.R. Perilla and T.B. Woolf,
                Zipping and Unzipping of Adenylate Kinase: Atomistic Insights
                into the Ensemble of Open ↔ Closed Transitions. J Mol Biol 394
                (2009), 160–176, doi:10.1016/j.jmb.2009.09.009


.. SeeAlso:: :mod:`MDAnalysis.analysis.contacts`

"""
import MDAnalysis.analysis.contacts
from MDAnalysis.tests.datafiles import PSF, DCD  # AdK example trajectory

try:
    import matplotlib

    matplotlib.use('agg')  # no interactive plotting, only save figures
    from pylab import savefig, clf

    have_matplotlib = True
except ImportError:
    have_matplotlib = False

C = MDAnalysis.analysis.contacts.ContactAnalysis(PSF, DCD, targetdir="./output")
C.run()

print "Data file was written to %r" % C.output_bz2

if have_matplotlib:
    matplotlib.rc('font', size=14)
    matplotlib.rc('figure', figsize=(3.4, 3.8))
    clf()
    C.plot(linewidth=3)

    savefig("./figures/nativecontacts.pdf")
    savefig("./figures/nativecontacts.png")

    print "Wrote figures to ./figures/nativecontacts.{pdf,png}"
