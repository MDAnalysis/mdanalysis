# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""GSD topology parser
=========================

.. versionadded:: 0.17.0

The :class:`GSDParser` generates a topology from HOOMD_ GSD topology/trajectory
files. The GSD file stores information on both the topology and the trajectory
in the same file, and allows for varying atom numbers/identities and topologies
during the course of the simulation. At the moment MDAnalysis can deal only
with the case in which there is no variation. The trajectory data are read with
the :class:`~MDAnalysis.coordinates.GSD.GSDReader` class.

.. _HOOMD: http://codeblue.umich.edu/hoomd-blue/index.html
.. _HOOMD GSD: https://github.com/glotzerlab/gsd


To load a GSD HOOMD file::

   import MDAnalysis as mda
   u = mda.Universe("example.gsd")


Classes
-------

.. autoclass:: GSDParser
   :members:
   :inherited-members:

"""
import os
import numpy as np

from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomtypes,
    Atomnames,
    Atomids,
    Angles,
    Bonds,
    Charges,
    Dihedrals,
    Impropers,
    Masses,
    Radii,
    Resids,
    Resnums,
    Resnames,
    Segids,
)

try:
    import gsd.hoomd
except ImportError:
    HAS_GSD = False
else:
    HAS_GSD = True


class GSDParser(TopologyReaderBase):
    """Parses a Hoomd GSD file to create a Topology

    Reads the following Attributes:
     - Atomtypes
     - Bonds
     - Angles
     - Dihedrals
     - Impropers
     - Radii
     - Masses

    The GSD file also stores a `body` property in the particles, and the parser
    uses this information to set the residue names and indices.

    NOTE: if the `body` index of any particle is negative, the parser will add
    an integer number (the absolute value of the minimum of all the body
    indices) to all the body indices. This is because MDAnalysis cannot handle
    negative residue indices. This means that in that case the residue index in
    the MDAnalysis.Universe will not correspond to the body index stored in the
    GSD file.

    """
    format = 'GSD'

    def __init__(self, filename):

        if not HAS_GSD:
            errmsg = ("GSDParser: To read a Topology from a Hoomd GSD "
                      "file, please install gsd")
            raise ImportError(errmsg)
        super(GSDParser, self).__init__(filename)

    def parse(self, **kwargs):
        """Parse Hoomd GSD file

        .. versionadded:: 0.17.0
        """
        attrs = {}

        with gsd.hoomd.open(self.filename, mode='r') as t :
            # Here it is assumed that the particle data does not change in the
            # trajectory.
            snap = t[0]

            natoms = snap.particles.N


            ptypes = snap.particles.types
            atypes = [ptypes[idx] for idx in snap.particles.typeid]
            if len(atypes) != natoms:
                raise IOError("Number of types does not equal natoms.")
            attrs['types'] = Atomtypes(np.array(atypes, dtype=object))

            # set radii, masses, charges
            p = snap.particles
            attrs['diameter'] = Radii(np.array(p.diameter / 2.,dtype=np.float32))
            attrs['mass'] = Masses(np.array(p.mass,dtype=np.float64))
            attrs['charge'] = Charges(np.array(p.charge,dtype=np.float32))

            # set bonds, angles, dihedrals, impropers
            for attrname, attr, in (
                    ('bonds', Bonds),
                    ('angles', Angles),
                    ('dihedrals', Dihedrals),
                    ('impropers', Impropers),
            ):
                try:
                    val = getattr(snap,attrname)
                    vals = [tuple(b_instance) for b_instance in val.group]
                except:
                    vals = []
                attrs[attrname] = attr(vals)

            # get body ids to set residue number and ids
            blist = snap.particles.body.astype(np.int64)
            bodies = np.unique(blist).astype(np.int32)
            # this fixes the fact that the Topology constructor gets stuck in an
            # infinite loop if any resid is negative.
            if (blist<0).any() :
                m = blist.min()
                blist += abs(m)
            bodies = np.unique(blist).astype(np.int32)
            nbodies = bodies.size

        attrs = list(attrs.values())
        attrs.append(Atomnames(np.array(atypes, dtype=object)))
        attrs.append(Atomids(np.arange(natoms) + 1))
        attrs.append(Resids(bodies))
        attrs.append(Resnums(bodies))
        attrs.append(Resnames(bodies))
        attrs.append(Segids(np.array(['SYSTEM'], dtype=object)))

        top = Topology(natoms, nbodies, 1,
                       attrs=attrs, atom_resindex=blist)

        return top
