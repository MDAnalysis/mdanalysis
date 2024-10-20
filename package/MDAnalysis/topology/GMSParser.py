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

"""
GAMESS Topology Parser
======================

.. versionadded:: 0.9.1

Reads a GAMESS_ output file (also Firefly_ and `GAMESS-UK`_) and pulls
element information from it.  Symmetrical assembly is read (not
symmetry element!).  Atom names are read from the GAMESS section.  Any
information about residues or segments will not be populated.

.. _GAMESS: http://www.msg.ameslab.gov/gamess/
.. _Firefly: http://classic.chem.msu.su/gran/gamess/index.html
.. _`GAMESS-UK`: http://www.cfs.dl.ac.uk/


Classes
-------

.. autoclass:: GMSParser
   :members:
   :inherited-members:

"""
import re
import numpy as np

from ..lib.util import openany
from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Resids,
    Resnums,
    Segids,
    AtomAttr,
)

class AtomicCharges(AtomAttr):
    attrname = 'atomiccharges'
    singular = 'atomiccharge'
    per_object = 'atom'


class GMSParser(TopologyReaderBase):
    """GAMESS_ topology parser.

    Creates the following Attributes:
     - names
     - atomic charges


    .. note::

        By default, atomtypes and masses will be guessed on Universe creation.
        This may change in release 3.0.
        See :ref:`Guessers` for more information.

    .. versionadded:: 0.9.1
    .. versionchanged:: 2.8.0
        Removed type and mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

    """
    format = 'GMS'

    def parse(self, **kwargs):
        """Read list of atoms from a GAMESS file."""
        names = []
        at_charges = []

        with openany(self.filename) as inf:
            while True:
                line = inf.readline()
                if not line:
                    raise EOFError
                if re.match(r'^\s+ATOM\s+ATOMIC\s+COORDINATES\s*\(BOHR\).*',\
                        line):
                    break
            line = inf.readline() # skip

            while True:
                line = inf.readline()
                _m = re.match(\
r'^\s*([A-Za-z_][A-Za-z_0-9]*)\s+([0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+)\s+(\-?[0-9]+\.[0-9]+).*',
                        line)
                if _m is None:
                    break
                name = _m.group(1)
                at_charge = int(float(_m.group(2)))

                names.append(name)
                at_charges.append(at_charge)
                #TODO: may be use coordinates info from _m.group(3-5) ??

        n_atoms = len(names)
        attrs = [
            Atomids(np.arange(n_atoms) + 1),
            Atomnames(np.array(names, dtype=object)),
            AtomicCharges(np.array(at_charges, dtype=np.int32)),
            Resids(np.array([1])),
            Resnums(np.array([1])),
            Segids(np.array(['SYSTEM'], dtype=object)),
        ]
        top = Topology(n_atoms, 1, 1,
                       attrs=attrs)

        return top
