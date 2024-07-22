# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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

"""MMTF trajectory reader --- :mod:`MDAnalysis.coordinates.MMTF`
================================================================

Reads coordinates data from the Macromolecular Transmission Format format
(MMTF_). This should generally be a quicker alternative to PDB.

.. versionadded:: 0.16.0
.. versionchanged:: 0.20.0
   Unit cell is now properly treated as optional

Classes
-------

.. autoclass:: MMTFReader
   :members:

.. _MMTF: https://mmtf.rcsb.org/

"""
import warnings
import mmtf

from . import base
from ..core.universe import Universe
from ..lib.util import cached, store_init_arguments
from ..due import due, Doi


def _parse_mmtf(fn):
    if fn.endswith('gz'):
        return mmtf.parse_gzip(fn)
    else:
        return mmtf.parse(fn)


class MMTFReader(base.SingleFrameReaderBase):
    """Coordinate reader for the Macromolecular Transmission Format format (MMTF_).


    .. deprecated:: 2.8.0
       The MMTF format is no longer supported / serviced by the
       Protein Data Bank. The Reader will be removed in version 3.0.
       Users are encouraged to instead use alternative PDB formats.
    """
    format = 'MMTF'

    @store_init_arguments
    def __init__(self, filename, convert_units=True, n_atoms=None, **kwargs):
        wmsg = ("The MMTF Reader is deprecated and will be removed in "
                "MDAnalysis version 3.0.0")
        warnings.warn(wmsg, DeprecationWarning)

        super(MMTFReader, self).__init__(
            filename, convert_units, n_atoms,
            **kwargs
        )

    @staticmethod
    def _format_hint(thing):
        """Can this Reader read *thing*?

        .. versionadded:: 1.0.0
        """
        return isinstance(thing, mmtf.MMTFDecoder)

    @due.dcite(
        Doi('10.1371/journal.pcbi.1005575'),
        description="MMTF Reader",
        path='MDAnalysis.coordinates.MMTF',
    )
    def _read_first_frame(self):
        # TOOD: Check units?
        if isinstance(self.filename, mmtf.MMTFDecoder):
            top = self.filename
        else:
            top = _parse_mmtf(self.filename)
        self.n_atoms = top.num_atoms

        self.ts = ts = self._Timestep(self.n_atoms,
                                      **self._ts_kwargs)
        ts._pos[:, 0] = top.x_coord_list
        ts._pos[:, 1] = top.y_coord_list
        ts._pos[:, 2] = top.z_coord_list
        if not top.unit_cell is None:
            # optional field
            ts.dimensions = top.unit_cell

        return ts
