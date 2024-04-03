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

"""
FHI-AIMS file format --- :mod:`MDAnalysis.coordinates.FHIAIMS`
==============================================================

Classes to read and write `FHI-AIMS`_ coordinate files.

The cell vectors are specified by the (optional) lines with the
``lattice_vector`` tag::

    lattice_vector x  y  z

where x, y, and z are expressed in ångström (Å).

.. Note::

   In the original `FHI-AIMS format`_, up to three lines with
   ``lattice_vector`` are allowed (order matters) where the absent line implies
   no periodicity in that direction.  In MDAnalysis, only the case of no
   ``lattice_vector`` or three ``lattice_vector`` lines are allowed.

Atomic positions and names are specified either by the ``atom`` or by the
``atom_frac`` tags::

    atom           x  y  z  name
    atom_frac      nx ny nz name

where x, y, and z are expressed in ångström, and nx, ny and nz are real numbers
in [0, 1] and are used to compute the atomic positions in units of the basic
cell.

Atomic velocities can be added on the line right after the corresponding
``atom`` in units of Å/ps using the ``velocity`` tag::

    velocity      vx vy vz


The field name is a string identifying the atomic species.  See also the
specifications in the official `FHI-AIMS format`_.

Classes
-------

.. autoclass:: FHIAIMSReader
   :members:
   :inherited-members:

.. autoclass:: FHIAIMSWriter
   :members:
   :inherited-members:

Developer notes: ``FHIAIMSWriter`` format strings
-------------------------------------------------

The :class:`FHIAIMSWriter` class has a :attr:`FHIAIMSWriter.fmt`
attribute, which is a dictionary of different strings for writing
lines in ``.in`` files.  These are as follows:

``xyz``
  An atom line without velocities.  Requires that the `name` and
  `pos` keys be supplied.  E.g.::

     fmt['xyz'].format(pos=(0.0, 1.0, 2.0), name='O')

``vel``
  An line that specifies velocities::

     fmt['xyz'].format(vel=(0.1, 0.2, 0.3))

``box_triclinic``
  The (optional) initial lines of the file which gives box dimensions.
  Requires the `box` keyword, as a length 9 vector.  This is a flattened
  version of the (3, 3) triclinic vector representation of the unit
  cell.


.. Links

.. _FHI-AIMS: https://aimsclub.fhi-berlin.mpg.de/
.. _FHI-AIMS format: https://doi.org/10.6084/m9.figshare.12413477.v1

"""
import re

import itertools
import warnings

import numpy as np

from . import base
from .core import triclinic_box, triclinic_vectors
from ..exceptions import NoDataError
from ..lib import util
from ..lib import mdamath


class FHIAIMSReader(base.SingleFrameReaderBase):
    """Reader for the FHIAIMS geometry format.

       Single frame reader for the `FHI-AIMS`_ input file format.  Reads
       geometry (3D and molecules only), positions (absolut or fractional),
       velocities if given, all according to the `FHI-AIMS format`_
       specifications

    """
    format = ['IN', 'FHIAIMS']
    units = {'time': 'ps', 'length': 'Angstrom', 'velocity': 'Angstrom/ps'}

    def _read_first_frame(self):
        with util.openany(self.filename, 'rt') as fhiaimsfile:
            relative, positions, velocities, lattice_vectors = [], [], [], []
            skip_tags = ["#", "initial_moment"]
            oldline = ''
            for line in fhiaimsfile:
                line = line.strip()
                if line.startswith("atom"):
                    positions.append(line.split()[1:-1])
                    relative.append('atom_frac' in line)
                    oldline = line
                    continue
                if line.startswith("velocity"):
                    if not 'atom' in oldline:
                        raise ValueError(
                            'Non-conforming line (velocity must follow atom): ({0})in FHI-AIMS input file {0}'.format(line, self.filename))
                    velocities.append(line.split()[1:])
                    oldline = line
                    continue
                if line.startswith("lattice"):
                    lattice_vectors.append(line.split()[1:])
                    oldline = line
                    continue
                if any([line.startswith(tag) for tag in skip_tags]):
                    oldline = line
                    continue
                raise ValueError(
                    'Non-conforming line: ({0})in FHI-AIMS input file {0}'.format(line, self.filename))

        # positions and velocities are lists of lists of strings; they will be
        # cast to np.arrays(..., dtype=float32) during assignment to ts.positions/ts.velocities
        lattice_vectors = np.asarray(lattice_vectors, dtype=np.float32)

        if len(velocities) not in (0, len(positions)):
            raise ValueError(
                'Found incorrect number of velocity tags ({0}) in the FHI-AIMS file, should be {1}.'.format(
                    len(velocities), len(positions)))
        if len(lattice_vectors) not in (0, 3):
            raise ValueError(
                'Found partial periodicity in FHI-AIMS file. This cannot be handled by MDAnalysis.')
        if len(lattice_vectors) == 0 and any(relative):
            raise ValueError(
                'Found relative coordinates in FHI-AIMS file without lattice info.')

        # create Timestep

        self.n_atoms = n_atoms = len(positions)
        self.ts = ts = self._Timestep(n_atoms, **self._ts_kwargs)
        ts.positions = positions

        if len(lattice_vectors) > 0:
            ts.dimensions = triclinic_box(*lattice_vectors)
            ts.positions[relative] = np.matmul(ts.positions[relative], lattice_vectors)

        if len(velocities) > 0:
            ts.velocities = velocities

        self.ts.frame = 0  # 0-based frame number

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Returns a FHIAIMSWriter for *filename*.

        Parameters
        ----------
        filename: str
            filename of the output FHI-AIMS file

        Returns
        -------
        :class:`FHIAIMSWriter`

        """
        if n_atoms is None:
            n_atoms = self.n_atoms
        return FHIAIMSWriter(filename, n_atoms=n_atoms, **kwargs)


class FHIAIMSWriter(base.WriterBase):
    """FHI-AIMS Writer.

       Single frame writer for the `FHI-AIMS`_ format.  Writes geometry (3D and
       molecules only), positions (absolut only), velocities if given, all
       according to the `FHI-AIMS format`_ specifications.

       If no atom names are given, it will set each atom name to "X".

    """

    format = ['IN', 'FHIAIMS']
    units = {'time': None, 'length': 'Angstrom', 'velocity': 'Angstrom/ps'}

    #: format strings for the FHI-AIMS file (all include newline)
    fmt = {
        # coordinates output format, see  https://doi.org/10.6084/m9.figshare.12413477.v1
        'xyz': "atom {pos[0]:12.8f} {pos[1]:12.8f} {pos[2]:12.8f} {name:<3s}\n",
        'vel': "velocity {vel[0]:12.8f} {vel[1]:12.8f} {vel[2]:12.8f}\n",
        # unitcell
        'box_triclinic': "lattice_vector {box[0]:12.8f} {box[1]:12.8f} {box[2]:12.8f}\nlattice_vector {box[3]:12.8f} {box[4]:12.8f} {box[5]:12.8f}\nlattice_vector {box[6]:12.8f} {box[7]:12.8f} {box[8]:12.8f}\n"
    }

    def __init__(self, filename, convert_units=True, n_atoms=None, **kwargs):
        """Set up the FHI-AIMS Writer

        Parameters
        -----------
        filename : str
            output filename
        n_atoms : int (optional)
            number of atoms

        """
        self.filename = util.filename(filename, ext='.in', keep=True)
        self.n_atoms = n_atoms

    def _write_next_frame(self, obj):
        """Write selection at current trajectory frame to file.

        Parameters
        -----------
        obj : AtomGroup or Universe


        .. versionchanged:: 2.0.0
           Support for deprecated Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
        # write() method that complies with the Trajectory API
        try:
            # make sure to use atoms (Issue 46)
            ag = obj.atoms
            # can write from selection == Universe (Issue 49)
        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        # Check for topology information
        missing_topology = []
        try:
            names = ag.names
        except (AttributeError, NoDataError):
            names = itertools.cycle(('X',))
            missing_topology.append('names')

        try:
            atom_indices = ag.ids
        except (AttributeError, NoDataError):
            atom_indices = range(1, ag.n_atoms+1)
            missing_topology.append('ids')

        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "".format(miss=', '.join(missing_topology)))

        positions = ag.positions
        try:
            velocities = ag.velocities
            has_velocities = True
        except (AttributeError, NoDataError):
            has_velocities = False

        with util.openany(self.filename, 'wt') as output_fhiaims:
            # Lattice
            try:  # for AtomGroup/Universe
                tri_dims = obj.universe.trajectory.ts.triclinic_dimensions
            except AttributeError:  # for Timestep
                tri_dims = obj.triclinic_dimensions
            # full output
            if tri_dims is not None:
                output_fhiaims.write(
                    self.fmt['box_triclinic'].format(box=tri_dims.flatten()))

            # Atom descriptions and coords
            # Dont use enumerate here,
            # all attributes could be infinite cycles!
            for atom_index, name in zip(
                    range(ag.n_atoms), names):
                output_fhiaims.write(self.fmt['xyz'].format(
                    pos=positions[atom_index],
                    name=name))
                if has_velocities:
                    output_fhiaims.write(self.fmt['vel'].format(
                        vel=velocities[atom_index]))
