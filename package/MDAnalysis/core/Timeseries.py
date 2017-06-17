# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""Compute observable timeseries from trajectories --- :mod:`MDAnalysis.core.Timeseries`
=======================================================================================

.. deprecated:: 0.16.2
   The Timeseries functionality (in particular correl) will be removed in the
   0.17 release. See issue `#1372
   <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more details.
   To extract coordinates efficiently one can use
   :class:`MDAnalysis.analysis.base.AnalysisFromFunction` as shown in
   :ref:`Creating an in-memory trajectory with AnalysisFromFunction
   <create-in-memory-trajectory-with-AnalysisFromFunction>`.


The collection of timeseries (such as :class:`Atom`, :class:`Bond`,
:class:`Dihedral`...)  can be computed from a trajectory in one go, foregoing
the need to iterate through the trajectory frame by frame in python. Inspired
by CHARMM's correl command.

The disadvantage is that the timeseries 'plugins' must be implemented
in C-code. Hence not all trajectory readers (see :mod:`MDAnalysis.coordinates`)
support them.


Basic classes
-------------

.. autoclass:: Timeseries
   :members:

.. autoclass:: TimeseriesCollection
   :members:


Timeseries of observables
-------------------------

.. autoclass:: Atom
.. autoclass:: Bond
.. autoclass:: Angle
.. autoclass:: Dihedral
.. autoclass:: Distance
.. autoclass:: CenterOfGeometry
.. autoclass:: CenterOfMass
.. autoclass:: WaterDipole

"""
from __future__ import division, absolute_import
import warnings

from numpy.lib.utils import deprecate

from . import groups

with warnings.catch_warnings():
    warnings.simplefilter('always', DeprecationWarning)
    warnings.warn(('The Timeseries module and TimeseriesCollection will be '
                   'removed in the 0.17 release. See issue #1372 '
                   'https://github.com/MDAnalysis/mdanalysis/issues/1373'),
                  DeprecationWarning)


class TimeseriesCollection(object):
    '''A collection of timeseries objects.

    The collection of timeseries (such as Atom, Bond, Dihedral...)  can be
    computed from a trajectory in one go, foregoing the need to iterate through
    the trajectory frame by frame in python. Inspired by CHARMM's correl
    command.

    The disadvantage is that the timeseries 'plugins' must be implemented in
    C-code. ::

       collection = TimeseriesCollection()
       collection.addTimeseries(Timeseries.Atom(...)) - add a new Timeseries object
       collection.compute(...)                        - compute the collection of timeseries from the trajectory
       collection.clear()                             - clear the collection
       collection[i]                                  - access the i'th timeseries
       len(collection)                                - return the number of Timeseries added to the collection


    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self):
        self.timeseries = []

    def __len__(self):
        return len(self.timeseries)

    def __getitem__(self, item):
        if (type(item) == int) or (type(item) == slice):
            return self.timeseries[item]
        else:
            raise IndexError

    def __repr__(self):
        if len(self) != 1:
            suffix = 's'
        else:
            suffix = ''
        return '<' + self.__class__.__name__ + ' with {0:d} timeseries object{1!s}>'.format(len(self), suffix)

    def addTimeseries(self, ts):
        '''add a Timeseries object to the collection'''
        if not isinstance(ts, Timeseries):
            raise Exception("Can only add Timeseries objects to TimeseriesCollection")
        self.timeseries.append(ts)

    def clear(self):
        '''clear the timeseries collection'''
        self.timeseries = []

    def compute(self, trj, start=None, stop=None, skip=None, step=None):
        '''Iterate through the trajectory *trj* and compute the time series.

         *trj*
            dcd trajectory object (i.e. :attr:`Universe.trajectory`)
         *start*
            First frame of trajectory to analyse, Default: None becomes 0.
         *stop*
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number,
            which means that the trajectory would be read until the end.
         *step*
            Step between frames to analyse, Default: None becomes 1.
         *Deprecated*
            Skip is deprecated in favor of step.

        '''
        if skip is not None:
            step = skip
            warnings.warn("Skip is deprecated and will be removed in"
                          "in 1.0. Use step instead.",
                          category=DeprecationWarning)

        start, stop, step = trj.check_slice_indices(start, stop, step)
        self.data = trj.correl(self, start, stop, step)
        # Now remap the timeseries data to each timeseries
        typestr = "|f8"
        start = 0
        for t in self.timeseries:
            finish = t.getDataSize()
            datasize = len(t.getFormatCode())
            subarray = self.data[start:start + finish]
            if datasize != 1:
                subarray.shape = (datasize, subarray.shape[0] // datasize, -1)
            t.__data__ = subarray
            t.__array_interface__ = subarray.__array_interface__
            start += finish

    def _getAtomList(self):
        self._atomlist = []
        for ts in self.timeseries:
            self._atomlist += ts.getAtomList()
        return self._atomlist

    def _getFormat(self):
        return ''.join([ts.getFormatCode() for ts in self.timeseries])

    def _getBounds(self):
        if not hasattr(self, '_atomlist'):
            self.getAtomList()
        if not len(self._atomlist) == 0:
            lowerb = min(self._atomlist)
            upperb = max(self._atomlist)
        else:
            lowerb = upperb = 0
        return lowerb, upperb

    def _getDataSize(self):
        return sum([ts.getDataSize() for ts in self.timeseries])

    def _getAtomCounts(self):
        atomCounts = []
        for ts in self.timeseries:
            atomCounts += ts.getAtomCounts()
        return atomCounts

    def _getAuxData(self):
        auxData = []
        for ts in self.timeseries:
            auxData += ts.getAuxData()
        return auxData


class Timeseries(object):
    '''Base timeseries class - define subclasses for specific timeseries computations.

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, code, atoms, dsize):
        if isinstance(atoms, groups.AtomGroup):
            self.atoms = atoms.atoms
        elif isinstance(atoms, list):
            self.atoms = atoms
        elif isinstance(atoms, groups.Atom):
            self.atoms = [atoms]
        else:
            raise TypeError("Invalid atoms passed to {0!s} timeseries".format(self.__class__.__name__))
        self.code = code
        self.n_atoms = len(self.atoms)
        self.dsize = dsize

    # Properties
    @property
    def shape(self):
        """shape tuple of the underlying numpy array"""
        return self.__data__.shape

    def __getitem__(self, item):
        return self.__data__[item]

    def __len__(self):
        return len(self.__data__)

    def __repr__(self):
        if hasattr(self, "__data__"):
            return '<' + self.__class__.__name__ + ' timeseries object is populated with data>\n{0!s}'.format( \
                                                   (repr(self.__data__)))
        else:
            return '<' + self.__class__.__name__ + ' timeseries object is not populated with data>'

    def getAtomList(self):
        return [a.index for a in self.atoms]

    def getFormatCode(self):
        return self.code

    def getDataSize(self):
        return self.dsize

    def getNumAtoms(self):
        return self.n_atoms

    def getAtomCounts(self):
        return [self.n_atoms]

    def getAuxData(self):
        return [0.] * self.n_atoms


class Atom(Timeseries):
    '''Create a timeseries that returns coordinate data for an atom or group of atoms ::

           t = Atom(code, atoms)

        *code*
          is one of 'x', 'y', 'z', or 'v' ('vector', which returns all three
          dimensions)

        *atoms*
          can be a single :class:`~MDAnalysis.core.groups.Atom` object,
          a list of :class:`~MDAnalysis.core.groups.Atom` objects, or an
          :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, code, atoms):
        if code not in ('x', 'y', 'z', 'v', 'w'):
            raise ValueError("Bad code")
        if code == 'v':
            size = 3
        else:
            size = 1
        if isinstance(atoms, groups.AtomGroup):
            n_atoms = len(atoms.atoms)
        elif isinstance(atoms, list):
            n_atoms = len(atoms)
        elif isinstance(atoms, groups.Atom):
            n_atoms = 1
        else:
            raise TypeError("Invalid atoms passed to {0!s} timeseries".format(self.__class__.__name__))
        Timeseries.__init__(self, code * n_atoms, atoms, size * n_atoms)

    def getAtomCounts(self):
        return [1, ] * self.n_atoms


class Bond(Timeseries):
    '''Create a timeseries that returns a timeseries for a bond

           t = Bond(atoms)

    *atoms* must contain 2 :class:`~MDAnalysis.core.groups.Atom` instances, either as a list or an
    :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        if not len(atoms) == 2:
            raise ValueError("Bond timeseries requires a 2 atom selection")
        Timeseries.__init__(self, 'r', atoms, 1)


class Angle(Timeseries):
    '''Create a timeseries that returns a timeseries for an angle

           t = Angle(atoms)

    atoms must contain 3 :class:`~MDAnalysis.core.groups.Atom` instances,
    either as a list or an :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        if not len(atoms) == 3:
            raise ValueError("Angle timeseries requires a 3 atom selection")
        Timeseries.__init__(self, 'a', atoms, 1)


class Dihedral(Timeseries):
    '''Create a timeseries that returns a timeseries for a dihedral angle

           t = Dihedral(atoms)

    atoms must contain 4 :class:`~MDAnalysis.core.groups.Atom` objects, either
    as a list or an :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        if not len(atoms) == 4:
            raise ValueError("Dihedral timeseries requires a 4 atom selection")
        Timeseries.__init__(self, 'h', atoms, 1)


class Distance(Timeseries):
    '''Create a timeseries that returns distances between 2 atoms

           t = Distance(code, atoms)

    code is one of 'd' (distance vector), or 'r' (scalar distance) atoms must
    contain 2 :class:`~MDAnalysis.core.groups.Atom` objects, either as a list
    or an :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, code, atoms):
        if code not in ('d', 'r'):
            raise ValueError("Bad code")
        if code == 'd':
            size = 3
        else:
            size = 1
        if not len(atoms) == 2:
            raise ValueError("Distance timeseries requires a 2 atom selection")
        Timeseries.__init__(self, code, atoms, size)


class CenterOfGeometry(Timeseries):
    '''Create a timeseries that returns the center of geometry of a group of atoms

           t = CenterOfGeometry(atoms)

    *atoms* can be a list of :class:`~MDAnalysis.core.groups.Atom` objects, or
    a :class:`~MDAnalysis.core.groups.AtomGroup`


    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        Timeseries.__init__(self, 'm', atoms, 3)

    def getAuxData(self):
        return [1.] * self.n_atoms


class CenterOfMass(Timeseries):
    '''Create a timeseries that returns the center of mass of a group of atoms

           t = CenterOfMass(atoms)

    *atoms* can be a list of :class:`~MDAnalysis.core.groups.Atom` objects or a
    :class:`~MDAnalysis.core.groups.AtomGroup`

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        Timeseries.__init__(self, 'm', atoms, 3)

    def getAuxData(self):
        return [a.mass for a in self.atoms]


class WaterDipole(Timeseries):
    r'''Create a Timeseries that returns a timeseries for the bisector vector of a 3-site water

           d = WaterDipole(atoms)

    *atoms* must contain 3 :class:`~MDAnalysis.core.groups.Atom` objects,
    either as a list or an :class:`~MDAnalysis.core.groups.AtomGroup`; the
    first one *must* be the oxygen, the other two are the hydrogens.

    The vector ``d``, multiplied by the partial charge on the oxygen atom
    (e.g. *q* = -0.0.834 for TIP3P water), gives the actual dipole moment.

    The vector is calculated from the positions of the oxygen atom
    (:math:`\mathbf{x}_{\text{O}}`) and the two hydrogen atoms
    (:math:`\mathbf{x}_{\text{H}_1}`, :math:`\mathbf{x}_{\text{H}_2}`) as

    .. math::

       \mathbf{d} = \mathbf{x}_{\text{O}} - \frac{1}{2}(\mathbf{x}_{\text{H}_1} + \mathbf{x}_{\text{H}_2})

    and the dipole moment vector is

    .. math::

       \boldsymbol{\mu} = q_{\text{O}} \mathbf{d}

    .. Note::

       This will only work for water models that have half of the oxygen charge
       on each hydrogen. The vector :math:`\mathbf{d}` has the opposite
       direction of the dipole moment; multiplying with the oxygen charge
       (:math:`q_{\text{O}}<0`) will flip the direction and produce the correct
       orientation.

       There are no sanity checks; *if the first atom in a water
       molecule is not oxygen then results will be wrong.*

    .. deprecated:: 0.16.2
       The Timeseries functionality (in particular correl) will be removed in
       the 0.17 release. See issue `#1372
       <https://github.com/MDAnalysis/mdanalysis/issues/1373>`_ for more
       details. Use :class:`MDAnalysis.analysis.base.AnalysisFromFunction` instead.

    '''

    @deprecate(message="This class will be removed in 0.17")
    def __init__(self, atoms):
        if not len(atoms) == 3:
            raise ValueError("WaterDipole timeseries requires a 3 atom selection")
        Timeseries.__init__(self, 'w', atoms, 3)
