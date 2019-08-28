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


"""Core functions of MDAnalysis
============================

The basic class is an :class:`~MDAnalysis.core.groups.AtomGroup`; the whole
simulation is called the
:class:`~MDAnalysis.core.universe.Universe`. Selections are computed on an
:class:`~MDAnalysis.core.groups.AtomGroup` and return another
:class:`~MDAnalysis.core.groups.AtomGroup`.

To get started, load the Universe::

  u = Universe(topology_file, trajectory_file)

A simple selection of all water oxygens within 4 A of the protein::

  water_shell = u.select_atoms('name OH2 and around 4.0 protein')
  water_shell.n_atoms           # how many waters were selected
  water_shell.total_mass()       # their total mass

:class:`~MDAnalysis.core.groups.AtomGroup` instances have various methods that
allow calculation of simple properties. For more complicated analysis, obtain
the coordinates as a numpy array ::

  coords = water_shell.positions

and write your own Python code.


.. _flags-label:

Flags
-----

.. deprecated:: 0.16.2
   The flags registry will be removed in release 1.0.
   Use keyword arguments for functions to obtain the desired behavior.
   See issue `#782 <https://github.com/MDAnalysis/mdanalysis/issues/782>`_
   for more details.

(This is an advanced topic and can probably be skipped by most people.)

There are a number flags that influence how MDAnalysis behaves. They are accessible
through the pseudo-dictionary :data:`MDAnalysis.core.flags`.

The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

 print(MDAnalysis.core.flags.doc())


List of MDAnalysis flags with default values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: flagsDocs


Classes
~~~~~~~

.. data:: flags

.. autoclass:: Flags
   :members:

.. autoclass:: Flag
   :members:

"""
from __future__ import absolute_import

import six
import warnings
_DEPMSG = {
    'use_KDTree_routines': ('This flag is obsolete, all selections now '
                            'automatically select the fastest method'),
    'use_periodic_selections': 'Use periodic=True/False to select_atoms',
    'convert_lengths': 'This behaviour has been deprecated',
    'length_unit': 'This behaviour has been deprecated',
    'time_unit': 'This behaviour has been deprecated',
    'speed_unit': 'This behaviour has been deprecated',
    'force_unit': 'This behaviour has been deprecated',
    'charge_unit': 'This behaviour has been deprecated',
    'use_pbc': 'Supply a pbc kwarg to the relevant methods',
}

__all__ = ['AtomGroup', 'Selection']


# set up flags for core routines (more convoluted than strictly necessary but
# should be clean to add more flags if needed)

class Flags(dict):
    """Global registry of flags. Acts like a dict for item access.

    There are a number flags defined that influence how MDAnalysis behaves. They are
    accessible through the pseudo-dictionary

      :data:`MDAnalysis.core.flags`

    The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
    raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

      print(MDAnalysis.core.flags.__doc__)

    New flags are added with the :meth:`Flags.register` method which takes a new :class:`Flag`
    instance as an argument.

    .. deprecated:: 0.16.2
       The flags registry will be removed in release 1.0.
       Use keyword arguments for functions to obtain the desired behavior.
       See issue `#782 <https://github.com/MDAnalysis/mdanalysis/issues/782>`_
       for more details.

    """

    def __init__(self, *args):
        """For DEVELOPERS: Initialize Flags registry with a *list* of :class:`Flag` instances."""
        super(Flags, self).__init__([(flag.name, flag) for flag in args])

    def get_flag(self, name):
        return super(Flags, self).__getitem__(name)

    def doc(self):
        """Shows doc strings for all flags."""
        return "\n\n".join([flag.__doc__ for flag in self._itervalues()])

    def register(self, flag):
        """Register a new Flag instance with the Flags registry."""
        super(Flags, self).__setitem__(flag.name, flag)

    def update(self, *flags):
        """Update Flags registry with a list of Flag instances."""
        super(Flags, self).update([(flag.name, flag) for flag in flags])

    def setdefault(self, k, d=None):
        raise NotImplementedError

    def __getitem__(self, name):
        return self.get_flag(name).get()

    def __setitem__(self, name, value):
        self.get_flag(name).set(value)

    def _itervalues(self):
        return six.itervalues(super(Flags, self))

    def _items(self):
        return super(Flags, self).items()

    def itervalues(self):
        for flag in self._itervalues():
            yield flag.value

    def iteritems(self):
        for flag in self._itervalues():
            yield flag.name, flag.value

    def values(self):
        return [flag.value for flag in self._itervalues()]

    def items(self):
        return [(flag.name, flag.value) for flag in self._itervalues()]

    def __repr__(self):
        return str(self.items())


class FlagsDynamicDocs(Flags):
    # docs are generated on the fly for interactive use; but because
    # this does not work well with the sphinx documentation system
    # ("AttributeError: 'property' object has no attribute
    # 'expandtabs'") we split the class...
    @property
    def __doc__(self):
        # generate dynamic docs on all flags
        return self.doc()


class IdentityMapping(dict):
    def __getitem__(self, key):
        return key


class Flag(object):
    """A Flag, essentially a variable that knows its default and legal values.

    .. deprecated:: 0.16.2
       The flags registry will be removed in release 1.0.
       Use keyword arguments for functions to obtain the desired behavior.
       See issue `#782 <https://github.com/MDAnalysis/mdanalysis/issues/782>`_
       for more details.
    """

    def __init__(self, name, default, mapping=None, doc=None):
        """Create a new flag which will be registered with Flags.

        Parameters
        ----------
        name : str
            name of the flag, must be a legal python name
        default
            default value
        mapping : dict
            dict that maps allowed input values to canonical values;
            if ``None`` then no argument checking will be performed and
            all values are directly set.
        doc : str
            doc string; may contain string interpolation mappings for::

                    %%(name)s        name of the flag
                    %%(default)r     default value
                    %%(value)r       current value
                    %%(mapping)r     mapping

            Doc strings are generated dynamically and reflect the current state.


        Example
        -------
        Create a new flag::

            newflag = Flag(name, default, mapping, doc)

        """
        self.name = name
        self.value = default
        self.default = default
        # {v1:v1,v2:v1,v3:v3, ...} mapping of allowed values to canonical ones
        self.mapping = mapping or IdentityMapping()
        self._doctemplate = "**%(name)s** = *%(value)r*\n" + (doc or "*undocumented flag*")

    def get(self):
        return self.value

    def set(self, value):
        warnings.warn('MDAnalysis.core.flags is deprecated and will be removed in version 1.0. '
                      '' + _DEPMSG[self.name],  # custom message for each flag too
                      DeprecationWarning)

        if value is not None:
            try:
                self.value = self.mapping[value]
            except KeyError:
                raise ValueError("flag must be None or one of " + str(self.mapping.keys()))
        return self.get()

    def prop(self):
        """Use this for ``property(**flag.prop())``"""
        return {'fget': self.get, 'fset': self.set, 'doc': self.__doc__}

    def __repr__(self):
        return """Flag('{name!s}',{value!r})""".format(**self.__dict__)


class _Flag(Flag):
    @property
    def __doc__(self):
        # generate dynamic docs with current values
        return self._doctemplate % self.__dict__


_flags = [
    _Flag(
        'use_periodic_selections',
        True,
        {True: True, False: False},
        """
        Determines if distance selections (AROUND, POINT) respect periodicity.

        >>> flags['%(name)s'] = value

        Values of flag:
         * True     - periodicity is taken into account if supported
         * False    - periodicity is ignored

        The MDAnalysis preset of this flag is %(default)r.
        """
    ),
    _Flag(
        'use_KDTree_routines',
        'fast',
        {True: 'fast', 'fast': 'fast',  # only KDTree if advantageous
         'always': 'always',  # always even if slower (for testing)
         False: 'never', 'never': 'never'},  # never, only use (slower) alternatives
        """
           Determines which KDTree routines are used for distance selections

           >>> flags['%(name)s'] = value

           Values for flag:

           * True, 'fast'   - only use KDTree routines that are typically faster than others
             -               POINT      uses distance matrix routines (with periodicity)
             -               AROUND     uses KDTree routines (always ignores periodicity)
           * 'always'       - always use KDTree routines where available (eg for benchmarking)
           * False, 'never' - always use alternatives

           The preset value for MDAnalysis is %(default)r.

           :mod:`MDAnalysis.lib.KDTree` routines are significantly faster for some distance
           selections. However, they cannot deal with periodic boxes and thus ignore
           periodicity; if periodicity is crucial, disable KDTree routines with

           >>> MDAnalysis.core.flags['use_KDTree_routines'] = False
           """
    ),
    _Flag(
        'convert_lengths',
        True,
        {True: True, False: False},
        """
           Determine if trajectory reader and writer converts length units between native and MDAnalysis default.

           >>> flags['%(name)s'] = value

           Some trajectories are in other length units than the MDAnalysis
           default (see the flag *length_unit*); e.g. Gromacs XTC and TRR
           trajectories are in nm. If ``True`` then coordinates are
           automatically converted, with ``False`` the coordinate values are
           presented as read from the trajectories.

           .. Note:: The conversion of lengths also affects conversion of velocities.
        """
    ),
    _Flag(
        'length_unit',
        'Angstrom',
        {
            'Angstrom': 'Angstrom', 'A': 'Angstrom',
            'nm': 'nm', 'nano meter': 'nm', 'nanometer': 'nm'
        },
        """
           Base unit for lengths (in particular coordinates in trajectories)

           >>> flags['%(name)s'] = value

           .. Warning:: Do not change, only Angstrom fully supported.
        """
    ),
    _Flag(
        'time_unit',
        'ps',
        {
            'ps': 'ps', 'pico second': 'ps', 'picosecond': 'ps',
            'ns': 'ns', 'nano second': 'ns', 'nanosecond': 'ns',
            'fs': 'fs', 'femto second': 'fs', 'femtosecond': 'fs',
            'AKMA': 'AKMA', 'Charmm': 'AKMA',
        },
        """
           Base unit for times (in particular time steps in trajectories)

           >>> flags['%(name)s'] = value

        """
    ),
    _Flag(
        'speed_unit',
        'Angstrom/ps',
        {
            'Angstrom/ps': 'Angstrom/ps', 'A/ps': 'Angstrom/ps',
            'nm/ps': 'nm/ps', 'nanometer/ps': 'nm/ps', 'nanometer/picosecond': 'nm/ps',
            'Angstrom/AKMA': 'Angstrom/AKMA', 'A/AKMA': 'Angstrom/AKMA',
            'pm/ps': 'pm/ps',
            'm/s': 'm/s',
        },
        """
           Base unit for speed (in particular velocities in trajectories)

           >>> flags['%(name)s'] = value

        """
    ),
    _Flag(
        'force_unit',
        'kJ/(mol*Angstrom)',
        {
            'kJ/(mol*Angstrom)': 'kJ/(mol*Angstrom)', 'kJ/(mol*A)': 'kJ/(mol*Angstrom)',
            'kJ/(mol*nm)': 'kJ/(mol*nm)',
            'kcal/(mol*Angstrom)': 'kcal/(mol*Angstrom)', 'kcal/(mol*A)': 'kcal/(mol*Angstrom)'
        },
        """
           Base unit for forces (in particular forces in trajectories)

           >>> flags['%(name)s'] = value

        """
    ),
    _Flag(
        'charge_unit',
        'e',
        {
            'e': 'e', 'electron charge': 'e',
            'Amber': 'Amber',
            'C': 'C', 'Coulomb': 'C', 'As': 'C',
        },
        """
           Base unit for charge

           >>> flags['%(name)s'] = value

        """
    ),
    _Flag(
        'use_pbc',
        False,
        {True: True, False: False},
        """
        Choose whether to consider periodic boundary conditions when
        performing many :class:`MDAnalysis.core.groups.AtomGroup` methods.
        This is set to ``False`` by default but can be enabled with:

        >>> MDAnalysis.core.flags['use_pbc'] = True

        Values for flag:

        * ``True`` - Move all atoms within the primary unit cell before calculation
        * ``False`` - Use coordinates as supplied

        .. Warning::

           Changing this to ``True`` changes the default behaviour of
           commonly used :class:`MDAnalysis.core.groups.AtomGroup` methods
           such as :meth:`MDAnalysis.core.groups.AtomGroup.center_of_mass`
           and :meth:`MDAnalysis.core.groups.AtomGroup.center_of_geometry`!
        """),
]

#: Global flag registry for :mod:`MDAnalysis.core`.
#: Can be accessed like a dictionary and appears to the casual user as such.
flags = FlagsDynamicDocs(*_flags)
del _flags


# only for sphinx docs
class flagsDocs(object):
    __doc__ = flags.doc()


from . import groups
from . import selection
from . import AtomGroup
