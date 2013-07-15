# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Core functions of MDAnalysis
============================

The basic class is an :class:`~MDAnalysis.core.AtomGroup.AtomGroup`;
the whole simulation is called the
:class:`~MDAnalysis.core.AtomGroup.Universe`. Selections are computed
on an :class:`~MDAnalysis.core.AtomGroup.AtomGroup` and return another
:class:`~MDAnalysis.core.AtomGroup.AtomGroup`.

:mod:`~MDAnalysis.Timeseries` are a convenient way to analyse trajectories.

To get started, load the Universe::

  u = Universe(psffilename,dcdfilename)

A simple selection of all water oxygens within 4 A of the protein::

  water_shell = u.selectAtoms('name OH2 and around 4.0 protein')
  water_shell.numberOfAtoms()   # how many waters were selected
  water_shell.totalMass()       # their total mass

:class:`AtomGroup` instances have various methods that allow
calculation of simple properties. For more complicated analysis,
obtain the coordinates as a numpy array ::

  coords = water_shell.coordinates()

and write your own Python code.


.. _flags-label:

Flags
-----

(This is an advanced topic and can probably be skipped by most people.)

There are a number flags that influence how MDAnalysis behaves. They are accessible
through the pseudo-dictionary :data:`MDAnalysis.core.flags`.

The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

 print MDAnalysis.core.flags.doc()


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

__all__ = ['AtomGroup', 'Selection', 'Timeseries']

# set up flags for core routines (more convoluted than strictly necessary but should
# be clean to add more flags if needed)
class Flags(dict):
    """Global registry of flags. Acts like a dict for item access.

    There are a number flags defined that influence how MDAnalysis behaves. They are
    accessible through the pseudo-dictionary

      :data:`MDAnalysis.core.flags`

    The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
    raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

      print MDAnalysis.core.flags.__doc__

    New flags are added with the :meth:`Flags.register` method which takes a new :class:`Flag`
    instance as an argument.
    """
    def __init__(self,*args):
        """For DEVELOPERS: Initialize Flags registry with a *list* of :class:`Flag` instances."""
        super(Flags,self).__init__([(flag.name,flag) for flag in args])
    def get_flag(self,name):
        return super(Flags,self).__getitem__(name)
    def doc(self):
        """Shows doc strings for all flags."""
        return "\n\n".join([flag.__doc__ for flag in self._itervalues()])
    def register(self,flag):
        """Register a new Flag instance with the Flags registry."""
        super(Flags,self).__setitem__(flag.name,flag)
    def update(self,*flags):
        """Update Flags registry with a list of Flag instances."""
        super(Flags,self).update([(flag.name,flag) for flag in flags])
    def setdefault(self,k,d=None):
        raise NotImplementedError
    def __getitem__(self,name):
        return self.get_flag(name).get()
    def __setitem__(self,name,value):
        self.get_flag(name).set(value)
    def _itervalues(self):
        return super(Flags,self).itervalues()
    def _items(self):
        return super(Flags,self).items()
    def itervalues(self):
        for flag in self._itervalues():
            yield flag.value
    def iteritems(self):
        for flag in self._itervalues():
            yield flag.name,flag.value
    def values(self):
        return [flag.value for flag in self._itervalues()]
    def items(self):
        return [(flag.name,flag.value) for flag in self._itervalues()]
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
    def __getitem__(self,key):
        return key

class Flag(object):
    """A Flag, essentially a variable that knows its default and legal values."""
    def __init__(self,name,default,mapping=None,doc=None):
        """Create a new flag which will be registered with FLags.

          newflag = Flag(name,default,mapping,doc)

        :Arguments:
         *name*
            name of the flag, must be a legal python name
         *default*
            default value
         *mapping*
            dict that maps allowed input values to canonical values;
            if ``None`` then no argument checking will be performed and
            all values are directly set.
         *doc*
            doc string; may contain string interpolation mappings for::

                    %%(name)s        name of the flag
                    %%(default)r     default value
                    %%(value)r       current value
                    %%(mapping)r     mapping

            Doc strings are generated dynamically and reflect the current state.
        """
        self.name = name
        self.value = default
        self.default = default
        # {v1:v1,v2:v1,v3:v3, ...} mapping of allowed values to canonical ones
        self.mapping = mapping or IdentityMapping()
        self._doctemplate = "**%(name)s** = *%(value)r*\n" + (doc or "*undocumented flag*")
    def get(self):
        return self.value
    def set(self,value):
        if value is not None:
            try:
                self.value = self.mapping[value]
            except KeyError:
                raise ValueError("flag must be None or one of "+str(self.mapping.keys()))
        return self.get()
    def prop(self):
        """Use this for ``property(**flag.prop())``"""
        return {'fget':self.get, 'fset':self.set, 'doc':self.__doc__}
    def __repr__(self):
        return """Flag('%(name)s',%(value)r)""" % self.__dict__

class _Flag(Flag):
    @property
    def __doc__(self):
        # generate dynamic docs with current values
        return  self._doctemplate % self.__dict__

_flags = [
    _Flag('use_periodic_selections',
         True,
         {True:True,False:False},
         """
            Determines if distance selections (AROUND, POINT) respect periodicity.

            >>> flags['%(name)s'] = value

            Values of flag:
             * True     - periodicity is taken into account if supported
             * False    - periodicity is ignored

            The MDAnalysis preset of this flag is %(default)r.

            Note that KD-tree based distance selections always ignore this flag. (For
            details see the docs for the 'use_KDTree_routines' flag.)
            """
         ),
    _Flag('use_KDTree_routines',
         'fast',
         {True:'fast','fast':'fast',   # only KDTree if advantageous
          'always':'always',           # always even if slower (for testing)
          False:'never','never':'never'},  # never, only use (slower) alternatives
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

            :mod:`MDAnalysis.KDTree` routines are significantly faster for some distance
            selections. However, they cannot deal with periodic boxes and thus ignore
            periodicity; if periodicity is crucial, disable KDTree routines with

            >>> MDAnalysis.core.flags['use_KDTree_routines'] = False
            """
         ),
    _Flag('convert_gromacs_lengths',
         True,
         {True:True, False:False},
         """
            Determine if the XTC and TRR trajectory reader and writer convert units.

            >>> flags['%(name)s'] = value

            Gromacs trajectories are in nm, MDAnalysis in Angstroem. If
            ``True`` then coordinates are automatically converted, with
            ``False`` the coordinate values are presented as read from the
            trajectories.
         """
         ),
    _Flag('length_unit',
         'Angstrom',
         {'Angstrom': 'Angstrom', 'A':'Angstrom',
          'nm': 'nm', 'nano meter': 'nm', 'nanometer': 'nm'
          },
         """
            Base unit for lengths (in particular coordinates in trajectories)

            >>> flags['%(name)s'] = value

            .. Warning:: Do not change, only Angstrom fully supported.
         """
         ),
    _Flag('time_unit',
         'ps',
         {'ps': 'ps', 'pico second': 'ps', 'picosecond': 'ps',
          'ns': 'ns', 'nano second': 'ns', 'nanosecond': 'ns',
          'fs': 'fs', 'femto second': 'fs', 'femtosecond': 'fs',
          'AKMA': 'AKMA', 'Charmm': 'AKMA',
          },
         """
            Base unit for times (in particular time steps in trajectories)

            >>> flags['%(name)s'] = value

         """
         ),
    _Flag('speed_unit',
         'Angstrom/ps',
         {'Angstrom/ps': 'Angstrom/ps', 'A/ps':'Angstrom/ps',
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
    _Flag('force_unit',
         'kJ/(mol*Angstrom)',
         {'kJ/(mol*Angstrom)': 'kJ/(mol*Angstrom)', 'kJ/(mol*A)': 'kJ/(mol*Angstrom)', 
          'kJ/(mol*nm)': 'kJ/(mol*nm)',
          },
         """
            Base unit for forces (in particular forces in trajectories)

            >>> flags['%(name)s'] = value

         """
         ),
    _Flag('charge_unit',
         'e',
         {'e': 'e', 'electron charge': 'e',
          'Amber': 'Amber',
          'C': 'C', 'Coulomb': 'C', 'As': 'C',
          },
         """
            Base unit for charge

            >>> flags['%(name)s'] = value

         """
         ),
    _Flag('permissive_pdb_reader',
         True,
         {'primitive': True, 'permissive': True, True: True,
          'Bio.PDB': False, 'biopython': False, False: False,
          },
         """
            Select the default reader for PDB Brookhaven databank files.

            >>> flags['%(name)s'] = value

            The Bio.PDB reader (value=``False``) can deal with 'proper' PDB
            files from the Protein Databank that contain special PDB features
            such as insertion codes and it can auto-correct some common
            mistakes; see :mod:`Bio.PDB` for details. However, Bio.PDB has been
            known to read some simulation system PDB files **incompletely**; a
            sure sign of problems is a warning that an atom has appeared twice
            in a residue.

            Therefore, the default for the PDB reader is ``True``, which
            selects the "primitive" (or "permissive") reader
            :class:`MDAnalysis.coordinates.PDB.PrimitivePDBReader`, which
            essentially just reads ATOM and HETATM lines and puts atoms in a
            list.

            One can manually switch between the two by providing the *permissive*
            keyword to :class:`MDAnalysis.Universe`.
         """
         ),
    ]

#: Global flag registry for :mod:`core`.
#: Can be accessed like a dictionary and appears to the casual user as such.
flags = FlagsDynamicDocs(*_flags)
del _flags

# only for sphinx docs
class flagsDocs(object):
    __doc__ = flags.doc()

import AtomGroup, Selection, Timeseries
