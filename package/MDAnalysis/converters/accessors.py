# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

"""AtomGroup accessors --- :mod:`MDAnalysis.lib.accessors`
=============================================================

The aim of this module is simply to link a wrapper class to an object in order
to avoid cluttering the namespace of the object.

Example
-------

.. code-block:: python

    >>> class SpeechWrapper:
    ...     def __init__(self, person):
    ...         self.person = person
    ...     def __call__(self, *args):
    ...         print(self.person.name, "says", *args)
    ...     def whoami(self):
    ...         print("I am %s" % self.person.name)
    ...
    >>> class Person:
    ...     def __init__(self, name):
    ...         self.name = name
    ...     say = Accessor(SpeechWrapper)
    ...
    >>> bob = Person("Bob")
    >>> bob.say("hello")
    Bob says hello
    >>> bob.say.whoami()
    I am Bob

"""

from functools import partial, update_wrapper

from .. import _CONVERTERS
from ..core._get_readers import get_converter_for


class Accessor:
    """Class that adds a wrapper to an object.

    Parameters
    ----------

    wrapper : object
        A wrapper class that implements the `__call__` method and some other
        custom methods.
    """

    def __init__(self, wrapper):
        self._wrapper = wrapper

    def __get__(self, obj, cls):
        return self._wrapper(obj)


class ConverterWrapper:
    """Wrapper class for AtomGroup converters. The converters are accessible
    to any AtomGroup through the `convert_to` property.

    `ag.convert_to` will return this ConverterWrapper, which can be called
    directly with the name of the destination package as a string (similarly
    to the old API), or through custom methods named after the package (in
    lowercase) that are automatically constructed thanks to metaclass magic.

    Parameters
    ----------

    ag : AtomGroup
        AtomGroup that will have access to the wrapper as a property
    """

    def __init__(self, ag):
        self._ag = ag
        for lib, converter in _CONVERTERS.items():
            method_name = lib.lower()
            # create partial function that passes ag to the converter
            fconvert = partial(converter().convert, self._ag)
            # copy docstring and metadata from the converter to fconvert
            update_wrapper(fconvert, converter().convert)
            setattr(self, method_name, fconvert)

    def __call__(self, package, **kwargs):
        """Convert :class:`AtomGroup` to a structure from another Python
        package.

        Example
        -------

        The code below converts a Universe to a :class:`parmed.structure.Structure`.

        .. code-block:: python

            >>> import MDAnalysis as mda
            >>> from MDAnalysis.tests.datafiles import GRO
            >>> u = mda.Universe(GRO)
            >>> parmed_structure = u.atoms.convert_to('PARMED')
            >>> parmed_structure
            <Structure 47681 atoms; 11302 residues; 0 bonds; PBC (triclinic); NOT parametrized>

        Parameters
        ----------
        package: str
            The name of the package to convert to, e.g. ``"PARMED"``
        **kwargs:
            Other parameters passed to the converter

        Returns
        -------
        output:
            An instance of the structure type from another package.

        Raises
        ------
        TypeError:
            No converter was found for the required package


        .. versionadded:: 1.0.0
        """
        converter = get_converter_for(package.upper())
        return converter().convert(self._ag, **kwargs)
