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

"""AtomGroup accessors --- :mod:`MDAnalysis.core.accessors`
====================================================================

This module provides classes for accessing and converting :class:`~MDAnalysis.core.groups.AtomGroup`
objects. It is used for the :meth:`~MDAnalysis.core.groups.AtomGroup.convert_to`
method to make it usable in two different ways: ``ag.convert_to("PACKAGE")`` or
``ag.convert_to.package()``

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
    ...     say = Accessor("say", SpeechWrapper)
    ...
    >>> bob = Person("Bob")
    >>> bob.say("hello")
    Bob says hello
    >>> bob.say.whoami()
    I am Bob

Classes
-------

.. autoclass:: Accessor
   :members:

.. autoclass:: ConverterWrapper
   :members:

"""
from functools import partial, update_wrapper

from .. import _CONVERTERS
from typing import TYPE_CHECKING
from typing import Type, Any

if TYPE_CHECKING:
    from .groups import AtomGroup

    from ..coordinates.base import ConverterBase



class Accessor:
    """Used to pass data between two classes

    Parameters
    ----------
    name : str
        Name of the property in the parent class
    accessor : class
        A class that needs access to its parent's instance

    Example
    -------
    If you want the property to be named "convert_to" in the AtomGroup class,
    use:

    .. code-block:: python

        >>> class AtomGroup:
        >>>     # ...
        >>>     convert_to = Accessor("convert_to", ConverterWrapper)

    And when calling ``ag.convert_to.rdkit()``, the "rdkit" method of the
    ConverterWrapper will be able to have access to "ag"


    .. versionadded:: 2.0.0
    """

    def __init__(self, name: str, accessor: Any) -> None:
        self._accessor = accessor
        self._name = name

    def __get__(self, obj, cls):
        if obj is None:
            # accessing from class instead of instance
            return self._accessor
        # instances the accessor class with the parent object as argument
        wrapped = self._accessor(obj)
        # replace the parent object's property with the wrapped instance
        # so we avoid reinstantiating the accessor everytime `obj.<_name>`
        # is called
        object.__setattr__(obj, self._name, wrapped)
        return wrapped


class ConverterWrapper:
    """Convert :class:`AtomGroup` to a structure from another Python
    package.

    The converters are accessible to any AtomGroup through the ``convert_to``
    property. `ag.convert_to` will return this ConverterWrapper, which can be
    called directly with the name of the destination package as a string
    (similarly to the old API), or through custom methods named after the
    package (in lowercase) that are automatically constructed thanks to
    metaclass magic.

    Example
    -------
    The code below converts a Universe to a :class:`parmed.structure.Structure`

    .. code-block:: python

        >>> import MDAnalysis as mda
        >>> from MDAnalysis.tests.datafiles import GRO
        >>> u = mda.Universe(GRO)
        >>> parmed_structure = u.atoms.convert_to('PARMED')
        >>> parmed_structure
        <Structure 47681 atoms; 11302 residues; 0 bonds; PBC (triclinic); NOT parametrized>

    You can also directly use ``u.atoms.convert_to.parmed()``

    Parameters
    ----------
    package: str
        The name of the package to convert to, e.g. ``"PARMED"``
    *args:
        Positional arguments passed to the converter
    **kwargs:
        Keyword arguments passed to the converter

    Returns
    -------
    output:
        An instance of the structure type from another package.

    Raises
    ------
    ValueError:
        No converter was found for the required package


    .. versionadded:: 1.0.0
    .. versionchanged:: 2.0.0
        Moved the ``convert_to`` method to its own class. The old API is still
        available and is now case-insensitive to package names, it also accepts
        positional and keyword arguments. Each converter function can also
        be accessed as a method with the name of the package in lowercase, i.e.
        `convert_to.parmed()`
    """
    _CONVERTERS: Type[ConverterBase] = {}

    def __init__(self, ag: "AtomGroup") -> None:
        """
        Parameters
        ----------
        ag : AtomGroup
            The AtomGroup to convert
        """
        self._ag = ag
        for lib, converter_cls in _CONVERTERS.items():
            method_name = lib.lower()
            # makes sure we always use the same instance of the converter
            # no matter which atomgroup instance called it
            try:
                converter = self._CONVERTERS[method_name]
            except KeyError:
                converter = converter_cls().convert
                # store in class attribute
                self._CONVERTERS[method_name] = converter
            # create partial function that passes ag to the converter
            convert = partial(converter, self._ag)
            # copy docstring and metadata to the partial function
            # note: it won't work with help()
            update_wrapper(convert, converter)
            setattr(self, method_name, convert)

    def __call__(self, package: str, *args: int, **kwargs: int) -> "ConverterBase":
        try:
            convert = getattr(self, package.lower())
        except AttributeError:
            raise ValueError(f"No {package!r} converter found. Available: "
                             f"{' '.join(self._CONVERTERS.keys())}") from None
        return convert(*args, **kwargs)
