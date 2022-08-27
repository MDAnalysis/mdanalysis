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
Base guesser classes --- :mod:`MDAnalysis.guesser.base`
================================================================

Derive topology reader classes from the base class in this module. All
topology readers raise :exc:`IOError` upon failing to read a topology
file and :exc:`ValueError` upon failing to make sense of the read data.

Classes
-------

.. autoclass:: GuesserBase
   :members:
   :inherited-members:

"""
from .. import _GUESSERS


class _GuesserMeta(type):
    """Internal: guesser classes registration

    When classes which inherit from GuesserBase are *defined*
    this metaclass makes it known to MDAnalysis.  'context'
    attribute  are read:
     - `context` defines the context of the guesser class for example: forcefield specific context
     as MartiniGuesser and file specific context as PDBGuesser.

    Eg::

      class FooGuesser(GuesserBase):
          format = 'foo'

  """
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
       
        _GUESSERS[classdict['context'].upper()] = cls
     

class GuesserBase(metaclass=_GuesserMeta):
    """Base class for context-aware guessers

    Parameters
    ----------
    universe : Universe, optional
        Supply a Universe to the guesser.  This then become the source of atom attributes
        to be used in guessing processes (more relevant with the universe guess_topologyAttributes API).
    **kwargs: to pass additional data to the guesser that can be used with different methos. 
    Eg van der Waals radii that is used with bond guessing methods
 
   """
    context = 'base'
    _guess = {}

    def __init__(self, universe=None, **kwargs):
        self._universe = universe
        self._kwargs = kwargs

    @property
    def universe(self):
        return self._universe
    
    @universe.setter
    def set_universe(self, u):
        self._universe = u

    def update_kwargs(self, **kwargs):
        self._kwargs.update(kwargs)

    def is_guessable(self, guess):
        """check that the passed atrributes in the to_guess
        list can be guessed by the class

        Parameters
        ----------
        guess: list of atrributes sto be guessed then added to the universe
        Returns
        -------
        Boolean value
        """
        for a in guess:
            if a.lower() not in self._guess:
                return False
        return True

    def guess_Attr(self, guess, **kwargs):
        """map the attribute to be guessed with the apporpiate guessing method

        Parameters
        ----------
        guess: an atrribute to be guessed then to be added to the universe
        
        Returns
        -------
        values: list of guessed values
        """
        self._kwargs.update(kwargs)
        return self._guess[guess]()


def get_guesser(context, u=None, **kwargs):
    """get an appropiate guesser to the universe and pass
       the atomGroup of the universe to the guesser

    Parameters
    ----------
    u: Universe
    to be passed to the guesser 

    context: string or Guesser class
    **kwargs: extra arguments are passed to the guesser.

    Returns
    -------
    Guesser class
    
    Raises
    ------
    * :exc:`KeyError` upon failing to return a guesser class

    """
    if isinstance(context, GuesserBase):
        context.set_universe = u
        context.update_kwargs(**kwargs)
        return context
    try:
        guesser = _GUESSERS[context.upper()](u, **kwargs)
    except KeyError:
        raise KeyError("Unidentified guesser type {0}".format(context))
    return guesser
