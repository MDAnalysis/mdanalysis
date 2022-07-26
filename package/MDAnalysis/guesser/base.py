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
from .. import _GUESSERS

class GuesserMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            _GUESSERS[classdict['context'].upper()] = cls
        except KeyError:
            pass
        
     





class GuesserBase(metaclass=GuesserMeta):
    context = 'base' 
    _guess = {}

    # give a rank to each atrribute based on its dependency on other attributes
    # to be guessed, so that the attribute with lesser dependcy will be guessed
    # first
    _rank = {}

    def __init__(self, atoms, **kwargs):
        self._atoms = atoms
        self._kwargs = kwargs

    def is_guessed(self, guess):
        """check that the passed atrributes in the to_guess
        list can be guessed by the class

        Parameters
        ----------
        to_guess: list of atrributes sto be guessed then added to the universe
        Returns
        -------
        True or ValueError
        """
        for a in guess:
            if a.lower() not in self._guess:
                return False
        return True

    def guess_Attr(self, guess):
        """map the attribute to be guessed with the apporpiate guessing method

        Parameters
        ----------
        guess: an atrribute to be guessed then added to the universe
        Returns
        -------
        values: list of guessed values
        """
        return self._guess[guess]()

    def rank_attributes(self, attrs):
       """give a rank to each atrribute based on
          its dependency on other attributes to be guessed,
          so that the attribute with lesser dependcy will
          be guessed first

       Parameters
       ----------
       attrs: attributes list to be sorted

       Returns
       -------
       sorted attributes list
       """
       to_rank = {a: self._rank[a] for a in attrs}
       ranked_attrs = sorted(to_rank, key=to_rank.get)
       return ranked_attrs
   
    
   
    
def get_guesser(context, atoms, **kwargs):
    """get an appropiate guesser to the universe and pass
       the atomGroup of the universe to the guesser

    Parameters
    ----------
    atoms: AtomGroup of the universe
    context: string or Guesser class
    **kwargs: extra arguments are passed to the guesser.

    Returns
    -------
    Guesser class
    """
    if isinstance(context, GuesserBase):
        return context
    try:
        guesser = _GUESSERS[context.upper()](atoms, **kwargs)
    except KeyError:
        raise TypeError("Unidentified guesser type {0}".format(context))
    return guesser

