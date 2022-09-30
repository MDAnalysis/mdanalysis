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

Derive context-specific guesser classes from the base class in this module. 

Classes
-------

.. autoclass:: GuesserBase
   :members:
   :inherited-members:

"""
from .. import _GUESSERS
import numpy as np
from ..  import _TOPOLOGY_ATTRS
import logging
from typing import Dict

logger = logging.getLogger("MDAnalysis.guesser.base")


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

    .. versionadded:: 2.4.0
    """
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
       
        _GUESSERS[classdict['context'].upper()] = cls
     

class GuesserBase(metaclass=_GuesserMeta):
    """Base class for context-specific guessers to inherit from

    Parameters
    ----------
    universe : Universe
        Universe
    universe : Universe, optional
        Supply a Universe to the Parser. This then become the source of atom attributes
        to be used in guessing processes. (this is relevant to how the universe's guess_topologyAttributes API works. See :ref:`guess_TopologyAttributes <guess_TopologyAttributes>`).
    **kwargs: to pass additional data to the guesser that can be used with different methos.
    
    .. versionadded:: 2.4.0

    """
    context = 'base'
    _guesser_methods : Dict = {}

    def __init__(self, universe=None, **kwargs):
        self._universe = universe
        self._kwargs = kwargs

    @property
    def universe(self):
        return self._universe

    def update_kwargs(self, **kwargs):
        self._kwargs.update(kwargs)

    def are_guessable(self, guess):
        """check if the passed atrributes can be guessed by the guesser class

        Parameters
        ----------
        guess: list
            Atrributes to be guessed then added to the universe

        Returns
        -------
        Boolean value
        """
        if guess:
            for a in guess:
                if a.lower() not in self._guesser_methods:
                    return False
            return True
        return False
    

    def guess_Attr(self, attr_to_guess, force_guess=False):
        """map the attribute to be guessed with the apporpiate guessing method

        Parameters
        ----------
        attr_to_guess: list
            an atrribute to be guessed then to be added to the universe
        force_guess: boolean 
            To indicate wether to only partialy guess the empty values of the attribute or to overwrite all existing values by guessed one

        Returns
        -------
        list of guessed values
        """
        
        # check if the topology already has the attribute to partially guess it
        if hasattr(self._universe.atoms, attr_to_guess) and not force_guess:
            attr = np.array(getattr(self._universe.atoms, attr_to_guess, None))

            emptyAttr = []
            try:
                top_attr = _TOPOLOGY_ATTRS[attr_to_guess]
                for a in attr:
                   emptyAttr.append(top_attr.is_value_missing(a))
            except:
                pass

            if True in emptyAttr:
                parentAttr =[]
                # check if the parent attribute from which we can begin guessing exists in the universe
                for p in self._guesser_methods[attr_to_guess]['parent']: 
                    if hasattr(self._universe.atoms, p):
                        parentAttr = np.array(getattr(self._universe.atoms, p))
                        break

                # if didn't find any parent attribute we try to guess one parent first
                if len(parentAttr) == 0:
                    for p in self._guesser_methods[attr_to_guess]['parent']:
                        try:
                            self._universe.guess_TopologyAttributes(context=self.context, to_guess=[p])
                            parentAttr = np.array(getattr(self._universe.atoms, p))
                            break
                        except ValueError:
                                pass
    
                if len(parentAttr) > 0:
                    attr[emptyAttr] = self._guesser_methods[attr_to_guess]['guesser'](parentAttr[emptyAttr])

                else:
                    raise ValueError(f'Failed to find/guess a parent attribute to guess {attr_to_guess} from')

                logger.info(
                    f'attribute {attr_to_guess} has been guessed successfully.')
                return attr

            else:
                logger.info(
                    f'There is no empty {attr_to_guess} values. Guesser did not guess any new values for {attr_to_guess} attribute')
                return None
        else:
            return np.array(self._guesser_methods[attr_to_guess]['guesser']())


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

    .. versionadded:: 2.4.0

    """
    if isinstance(context, GuesserBase):
        context._universe = u
        context.update_kwargs(**kwargs)
        return context
    try:
        guesser = _GUESSERS[context.upper()](u, **kwargs)
    except KeyError:
        raise KeyError("Unidentified guesser type {0}".format(context))
    return guesser
