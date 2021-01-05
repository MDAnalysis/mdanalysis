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

import numpy as np

from .universe import Universe
from .groups import AtomGroup
from ..lib import util
from ..coordinates.base import ProtoReader


class Ensemble(object):
    """
    Should be able to plug this into any AnalysisBase class.
    """

    def __init__(self, universes, select=None, labels=None, frames=None):
        # set universe info
        universes = util.asiterable(universes)
        try:
            self.universes = [u.universe for u in universes]
        except AttributeError:
            raise ValueError('universes must be a list of Universes '
                             'or AtomGroups')
        self.n_universes = len(universes)
        if frames is None:
            frame_ix = [np.arange(len(u.trajectory)) for u in self.universes]
            frames = []
            for ix in frame_ix:
                frames.extend(ix+len(frames))
        self.frames = np.asarray(frames, dtype=int)
        self.n_frames = len(self.frames)
        self.traj_frames = tuple([len(u.trajectory) for u in self.universes])
        self._frame_edges = np.r_[0, np.cumsum(self.traj_frames)]
        self._universe_frames = np.zeros((self.n_frames, 2), dtype=int)
        for i, f in enumerate(self.frames):
            n_u = self._get_universe_from_frame(f)
            self._universe_frames[i] = n_u, f-self._frame_edges[n_u]

        # set atom info
        if select is None:
            n = min(len(ag.atoms) for ag in universes)
            self._ags = [ag.atoms[:n] for ag in universes]
        elif isinstance(select, str):
            self._ags = [u.select_atoms(select) for u in universes]
        elif util.isiterable(select) and len(select) == self.n_universes:
            self._ags = [u.select_atoms(s) for u, s in zip(universes, select)]
        else:
            err = ('select must be None, a string, or an iterable of strings '
                   'with the same length as universes')
            raise ValueError(err)

        self._ag_indices = [ag.ix for ag in self._ags]
        self.n_atoms = len(self._ags[0])

        # label universes with filename or user-given name
        labels = util.asiterable(labels)
        if len(labels) == 1:
            labels = [labels[0]]*self.n_universes
        if len(labels) != self.n_universes:
            raise ValueError('labels must be None, a string, or an iterable '
                             'of strings with the same length as universes')
        self.labels = []
        self._universe_labels = {}
        _ctr = {}
        for i, (l, u) in enumerate(zip(labels, self.universes)):
            if l is None:
                try:
                    l = u.trajectory.filename
                except AttributeError:
                    l = 'Universe'

            if l in self._universe_labels:
                if l not in _ctr:
                    _ctr[l] = 1
                _ctr[l] += 1
                l += ' {}'.format(_ctr[l])
                
            self.labels.append(l)
            self._universe_labels[l] = i

        # pretend to be Universe
        self.trajectory = self
        self._ts_u = 0

    def __len__(self):
        return self.n_frames

    def __getitem__(self, i):
        if isinstance(i, (int, np.integer)):
            n_u, n_frame = self._get_relative_frame(i)
            u = self.universes[n_u]
            u.trajectory[n_frame]
            self._ts_u = n_u
        elif isinstance(i, str):
            self._ts_u = self._universe_labels[i]
        elif util.iterable(i) or isinstance(i, slice):
            return self._slice_universe(i)
        else:
            raise ValueError('You can only index or slice Ensembles with '
                             'integers, slices, arrays, or string labels')
        return self.ts

    def __iter__(self):
        for n_u, f in self._universe_frames:
            self._ts_u = n_u
            self.universes[n_u].trajectory[f]
            yield self.ts

    @property
    def atoms(self):
        return self._ags[self._ts_u]

    @property
    def ts(self):
        return self.universes[self._ts_u].trajectory.ts

    @property
    def universe(self):
        return self.universes[self._ts_u]

    @property
    def positions(self):
        return self.atoms.positions

    check_slice_indices = ProtoReader.check_slice_indices

    def iterate_over_atomgroups(self, start=None, stop=None, step=None,
                                frames=None):
        if frames is None:
            if start is None:
                start = 0
            if stop is None:
                stop = self.n_frames-1
            stop = max(self.n_frames-1, stop)
            if step is None:
                step = 1
            frames = range(start, stop, step)

        for i in frames:
            n_u, f = self._universe_frames[i]
            self._ts_u = n_u
            self.universe.trajectory[f]
            yield self._ags[n_u]

    def _get_relative_frame(self, i):
        if not isinstance(i, (int, np.integer)):
            raise IndexError('only integers are valid indices. '
                             'Given: {}, type {}'.format(i, type(i).__name__))

        err = ('{} is out of range for ensemble '
               'with length {}').format(i, self.n_frames)

        return self._universe_frames[i]

    def _get_universe_from_frame(self, i):
        n_u = np.searchsorted(self._frame_edges, i, side='right')
        return n_u-1

    def _slice_universe(self, sl):
        frames = self.frames[sl]
        return type(self)(self._ags, labels=self.labels, frames=frames)

    def select_atoms(self, *sel):
        sel = [s for s in sel if s is not None]
        ags = self._ags
        if sel:
            ags = [ag.select_atoms(*sel) for ag in ags]
        return type(self)(ags, select=None, labels=self.labels)
