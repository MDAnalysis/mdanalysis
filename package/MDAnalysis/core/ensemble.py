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


class AtomsWrapper:
    def __init__(self, ensemble):
        self.ensemble = ensemble
        self.universe = ensemble

    def __len__(self):
        return len(self.ensemble._atoms)

    def __getattr__(self, attr):
        try:      
            return getattr(self.ensemble._atoms, attr)
        except AttributeError:
            return getattr(self.ensemble, attr)

    def select_atoms(self, *args, **kwargs):
        return self.ensemble.select_atoms(*args, **kwargs).atoms

    @property
    def positions(self):
        return self.ensemble._atoms.positions
    
    @positions.setter
    def positions(self, value):
        self.ensemble._atoms.positions = value


class Ensemble(object):
    """
    Should be able to plug this into any AnalysisBase class.
    """

    def __init__(self, universes, select=None, labels=None, frames=None,
                 _ts_u=0, links=None):
        self.ensemble = self
        self.universe = self
        # set universe info
        universes = util.asiterable(universes)
        self.atoms = AtomsWrapper(self)
        try:
            self.universes = [u.universe for u in universes]
        except AttributeError:
            raise ValueError('universes must be a list of Universes '
                             'or AtomGroups')
        self.n_universes = len(universes)
        self._set_frames(frames)

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
        self._ags = tuple(self._ags)
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
                    pass
            if l is None:
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
        self._ts_u = _ts_u
        
        if links is None:
            links = set()
        self.links = links
        self.links.add(self)

    def __eq__(self, other):
        return self._ags == other._ags

    def __hash__(self):
        return hash(id(self))

    def __getattr__(self, attr):
        try:
            return getattr(self.atoms, attr)
        except AttributeError:
            try:
                return getattr(self.universe, attr)
            except AttributeError:
                pass
        raise AttributeError(f"{type(self).__name__} does not have "
                             f"attribute {attr}")

    def __repr__(self):
        return (f"<{type(self).__name__} with {self.n_atoms} atoms, "
                f"{len(self)} frames, frame={self.frame} at {hex(id(self))}>")

    def __len__(self):
        return self.n_frames

    def __getitem__(self, i):
        if isinstance(i, (int, np.integer)):
            n_u, n_frame = self._get_relative_frame(i)
            self._set_frame(ts_u=n_u, n_frame=n_frame)
        elif isinstance(i, str):
            self._set_frame(ts_u=self._universe_labels[i])
        elif util.iterable(i) or isinstance(i, slice):
            return self._slice_universe(i)
        else:
            raise ValueError('You can only index or slice Ensembles with '
                             'integers, slices, arrays, or string labels')
        return self.ts

    def _set_frame(self, ts_u=None, n_frame=None):
        if ts_u is not None:
            for ens in self.links:
                ens._ts_u = ts_u
        
        if n_frame is not None:
            self._universe.trajectory[n_frame]
            for ens in self.links:
                ens._universe.trajectory[n_frame]

    @property
    def frame(self):
        ix = self._universe_frames
        u = self.universes[self._ts_u]
        u_ix = np.where(ix[:, 0] == self._ts_u)[0]
        frame = u.trajectory.frame
        fr_ix = np.where(ix[:, 1] == frame)[0]
        try:
            return u_ix[np.in1d(u_ix, fr_ix)][0]
        except IndexError:
            return None

    @property
    def positions(self):
        return self._atoms.positions

    @positions.setter
    def positions(self, value):
        self._atoms.positions = value


    def __iter__(self):
        for n_u, f in self._universe_frames:
            self._set_frame(ts_u=n_u, n_frame=f)
            yield self.ts
        self._set_frame(ts_u=0)

    def timeseries(self, atomgroup=None, start=None, stop=None, step=None,
                   frames=None, order="afc"):
        n_atoms = self.n_atoms
        frames = np.array(self._prepare_frames(start=start, stop=stop,
                                               step=step, frames=frames))
        n_frames = len(frames)
        fac = np.zeros((n_frames, n_atoms, 3))
        for i, frame in enumerate(frames):
            self[frame]
            np.copyto(fac[i], self.atoms.positions)
        if order == "afc":
            np.swapaxes(fac, 0, 1)
        elif order == "caf":
            np.swapaxes(fac, 2, 0)
        elif order == "fca":
            np.swapaxes(fac, 2, 1)
        elif order == "cfa":
            fac = fac.transpose((2, 0, 1))
        elif order == "acf":
            fac = fac.transpose((1, 2, 0))
        return fac
        

    @property
    def filename(self):
        return self.universe.trajectory.filename

    @property
    def _atoms(self):
        return self._ags[self._ts_u]

    @property
    def ts(self):
        return self.universes[self._ts_u].trajectory.ts

    @property
    def _universe(self):
        return self.universes[self._ts_u]

    check_slice_indices = ProtoReader.check_slice_indices

    def _set_frames(self, frames=None):
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

    def _prepare_frames(self, start=None, stop=None, step=None,
                        frames=None):
        if frames is None:
            if start is None:
                start = 0
            if stop is None:
                stop = self.n_frames
            stop = max(self.n_frames-1, stop)
            if step is None:
                step = 1
            frames = range(start, stop, step)
        return frames

    def _prepare_frames_by_universe(self, start=None, stop=None, step=None,
                                    frames=None):
        frames = np.array(self._prepare_frames(start=start, stop=stop,
                                               step=step, frames=frames))
        u_frames = self._universe_frames[frames]
        splix = np.where(np.ediff1d(u_frames[:, 0]))[0] + 1
        return np.split(u_frames[:, 1], splix)

        # bins = np.searchsorted(self._frame_edges[1:], frames)
        # for i in range(self.n_universes):
        #     yield frames[bins == i]

    def iterate_over_atomgroups(self, start=None, stop=None, step=None,
                                frames=None):
        frames = self._prepare_frames(start=start, stop=stop, step=step,
                                      frames=frames)

        for i in frames:
            n_u, f = self._universe_frames[i]
            self._set_frame(n_u, f)
            yield self._ags[n_u]

    def iterate_over_universes(self, start=None, stop=None, step=None,
                               frames=None):
        frames_ = self._prepare_frames_by_universe(start=start, stop=stop,
                                                   step=step, frames=frames)
        for ag, label, fr in zip(self._ags, self.labels, frames_):
            yield type(self)([ag], labels=[label], frames=fr, _ts_u=self._ts_u,
                             links=self.links)


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
        return type(self)(self._ags, labels=self.labels, frames=frames,
                          _ts_u=self._ts_u, links=self.links)

    def select_atoms(self, *sel):
        sel = [s for s in sel if s is not None]
        ags = self._ags
        if sel:
            ags = [ag.select_atoms(*sel) for ag in ags]
        return type(self)(ags, select=None, labels=self.labels,
                          _ts_u=self._ts_u, links=self.links)

    def split_array(self, arr):
        for i in range(self.n_universes):
            yield arr[self._frame_edges[i]:self._frame_edges[i+1]]
    
    def transfer_to_memory(self, start=None, stop=None, step=None,
                           frames=None):

        frames_ = self._prepare_frames_by_universe(start=start, stop=stop,
                                                   step=step, frames=frames)
        frames = np.array(list(frames_))
        for u, fr in zip(self.universes, frames):
            u.transfer_to_memory(frames=fr)
        self._set_frames()