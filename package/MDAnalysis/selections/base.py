# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Base classes for the selection writers
======================================

Specialized SelectionWriters are derived from
:class:`SelectionWriter`. Override the :meth:`~SelectionWriter._write_head`,
:meth:`~SelectionWriter._translate`, and :meth:`~SelectionWriter._write_tail`
methods.

.. autoclass:: SelectionWriter
   :members: __init__, write, _translate, _write_head, _write_tail, comment

.. autofunction:: join

"""
from __future__ import absolute_import

from six.moves import range

import os.path

from ..lib import util

def join(seq, string="", func=None):
    """Create a list from sequence.

    *string* is appended to each element but the last.

    *func* is applied to every element before appending *string*.
    """
    if func is None:
        func = lambda x: x
    return [func(x) + string for x in seq[:-1]] + [func(seq[-1])]


class SelectionWriter(object):
    """Export a selection in MDAnalysis to a format usable in an external package.

    The :class:`SelectionWriter` writes a selection string to a file
    that can be used in another package such as `VMD`_, `PyMOL`_,
    `Gromacs`_ or `CHARMM`_. In this way, analysis and visualization
    can be done with the best or most convenient tools at hand.

    :class:`SelectionWriter` is a base class and child classes are
    derived with the appropriate customizations for the package file
    format.

    .. _VMD: http://www.ks.uiuc.edu/Research/vmd/
    .. _PyMol: http://www.pymol.org/
    .. _CHARMM:  http://www.charmm.org/
    .. _Gromacs: http://www.gromacs.org/

    .. versionchanged:: 0.11.0
       Can now also write to a :class:`~MDAnalysis.lib.util.NamedStream` instead
       of a normal file (using :class:`~MDAnalysis.lib.util.openany`).
    """
    #: Name of the format.
    format = None
    #: Extension of output files.
    ext = None
    #: Special character to continue a line across a newline.
    continuation = ''
    #: Comment format string; should contain '%s' or ``None`` for no comments.
    commentfmt = None
    default_numterms = 8

    def __init__(self, filename, mode="wa", numterms=None, preamble=None, **kwargs):
        """Set up for writing to *filename*.

        :Arguments:
           *filename*
               output file
           *mode*
               overwrite ("w") for every write, append ("a") to existing
               file, or overwrite an existing file and then append ("wa") ["wa"]
           *numterms*
               number of individual index numbers per line for output
               formats that write multiple entries in one line. If set
               to 0 or ``False`` then no special formatting is done  [8]
           *preamble*
               string that is written as a comment at the top of the file []
           *kwargs*
               use as defaults for :meth:`write`
        """
        self.filename = util.filename(filename, ext=self.ext)
        if not mode in ('a', 'w', 'wa'):
            raise ValueError("mode must be one of 'w', 'a', 'wa', not {0!r}".format(mode))
        self.mode = mode
        self._current_mode = mode[0]
        if numterms is None or numterms < 0:
            self.numterms = self.default_numterms
        elif numterms is False:
            self.numterms = 0
        else:
            self.numterms = numterms
        self.preamble = preamble
        self.otherargs = kwargs  # hack
        self.number = 0

        self.write_preamble()

    def comment(self, s):
        """Return string *s* interpolated into the comment format string.

        If no :attr:`SelectionWriter.commentfmt` is defined (None) then the
        empty string is returned because presumably there is no way to enter
        comments into the file.

        A newline is appended to non-empty strings.
        """
        if self.commentfmt is None:
            return ''
        return self.commentfmt % s + '\n'

    def write_preamble(self):
        """Write a header, depending on the file format."""
        if self.preamble is None:
            return
        with util.openany(self.filename, self._current_mode) as out:
            out.write(self.comment(self.preamble))
        self._current_mode = 'a'

    def write(self, selection, number=None, name=None, frame=None, mode=None):
        """Write selection to the output file.

        :Arguments:
           *selection*
               a :class:`MDAnalysis.core.AtomGroup.AtomGroup`
           *number*
               selection will be named "mdanalysis<number>"
               (``None`` auto increments between writes; useful
               when appending) [``None``]
           *name*
               selection will be named *name* (instead of numbered)
               [``None``]
           *frame*
               write selection of this frame (or the current one if
               ``None`` [``None``]
        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            try:
                frame = u.trajectory.ts.frame
            except AttributeError:
                frame = 1  # should catch cases when we are analyzing a single PDB (?)
        name = name or self.otherargs.get('name', None)
        if name is None:
            if number is None:
                self.number += 1
                number = self.number
            name = "mdanalysis{number:03d}".format(**vars())
        # build whole selection in one go (cleaner way to deal with
        # to deal with line breaks after self.numterms entries)
        # selection_list must contain entries to be joined with spaces or linebreaks
        selection_terms = self._translate(selection.atoms)
        step = self.numterms or len(selection.atoms)
        with util.openany(self.filename, self._current_mode) as out:
            self._write_head(out, name=name)
            for iatom in range(0, len(selection.atoms), step):
                line = selection_terms[iatom:iatom + step]
                out.write(" ".join(line))
                if len(line) == step and not iatom + step == len(selection.atoms):
                    out.write(' ' + self.continuation + '\n')
            out.write(' ')  # safe so that we don't have to put a space at the start of tail
            self._write_tail(out)
            out.write('\n')  # always terminate with newline

            if self.mode == 'wa':
                self._current_mode = 'a'  # switch after first write
            elif self.mode == 'w':
                self._current_mode = 'w'  # switch back after eg preamble

    def _translate(self, atoms, **kwargs):
        """Translate atoms into a list of native selection terms.

        - build list of ALL selection terms as if this was a single line, e.g.
          ``['index 12 |', 'index 22 |', 'index 33']``
        - only one term per atom!!
        - terms *must* be strings
        - something like::
             " ".join(terms)

          must work
        """
        raise NotImplementedError

    def _write_head(self, out, **kwargs):
        """Initial output to open file object *out*."""
        pass

    def _write_tail(self, out, **kwargs):
        """Last output to open file object *out*."""
        pass
