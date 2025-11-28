# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
Base classes for the selection writers
======================================

Specialized SelectionWriters are derived from
:class:`SelectionWriterBase`. Override the :meth:`~SelectionWriterBase._write_head`,
:meth:`~SelectionWriterBase._translate`, and :meth:`~SelectionWriterBase._write_tail`
methods.

.. autoclass:: SelectionWriterBase
   :members: __init__, write, _translate, _write_head, _write_tail, comment

.. autofunction:: join

"""
import os.path

from ..lib import util
from . import _SELECTION_WRITERS


def join(seq, string="", func=None):
    """Create a list from sequence.

    *string* is appended to each element but the last.

    *func* is applied to every element before appending *string*.
    """
    if func is None:
        func = lambda x: x
    return [func(x) + string for x in seq[:-1]] + [func(seq[-1])]


class _Selectionmeta(type):
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = util.asiterable(classdict["format"])
        except KeyError:
            pass
        else:
            for f in fmt:
                if f is None:
                    continue
                f = f.upper()
                _SELECTION_WRITERS[f] = cls


class SelectionWriterBase(metaclass=_Selectionmeta):
    """Export a selection in MDAnalysis to a format usable in an external package.

    The :class:`SelectionWriterBase` writes a selection string to a file
    that can be used in another package such as `VMD`_, `PyMOL`_,
    `Gromacs`_ or `CHARMM`_. In this way, analysis and visualization
    can be done with the best or most convenient tools at hand.

    :class:`SelectionWriterBase` is a base class and child classes are
    derived with the appropriate customizations for the package file
    format.

    .. _VMD: http://www.ks.uiuc.edu/Research/vmd/
    .. _PyMol: http://www.pymol.org/
    .. _CHARMM:  http://www.charmm.org/
    .. _Gromacs: http://www.gromacs.org/

    .. versionchanged:: 0.11.0
       Can now also write to a :class:`~MDAnalysis.lib.util.NamedStream` instead
       of a normal file (using :class:`~MDAnalysis.lib.util.openany`).

    .. versionchanged:: 0.16.0
       Remove the `wa` mode. The file is now open when the instance is created
       and closed with the :meth:`close` method or when exiting the `with`
       statement.
    """

    #: Name of the format.
    format = None
    #: Extension of output files.
    ext = None
    #: Special character to continue a line across a newline.
    continuation = ""
    #: Comment format string; should contain '%s' or ``None`` for no comments.
    commentfmt = None
    default_numterms = 8

    def __init__(
        self, filename, mode="w", numterms=None, preamble=None, **kwargs
    ):
        """Set up for writing to *filename*.

        Parameters
        ----------
        filename:
            output file
        mode:
            create a new file ("w"), or append ("a") to existing file ["w"]
        numterms:
            number of individual index numbers per line for output
            formats that write multiple entries in one line. If set
            to 0 or ``False`` then no special formatting is done  [8]
        preamble:
            string that is written as a comment at the top of the file []
        kwargs:
            use as defaults for :meth:`write`
        """
        self.filename = util.filename(filename, ext=self.ext)
        if not mode in ("a", "w"):
            raise ValueError(
                "mode must be one of 'w', 'a', not {0!r}".format(mode)
            )
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

        self._outfile = util.anyopen(self.filename, mode=self._current_mode)

        self.write_preamble()

    def comment(self, s):
        """Return string *s* interpolated into the comment format string.

        If no :attr:`SelectionWriterBase.commentfmt` is defined (None) then the
        empty string is returned because presumably there is no way to enter
        comments into the file.

        A newline is appended to non-empty strings.
        """
        if self.commentfmt is None:
            return ""
        return self.commentfmt % s + "\n"

    def write_preamble(self):
        """Write a header, depending on the file format."""
        if self.preamble is None:
            return
        self._outfile.write(self.comment(self.preamble))

    def write(self, selection, number=None, name=None, frame=None, mode=None):
        """Write selection to the output file.

        Parameters
        ----------
        selection:
            a :class:`MDAnalysis.core.groups.AtomGroup`
        number:
            selection will be named "mdanalysis<number>"
            (``None`` auto increments between writes; useful
            when appending) [``None``]
        name:
            selection will be named *name* (instead of numbered) [``None``]
        frame:
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
        name = name or self.otherargs.get("name", None)
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

        out = self._outfile
        self._write_head(out, name=name)
        for iatom in range(0, len(selection.atoms), step):
            line = selection_terms[iatom : iatom + step]
            out.write(" ".join(line))
            if len(line) == step and not iatom + step == len(selection.atoms):
                out.write(" " + self.continuation + "\n")
        out.write(
            " "
        )  # safe so that we don't have to put a space at the start of tail
        self._write_tail(out)
        out.write("\n")  # always terminate with newline

    def close(self):
        """Close the file

        .. versionadded:: 0.16.0
        """
        self._outfile.close()

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
        pass  # pylint: disable=unnecessary-pass

    def _write_tail(self, out, **kwargs):
        """Last output to open file object *out*."""
        pass  # pylint: disable=unnecessary-pass

    # Context manager support to match Coordinate writers
    # all file handles use a with block in their write method, so these do nothing special
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
