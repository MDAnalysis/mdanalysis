
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
"""
import io


class FileIOPicklable(io.FileIO):
    """Stream for read a file

    Picklable FiloIO class that only support read mode.

    """
    def __getstate__(self):
        return self.name, self.tell()

    def __setstate__(self, args):
        name = args[0]
        super().__init__(name)
        self.seek(args[1])


class BufferIOPicklable(io.BufferedReader):
    """A picklable buffer for a readable FilIO object

    Wrap raw FileIOPicklable inside

    """
    def __init__(self, raw):
        super().__init__(raw)
        self.raw_class = raw.__class__

    def __getstate__(self):
        return self.raw_class, self.name, self.tell()

    def __setstate__(self, args):
        raw_class = args[0]
        name = args[1]
        raw = raw_class(name)
        super().__init__(raw)
        self.seek(args[2])


class TextIOPicklable(io.TextIOWrapper):
    """Character and line based layer over a picklable FileIO based object.

    Example
    -------
        file = FileIOPicklable('filename')
        text_wrapped = TextIOPicklable(file)
    """
    def __init__(self, raw):
        super().__init__(raw)
        self.raw_class = raw.__class__

    def __getstate__(self):
        return self.raw_class, self.name, self.tell()

    def __setstate__(self, args):
        raw_class = args[0]
        name = args[1]
        # raw_class is used for further expansion this functionality to
        # GZip files, which also requires a text wrapper.
        raw = raw_class(name)
        super().__init__(raw)
        self.seek(args[2])


def pickle_open(name, mode='rt'):
    """Open file and return a stream with pickle function implemented.

    Parameters
    ----------
    name : str;
        a filename given a text or byte string.
    mode: {'r', 'rt', 'rb'} (optional)
        'r':  open for reading in text mode;
        'rt': read in text mode (default);
        'rb': read in binary mode;
        raise ValueError with other modes.

    Returns
    -------
    stream-like object

    See Also
    --------
    :func:`anyopen`

    Warning
    -------
    Should be only used with read mode.
    """
    if mode not in {'r', 'rt', 'rb'}:
        raise ValueError("Only read mode ('r', 'rt', 'rb') files can be pickled.")
    raw = FileIOPicklable(name)
    if mode == 'rb':
        return BufferIOPicklable(raw)
    elif mode in {'r', 'rt'}:
        return TextIOPicklable(raw)
