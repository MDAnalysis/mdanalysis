.. -*- coding: utf-8 -*-
.. _formats:

====================
Format overview
====================

MDAnalysis can read topology or coordinate information from a wide variety of file formats. The emphasis is on formats used in popular simulation packages. By default, MDAnalysis figures out formats by looking at the extension, unless the format is :ref:`explicitly specified <universe-loading>` with the ``format`` or ``topology_format`` keywords.

Below is :ref:`a table of formats in MDAnalysis <format-overview>`, and which information can be read from them. A topology file supplies the list of atoms in the system, their connectivity and possibly additional information such as B-factors, partial charges, etc. The details depend on the file format and not every topology file provides all (or even any) additional data.

.. important ::

    File formats are complicated and not always well defined. MDAnalysis tries to follow published standards but this can sometimes surprise users. It is *highly* recommended that you read the page for your data file format instead of assuming certain behaviour. If you encounter problems with a file format, please :ref:`get in touch with us <contributing>`.


As a minimum, all topology parsers will provide atom ``ids``, atom ``types``, ``masses``, ``resids``, ``resnums``, and ``segids``. They will also assign all Atoms to Residues and all Residues to Segments. For systems without residues and segments, this results in there being a single Residue and Segment to which all Atoms belong. See :ref:`topology-attributes` for more topology attributes.

Often when data is not provided by a file, it will be guessed based on other data in the file. In this scenario, MDAnalysis will issue a warning. See :ref:`guessing` for more information.

If a trajectory is loaded without time information, MDAnalysis will set a default timestep of 1.0 ps, where the first frame starts at 0.0 ps. In order to change these, :ref:`pass the following optional arguments to Universe <universe-kwargs>`:

    * ``dt``: the timestep
    * ``time_offset``: the starting time from which to calculate the time of each frame

.. _format-overview:

.. table:: Table of all supported formats in MDAnalysis

    .. include:: format_overview.txt

.. _topology-parsers:

.. include:: topology.rst

.. include:: coordinates.rst