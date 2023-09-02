.. -*- coding: utf-8 -*-
.. note: make sure that no lines accidentally start with a single character
..       followed by a period: reST interprets it as an enumerated list and
..       messes up the formatting

.. The references are accessible globally; you can cite these papers anywhere
.. in the docs.

.. _references:

************
 References
************

MDAnalysis and the included algorithms are scientific software that
are described in academic publications. **Please cite these papers when you use
MDAnalysis in published work.**

It is possible to :ref:`automatically generate a list of references
<citations-using-duecredit>` for any program that uses
MDAnalysis. This list (in common reference manager formats) contains
the citations associated with the specific algorithms and libraries
that were used in the program.


Citations for the whole MDAnalysis library
==========================================

When using MDAnalysis in published work, please cite
:cite:p:`Michaud-Agrawal2011` and :cite:p:`Gowers2016`.

(We are currently asking you to cite *both* papers if at all possible
because the 2016 paper describes many updates to the original 2011
paper and neither paper on its own provides a comprehensive
description of the library. We will publish a complete self-contained
paper with the upcoming 1.0 release of MDAnalysis, which will then
supersede these two citations.)

.. bibliography::
   :filter: False
   :style: MDA

   Michaud-Agrawal2011
   Gowers2016


.. _references-components:

Citations for included algorithms and modules
=============================================

If you use the RMSD calculation (:mod:`MDAnalysis.analysis.rms`) or alignment
code (:mod:`MDAnalysis.analysis.align`) that uses the
:mod:`~MDAnalysis.lib.qcprot` module please also cite :cite:p:`a-Theobald2005`
and :cite:p:`a-Liu2010`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Liu2010
   Theobald2005

If you use the HELANAL_ algorithm in
:mod:`MDAnalysis.analysis.helix_analysis` please cite :cite:p:`a-Bansal2000,a-Sugeta1967`.

.. _HELANAL: https://web.archive.org/web/20100818185943/http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.html

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Bansal2000
   Sugeta1967

If you use the GNM trajectory analysis code in
:mod:`MDAnalysis.analysis.gnm` please cite :cite:p:`a-Hall2007`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Hall2007

If you use the water analysis code in
:mod:`MDAnalysis.analysis.waterdynamics` please cite :cite:p:`a-ArayaSecchi2014`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   ArayaSecchi2014


If you use the Path Similarity Analysis (PSA) code in
:mod:`MDAnalysis.analysis.psa` please :cite:p:`a-Seyler2015`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Seyler2015

If you use the implementation of the ENCORE ensemble analysis in
:mod:`MDAnalysis.analysis.encore` please cite :cite:p:`a-Tiberti2015`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Tiberti2015

If you use the streamline visualization in
:mod:`MDAnalysis.visualization.streamlines` and
:mod:`MDAnalysis.visualization.streamlines_3D` please cite :cite:p:`a-Chavent2014`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Chavent2014

If you use the hydrogen bond analysis code in
:mod:`MDAnalysis.analysis.hydrogenbonds.hbond_analysis` please cite :cite:p:`a-Smith2019`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Smith2019

If you use :meth:`~MDAnalysis.analysis.pca.PCA.rmsip` or
:func:`~MDAnalysis.analysis.pca.rmsip` please cite :cite:p:`a-Amadei1999` and
:cite:p:`a-Leo-Macias2005`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Amadei1999
   Leo-Macias2005

If you use :meth:`~MDAnalysis.analysis.pca.PCA.cumulative_overlap` or
:func:`~MDAnalysis.analysis.pca.cumulative_overlap` please cite
:cite:p:`a-Yang2008`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Yang2008

If you use the Mean Squared Displacement analysis code in
:mod:`MDAnalysis.analysis.msd` please cite :cite:p:`a-Calandrini2011` and
:cite:p:`a-Buyl2018`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Calandrini2011
   Buyl2018

If you calculate shape parameters using
:meth:`~MDAnalysis.core.group.AtomGroup.shape_parameter`,
:meth:`~MDAnalysis.core.group.ResidueGroup.shape_parameter`,
:meth:`~MDAnalysis.core.group.SegmentGroup.shape_parameter`,
or if you calculate asphericities using
:meth:`~MDAnalysis.core.group.AtomGroup.asphericity`,
:meth:`~MDAnalysis.core.group.ResidueGroup.asphericity`,
:meth:`~MDAnalysis.core.group.SegmentGroup.asphericity`,
please cite :cite:p:`a-Dima2004`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Dima2004

If you use use the dielectric analysis code in
:class:`~MDAnalysis.analysis.dielectric.DielectricConstant` please cite
:cite:p:`a-Neumann1983`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Neumann1983

If you use H5MD files with
:mod:`MDAnalysis.coordinates.H5MD`, please cite :cite:p:`a-deBuyl2014` and
:cite:p:`a-Jakupovic2021`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   deBuyl2014
   Jakupovic2021

.. _citations-using-duecredit:


Citations using Duecredit
=========================

Citations can be automatically generated using duecredit_, depending on the
packages used. Duecredit is easy to install via ``pip``. Simply type:

.. code-block:: bash

   pip install duecredit

duecredit_ will remain an optional dependency, i.e. any code using
MDAnalysis will work correctly even without duecredit installed.

A list of citations for ``yourscript.py`` can be obtained using simple
commands.

.. code-block:: bash

   cd /path/to/yourmodule
   python -m duecredit yourscript.py

or set the environment variable :envvar:`DUECREDIT_ENABLE`

.. code-block:: bash

   DUECREDIT_ENABLE=yes python yourscript.py

Once the citations have been extracted (to a hidden file in the
current directory), you can use the :program:`duecredit` program to
export them to different formats. For example, one can display them in
BibTeX format, using:

.. code-block:: bash

   duecredit summary --format=bibtex


**Please cite your use of MDAnalysis and the packages and algorithms
that it uses. Thanks!**


.. _duecredit: https://github.com/duecredit/duecredit
