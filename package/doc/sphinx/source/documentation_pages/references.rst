.. -*- coding: utf-8 -*-
.. note: make sure that no lines accidentaly start with a single character
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
:mod:`~MDAnalysis.core.qcprot` module please also cite [Theobald2005b]_
and :cite:p:`a-Liu2010`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Liu2010

.. [Theobald2005b] Douglas L. Theobald. Rapid calculation of RMSD using a
   quaternion-based characteristic polynomial. *Acta Crystallographica A*
   **61** (2005), 478-480.

If you use the helix analysis algorithm HELANAL_ in
:mod:`MDAnalysis.analysis.helix_analysis` please cite :cite:p:`a-Bansal2000`.

.. bibliography::
   :filter: False
   :style: MDA
   :keyprefix: a-
   :labelprefix: ᵃ

   Bansal2000

.. _HELANAL: http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.html

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
:mod:`MDAnalysis.analysis.encore` please cite [Tiberti2015b]_.

.. [Tiberti2015b] M. Tiberti, E. Papaleo, T. Bengtsen, W. Boomsma,
   and K. Lindorff-Larsen. ENCORE: Software for quantitative ensemble
   comparison. *PLoS Comput Biol*, **11** (2015), e1004415.  doi:
   `10.1371/journal.pcbi.1004415`_

.. _`10.1371/journal.pcbi.1004415`: http://doi.org/10.1371/journal.pcbi.1004415

If you use the streamline visualization in
:mod:`MDAnalysis.visualization.streamlines` and
:mod:`MDAnalysis.visualization.streamlines_3D` please cite [Chavent2014b]_.

.. [Chavent2014b] Chavent, M., Reddy, T., Dahl, C.E., Goose, J., Jobard, B.,
   and Sansom, M.S.P. Methodologies for the analysis of instantaneous lipid
   diffusion in MD simulations of large membrane systems.  *Faraday
   Discussions* **169** (2014), 455–475. doi: `10.1039/c3fd00145h`_

.. _`10.1039/c3fd00145h`: https://doi.org/10.1039/c3fd00145h

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
:func:`~MDAnalysis.analysis.pca.cumulative_overlap` please cite [Yang2008]_ .

.. [Yang2008] Yang, L., Song, G., Carriquiry, A. & Jernigan, R. L.
   Close Correspondence between the Motions from Principal Component Analysis of Multiple HIV-1 Protease Structures and Elastic Network Modes.
   *Structure* **16**, 321–330 (2008). doi: `10.1016/j.str.2007.12.011`_

.. _`10.1016/j.str.2007.12.011`: https://dx.doi.org/10.1016/j.str.2007.12.011

If you use the Mean Squared Displacement analysis code in
:mod:`MDAnalysis.analysis.msd` please cite [Calandri2011]_ and [Buyl2018]_.

.. [Calandri2011] Calandrini, V., Pellegrini, E., Calligari, P., Hinsen, K., Kneller, G. R.
   NMoldyn-Interfacing Spectroscopic Experiments, Molecular Dynamics Simulations and Models for Time Correlation Functions.
   *Collect. SFN*, **12**, 201–232 (2011). doi: `10.1051/sfn/201112010`_

.. _`10.1051/sfn/201112010`: https://doi.org/10.1051/sfn/201112010

.. [Buyl2018] Buyl, P. tidynamics: A tiny package to compute the dynamics of stochastic and molecular simulations. Journal of Open Source Software,
   3(28), 877 (2018). doi: `10.21105/joss.00877`_

.. _`10.21105/joss.00877`: https://doi.org/10.21105/joss.00877

If you calculate shape parameters using
:meth:`~MDAnalysis.core.group.AtomGroup.shape_parameter`,
:meth:`~MDAnalysis.core.group.ResidueGroup.shape_parameter`,
:meth:`~MDAnalysis.core.group.SegmentGroup.shape_parameter`
please cite [Dima2004a]_.

.. [Dima2004a] Dima, R. I., & Thirumalai, D. (2004). Asymmetry
   in the shapes of folded and denatured states of
   proteins. *J Phys Chem B*, 108(21),
   6564-6570. doi:`10.1021/jp037128y
   <https://doi.org/10.1021/jp037128y>`_

If you calculate asphericities using
:meth:`~MDAnalysis.core.group.AtomGroup.asphericity`,
:meth:`~MDAnalysis.core.group.ResidueGroup.asphericity`,
:meth:`~MDAnalysis.core.group.SegmentGroup.asphericity`
please cite [Dima2004b]_.

.. [Dima2004b] Dima, R. I., & Thirumalai, D. (2004). Asymmetry
   in the shapes of folded and denatured states of
   proteins. *J Phys Chem B*, 108(21),
   6564-6570. doi:`10.1021/jp037128y
   <https://doi.org/10.1021/jp037128y>`_

If you use use the dielectric analysis code in
:class:`~MDAnalysis.analysis.dielectric.DielectricConstant` please cite [Neumann1983]_.

.. [Neumann1983] Neumann, M. (1983). Dipole
   Moment Fluctuation Formulas in Computer Simulations of Polar Systems.
   *Molecular Physics* **50**, no. 4, 841–858. doi: `10.1080/00268978300102721`_

.. _`10.1080/00268978300102721`: http://doi.org/10.1080/00268978300102721

If you use H5MD files using
:mod:`MDAnalysis.coordinates.H5MD.py`, please cite [Buyl2013]_ and
[Jakupovic2021]_.

.. [Buyl2013] Buyl P., Colberg P., and Höfling F.(2013).
   H5MD: A structured, efficient, and portable file format for molecular data.
   *Computer Physics Communications*, 185. doi:`10.1016/j.cpc.2014.01.018.
   <https://doi.org/10.1016/j.cpc.2014.01.018>`_

.. [Jakupovic2021] Jakupovic E. and Beckstein O., MPI-parallel Molecular
   Dynamics Trajectory Analysis with the H5MD Format in the MDAnalysis
   Python Package, in *Proceedings of the 20th Python in Science Conference*,
   (Meghann Agarwal, Chris Calloway, Dillon Niederhut, and David Shupe, eds.),
   pp. 18 – 26, 2021. doi:`10.25080/majora-1b6fd038-005.
   <https://www.doi.org/10.25080/majora-1b6fd038-005>`_


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

   DUECREDIT-ENABLE=yes python yourscript.py

Once the citations have been extracted (to a hidden file in the
current directory), you can use the :program:`duecredit` program to
export them to different formats. For example, one can display them in
BibTeX format, using:

.. code-block:: bash

   duecredit summary --format=bibtex


**Please cite your use of MDAnalysis and the packages and algorithms
that it uses. Thanks!**


.. _duecredit: https://github.com/duecredit/duecredit
