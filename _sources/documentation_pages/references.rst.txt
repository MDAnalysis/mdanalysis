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


Citations for the whole MDAnalysis library
==========================================

When using MDAnalysis in published work, please cite
[Michaud-Agrawal2011]_ and [Gowers2016]_.

.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf,
   and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics
   Simulations. *J. Comput. Chem.* **32** (2011),
   2319–2327. doi:`10.1002/jcc.21787`_

.. [Gowers2016] R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N.
   Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney,
   and O. Beckstein. `MDAnalysis: A Python package for the rapid analysis of
   molecular dynamics simulations`_. In S. Benthall and S. Rostrup, editors,
   *Proceedings of the 15th Python in Science Conference*, pages 102 – 109,
   Austin, TX, 2016. SciPy.

.. _`10.1002/jcc.21787`: http://dx.doi.org/10.1002/jcc.21787

.. _`MDAnalysis: A Python package for the rapid analysis of molecular
   dynamics simulations`:
   http://conference.scipy.org/proceedings/scipy2016/oliver_beckstein.html


.. _references-components:

Citations for included algorithms and modules
=============================================

If you use the RMSD calculation (:mod:`MDAnalysis.analysis.rms`) or alignment
code (:mod:`MDAnalysis.analysis.align`) that uses the
:mod:`~MDAnalysis.core.qcprot` module please also cite [Theobald2005b]_ and
[Liu2010b]_.

.. [Theobald2005b] Douglas L. Theobald. Rapid calculation of RMSD using a
   quaternion-based characteristic polynomial. *Acta Crystallographica A*
   **61** (2005), 478-480.

.. [Liu2010b] Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald. Fast
   determination of the optimal rotational matrix for macromolecular
   superpositions. *J. Comput. Chem.* **31** (2010), 1561–1563.

If you use the helix analysis algorithm HELANAL_ in
:mod:`MDAnalysis.analysis.helanal` please cite [Bansal2000b]_.

.. [Bansal2000b] Bansal M, Kumar S, Velavan R. HELANAL — A program to
   characterise helix geometry in proteins. *J. Biomol. Struct. Dyn.* **17**
   (2000), 811–819

.. _HELANAL: http://www.ccrnp.ncifcrf.gov/users/kumarsan/HELANAL/helanal.html

If you use the GNM trajectory analysis code in
:mod:`MDAnalysis.analysis.gnm` please cite [Hall2007b]_.

.. [Hall2007b] Benjamin A. Hall, Samantha L. Kaye, Andy Pang, Rafael Perera, and
   Philip C. Biggin. Characterization of Protein Conformational States by
   Normal-Mode Frequencies. *JACS* **129** (2007), 11394–11401.

If you use the water analysis code in
:mod:`MDAnalysis.analysis.waterdynamics` please cite [Araya-Secchi2014b]_.

.. [Araya-Secchi2014b] R. Araya-Secchi., Tomas Perez-Acle, Seung-gu Kang, Tien
   Huynh, Alejandro Bernardin, Yerko Escalona, Jose-Antonio Garate,
   Agustin D. Martinez, Isaac E. Garcia, Juan C. Saez, Ruhong
   Zhou. Characterization of a novel water pocket inside the human Cx26
   hemichannel structure. *Biophysical Journal*, **107** (2014), 599-612.

If you use the Path Similarity Analysis (PSA) code in
:mod:`MDAnalysis.analysis.psa` please cite [Seyler2015b]_.

.. [Seyler2015b] Seyler SL, Kumar A, Thorpe MF, Beckstein O. Path Similarity
  Analysis: A Method for Quantifying Macromolecular Pathways. *PLoS
  Comput Biol* **11** (2015), e1004568. doi: `10.1371/journal.pcbi.1004568`_

.. _`10.1371/journal.pcbi.1004568`: http://doi.org/10.1371/journal.pcbi.1004568

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


Thanks!


