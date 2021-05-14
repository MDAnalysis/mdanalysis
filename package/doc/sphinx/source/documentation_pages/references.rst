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
[Michaud-Agrawal2011]_ and [Gowers2016]_.

(We are currently asking you to cite *both* papers if at all possible
because the 2016 paper describes many updates to the original 2011
paper and neither paper on its own provides a comprehensive
description of the library. We will publish a complete self-contained
paper with the upcoming 1.0 release of MDAnalysis, which will then
supersede these two citations.)


.. [Michaud-Agrawal2011] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf,
   and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics
   Simulations. *J. Comput. Chem.* **32** (2011),
   2319–2327. doi:`10.1002/jcc.21787`_

.. [Gowers2016] R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N.
   Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney,
   and O. Beckstein. `MDAnalysis: A Python package for the rapid analysis of
   molecular dynamics simulations`_. In S. Benthall and S. Rostrup, editors,
   *Proceedings of the 15th Python in Science Conference*, pages 98-105,
   Austin, TX, 2016. SciPy. doi:`10.25080/Majora-629e541a-00e`_

.. _`10.1002/jcc.21787`: http://dx.doi.org/10.1002/jcc.21787
.. _`10.25080/Majora-629e541a-00e`:
   https://doi.org/10.25080/Majora-629e541a-00e

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

If you use the hydrogen bond analysis code in
:mod:`MDAnalysis.analysis.hydrogenbonds.hbond_analysis` please cite [Smith2019]_.

.. [Smith2019] P. Smith, R. M. Ziolek, E. Gazzarrini, D. M. Owen, and C. D. Lorenz.
   On the interaction of hyaluronic acid with synovial fluid lipid membranes. *PCCP*
   **21** (2019), 9845-9857. doi:  `10.1039/C9CP01532A`_

.. _`10.1039/C9CP01532A`: http://dx.doi.org/10.1039/C9CP01532A

If you use :meth:`~MDAnalysis.analysis.pca.PCA.rmsip` or 
:func:`~MDAnalysis.analysis.pca.rmsip` please cite [Amadei1999]_ and 
[Leo-Macias2004]_ .

.. [Amadei1999] Amadei, A., Ceruso, M. A. & Nola, A. D. 
   On the convergence of the conformational coordinates basis set obtained by the essential dynamics analysis of proteins’ molecular dynamics simulations. 
   *Proteins: Structure, Function, and Bioinformatics* **36**, 419–424 (1999).
   doi: `10.1002/(SICI)1097-0134(19990901)36:4<419::AID-PROT5>3.0.CO;2-U`_

.. _`10.1002/(SICI)1097-0134(19990901)36:4<419::AID-PROT5>3.0.CO;2-U`: https://doi.org/10.1002/(SICI)1097-0134(19990901)36:4%3C419::AID-PROT5%3E3.0.CO;2-U

.. [Leo-Macias2004] Leo-Macias, A., Lopez-Romero, P., Lupyan, D., Zerbino, D. & Ortiz, A. R. 
   An Analysis of Core Deformations in Protein Superfamilies. 
   *Biophys J* **88**, 1291–1299 (2005). doi: `10.1529/biophysj.104.052449`_

.. _`10.1529/biophysj.104.052449`: https://dx.doi.org/10.1529%2Fbiophysj.104.052449

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

