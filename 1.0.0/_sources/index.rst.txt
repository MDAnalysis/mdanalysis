.. MDAnalysis User Guide documentation master file, created by
   sphinx-quickstart on Wed Aug 14 14:47:38 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=================================================
Welcome to MDAnalysis User Guide's documentation!
=================================================

**MDAnalysis version:** |MDAnalysis_version|

**Last updated:** |today|

**MDAnalysis** (`www.mdanalysis.org`_) is a Python
toolkit to analyse molecular dynamics files and trajectories in :ref:`many popular formats <formats>`. MDAnalysis can write
most of these formats, too, together with atom selections for use in :ref:`visualisation tools or other analysis programs <selection-exporters>`.
It provides a fast framework for :ref:`complex analysis tasks <analysis>`, 
as well as flexible tooling to construct your own analyses.

.. _`www.mdanalysis.org`: https://www.mdanalysis.org

Why MDAnalysis?
===============

The typical use case for MDAnalysis is to manipulate or analyse molecular dynamics trajectories. The library focuses on two key features:

   * **Memory efficiency.**
     The size of trajectory data can quickly overwhelm the memory resources of your computer. 
     MDAnalysis typically accesses your trajectory by only loading data for one frame at a time. 
     This allows you to work with trajectories of any length without difficulty.

   * **Flexibility.**
     MDAnalysis is constructed to be easily extensible.
     If an analysis method is not already available in MDAnalysis, 
     you can write your own custom trajectory analysis with the building blocks provided. 
     If you need to read in a custom file format, you can construct your own Reader or Parser that will automatically get picked up when MDAnalysis is constructing a Universe from files. You can create and add your own labels for atoms, residues, or segments (called :ref:`topology attributes <topology-attributes>`) and relationships between atoms (e.g. bonds, angles). 


Participating
=============

MDAnalysis welcomes all contributions from its users. There are many ways you can help improve MDAnalysis, from asking questions on the `mdnalysis-discussion`_ mailing list, to raising issues on the `Issue Tracker`_, to adding your own code. Please see :ref:`contributing` for an introduction and guide to contributing to the code and documentation. 

.. important::

   **Ground rules and expectations**

   The MDAnalysis community subscribes to a `Code of Conduct`_. By participating in this project and community, you agree to abide by its terms. Please read it.

   In general, we expect you to **be kind and thoughtful in your conversations around this project.** We all come from different backgrounds and projects, which means we will not always agree. Try to listen and understand why others hold their viewpoints in discussions. Rather than blaming each other, focus on helping to resolve issues and learning from mistakes.



--------------
Communications
--------------

Questions and discussions about MDAnalysis take place on the mailing lists and this repository’s `Issue Tracker`_. Anybody is welcome to join these conversations. Please ask questions about the usage of MDAnalysis on the `mdnalysis-discussion`_ mailing list, and report problems on the `Issue Tracker`_.

Wherever possible, do not take these conversations to private channels, including contacting the maintainers directly. Keeping communication public means everybody can benefit and learn from the conversation.

.. _`mdnalysis-discussion`:
   https://groups.google.com/group/mdnalysis-discussion
.. _`Code of Conduct`: https://www.mdanalysis.org/pages/conduct/
.. _`Issue Tracker`: https://github.com/MDAnalysis/mdanalysis/issues

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:

   installation
   examples/quickstart
   faq
   examples/README

.. toctree::
   :maxdepth: 1
   :caption: Data structures
   :hidden:

   universe
   atomgroup
   groups_of_atoms
   selections
   topology_system

.. toctree::
   :maxdepth: 1
   :caption: Trajectories
   :hidden:

   trajectories/trajectories
   trajectories/slicing_trajectories
   trajectories/transformations
   units

.. toctree::
   :maxdepth: 1
   :caption: Input/output
   :hidden:

   reading_and_writing
   formats/index
   formats/guessing
   formats/auxiliary
   formats/selection_exporters
   formats/format_reference

.. toctree::
   :maxdepth: 1
   :caption: Analysis
   :hidden:

   examples/analysis/README
   examples/analysis/custom_trajectory_analysis

.. toctree::
   :maxdepth: 1
   :caption: Advanced
   :hidden:

   standard_selections
   advanced_topology
   datasets

.. toctree::
   :maxdepth: 1
   :caption: Contributing
   :hidden:

   contributing
   contributing_code
   contributing_docs
   preparing_releases_and_hotfixes
   module_imports
   testing
   references
