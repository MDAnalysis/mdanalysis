=========================
 MDAnalysis Applications
=========================

The applications/ directory contains self-contained Python packages
that make use of MDAnalysis as their primary dependency. Each package
should be installable on its own but make use of dependency mechanisms
to automatically install MDAnalysis (and other required packages). 

The idea is that an interested use can easily install a particular
application without having to worry too much about MDAnalysis itself.

Each application should contain

 * code in a Python module or package
 * setup.py
 * installable scripts (this is how a user would typically make use of
   the application)
 * documentation
   - usage information and an example
   - data to execute the example
   - citation information (if applicable)



If you want your application hosted with MDAnalysis then get in touch
on the developer mailing list
https://groups.google.com/forum/?fromgroups#!forum/mdnalysis-devel


