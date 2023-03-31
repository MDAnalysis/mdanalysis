"""
tidynamics
==========

A tiny package to compute the dynamics of stochastic and molecular simulations.

Documentation: http://lab.pdebuyl.be/tidynamics/

When using tidynamics in a publication, please cite the following paper:
    Pierre de Buyl (2018), *tidynamics: A tiny package to compute the dynamics
    of stochastic and molecular simulations*, The Journal of Open Source
    Software https://doi.org/10.21105/joss.00877


"""
from ._correlation import acf, msd, cross_displacement, correlation
import os.path

with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
    __version__ = f.read().strip()
