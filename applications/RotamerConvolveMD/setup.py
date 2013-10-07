#! /usr/bin/python
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

if __name__ == '__main__':
    RELEASE = "1.0"
    with open("README") as summary:
        LONG_DESCRIPTION = summary.read()
    CLASSIFIERS = ['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ]
    setup(name              = 'rotamers',
          version           = RELEASE,
          description       = 'Analysis of spin label distances over structural ensembles',
          author            = 'Philip W. Fowler',
          author_email      = 'philip.fowler@bioch.ox.ac.uk',
          url               = 'http://mdanalysis.googlecode.com/',
          requires          = ['numpy (>=1.6)', 'MDAnalysis (>0.7.7)'],
          provides          = ['rotamers'],
          license           = 'GPL 2',
          packages          = find_packages(exclude=['scripts', 'rotamers/data']),
          package_data      = {'rotamers': ['data/*.pdb', 'data/*.dcd' ,'data/*.dat']},
          scripts           = ['scripts/convolve-mtss-rotamers.py'],
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          zip_safe = True,
          )

