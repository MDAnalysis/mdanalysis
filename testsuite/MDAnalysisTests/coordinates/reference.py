# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
import numpy as np

from MDAnalysisTests import datafiles
from MDAnalysisTests.datafiles import (LAMMPSdata,
                                       LAMMPSdata2, LAMMPSdcd2,
                                       LAMMPSdata_mini,
                                       )


class RefAdKSmall(object):
    """Mixin class to provide comparison numbers.

    Based on small PDB with AdK (:data:`PDB_small`).

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    filename = datafiles.PDB_small
    ref_coordinates = {
        # G11:CA, copied frm adk_open.pdb
        'A10CA': np.array([-1.198, 7.937, 22.654]),
    }
    ref_distances = {'endtoend': 11.016959}
    ref_E151HA2_index = 2314
    ref_n_atoms = 3341
    ref_charmm_totalcharge = -4.0
    ref_charmm_Hcharges = [0.33] + 203 * [0.31]
    ref_charmm_ArgCAcharges = 13 * [0.07]
    ref_charmm_ProNcharges = 10 * [-0.29]
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    ref_volume = 0.0


class RefAdK(object):
    """Mixin class to provide comparison numbers.

    Based on PDB/GRO with AdK in water + Na+ (:data:`PDB`).

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    filename = datafiles.PDB
    ref_coordinates = {
        # Angstroem as MDAnalysis unit!!
        'A10CA': np.array([62.97600174, 62.08800125, 20.2329998]),
    }
    ref_distances = {'endtoend': 9.3513174}
    ref_E151HA2_index = 2314
    ref_n_atoms = 47681
    ref_Na_sel_size = 4
    # CRYST1 80.017   80.017   80.017  60.00  60.00  90.00
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    #ref_volume = 362270.0  # computed with Gromacs ## NOT EXACT!
    ref_volume = 362269.520669292


class Ref2r9r(object):
    """Mixin class to provide comparison numbers.

    Based on S6 helices of chimeric Kv channel

    .. Note::

       All distances must be in ANGSTROEM as this is the MDAnalysis
       default unit. All readers must return Angstroem by default.
    """
    ref_n_atoms = 1284
    ref_sum_centre_of_geometry = -98.24146
    ref_n_frames = 10


class RefACHE(object):
    """Mixin class to provide comparison numbers.

    ACHE peptide

    # COM check in VMD::

        set p [atomselect top "not water"]
        set total {0 0 0};
        for {set i 0} {$i < 11} {incr i} {
           $p frame $i; set total [vecadd $total [measure center $p]]}

        puts [vecsum $total]
        # 472.2592159509659

    """
    ref_n_atoms = 252
    ref_proteinatoms = ref_n_atoms
    ref_sum_centre_of_geometry = 472.2592159509659  # 430.44807815551758
    ref_n_frames = 11
    ref_periodic = False


class RefCappedAla(object):
    """Mixin class to provide comparison numbers.

    Capped Ala in water

    # COM check in VMD (load trajectory as *AMBER with periodic box*!)::

        set p [atomselect top "not water"]
        set total {0 0 0};
        for {set i 0} {$i < 11} {incr i} {
           $p frame $i; set total [vecadd $total [measure center $p]]}

        puts [vecsum $total]
        # 686.276834487915

    """
    ref_n_atoms = 5071
    ref_proteinatoms = 22
    ref_sum_centre_of_geometry = 686.276834487915
    ref_n_frames = 11
    ref_periodic = True


class RefVGV(object):
    """Mixin class to provide comparison numbers.

    Computed from bala.trj::

      w = MDAnalysis.Universe(PRMncdf, TRJncdf)
      ref_n_atoms = len(w.atoms) ref_proteinatoms = len(w.select_atoms("protein"))
      ref_sum_centre_of_geometry = np.sum([protein.center_of_geometry()
                                           for ts in w.trajectory])
    """
    topology = datafiles.PRMncdf
    filename = datafiles.NCDF
    ref_n_atoms = 2661
    ref_proteinatoms = 50
    ref_sum_centre_of_geometry = 1552.9125
    ref_n_frames = 30
    ref_periodic = True

class RefTZ2(object):
    """Reference values for the cpptraj testcase tz2.truncoct.nc

    Used under the GPL v3.
    """
    topology = datafiles.PRM7
    filename = datafiles.NCDFtruncoct
    ref_n_atoms = 5827
    ref_proteinatoms = 217
    ref_sum_centre_of_geometry = -68.575745
    ref_n_frames = 10
    ref_periodic = True
    

class RefTRZ(object):
    #    ref_coordinates = {}
    #    ref_distances = {'endtoend': }
    ref_n_atoms = 8184
    ref_dimensions = np.array([55.422830581665039, 55.422830581665039,
                               55.422830581665039, 90., 90., 90.],
                              dtype=np.float32)
    ref_volume = 170241.762765
    ref_n_frames = 6
    ref_coordinates = np.array([72.3163681, -130.31130981, 19.97969055],
                               dtype=np.float32)
    ref_velocities = np.array([[14.83297443, 18.02611542, 6.07733774]],
                              dtype=np.float32)
    ref_delta = 0.001
    ref_time = 0.01
    ref_title = ('ABCDEFGHIJKLMNOPQRSTUVWXYZ12345678901234'
                 'ABCDEFGHIJKLMNOPQRSTUVWXYZ12345678901234')


class RefLAMMPSData(object):
    filename = LAMMPSdata
    n_atoms = 18364
    pos_atom1 = np.array([11.89985657, 48.4455719, 19.09719849],
                         dtype=np.float32)
    vel_atom1 = np.array([-0.005667593, 0.00791380978, -0.00300779533],
                         dtype=np.float32)
    dimensions = np.array([55.42282867, 55.42282867, 55.42282867, 90., 90., 90.
                           ],
                          dtype=np.float32)

class RefLAMMPSDataDCD(object):
    format = "LAMMPS"
    topology = LAMMPSdata2
    trajectory = LAMMPSdcd2
    n_atoms = 12421
    n_frames = 5
    dt = 0.5  # ps per frame
    mean_dimensions = np.array(
        [ 50.66186142,  47.18824387,  52.33762741,
          90.        ,  90.        ,  90.        ], dtype=np.float32)


class RefLAMMPSDataMini(object):
    filename = LAMMPSdata_mini
    n_atoms = 1
    pos_atom1 = np.array([11.89985657, 48.4455719, 19.09719849],
                         dtype=np.float32)
    vel_atom1 = np.array([-0.005667593, 0.00791380978, -0.00300779533],
                         dtype=np.float32)
    dimensions = np.array([60., 50., 30., 90., 90., 90.], dtype=np.float32)


class RefLAMMPSDataAdditionalColumns(object):
    q = np.array([2.58855e-03, 6.91952e-05, 1.05548e-02, 4.20319e-03,
                  9.19172e-03, 4.79777e-03, 6.36864e-04, 5.87125e-03,
                  -2.18125e-03, 6.88910e-03])
    p = np.array(5 * [1.1, 1.2])
