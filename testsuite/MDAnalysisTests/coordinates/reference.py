import numpy as np

from MDAnalysisTests import datafiles
from MDAnalysisTests.datafiles import (PDB_small, PDB, PDB_full, LAMMPSdata,
                                       LAMMPSdata2, LAMMPSdcd2,
                                       LAMMPSdata_mini, PSF_TRICLINIC,
                                       DCD_TRICLINIC, PSF_NAMD_TRICLINIC,
                                       DCD_NAMD_TRICLINIC)


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
    ref_volume = 362270.0  # computed with Gromacs


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


class Ref4e43(object):
    """Mixin class for a clean Protein Databank PDB file"""
    filename = datafiles.PDB_full
    header = 'HYDROLASE                               11-MAR-12   4E43'
    title = ['HIV PROTEASE (PR) DIMER WITH ACETATE IN EXO SITE AND PEPTIDE '
             'IN ACTIVE', '2 SITE']
    compnd = ['MOL_ID: 1;',
              '2 MOLECULE: PROTEASE;',
              '3 CHAIN: A, B;',
              '4 ENGINEERED: YES;',
              '5 MUTATION: YES;',
              '6 MOL_ID: 2;',
              '7 MOLECULE: RANDOM PEPTIDE;',
              '8 CHAIN: C;',
              '9 ENGINEERED: YES;',
              '10 OTHER_DETAILS: UNKNOWN IMPURITY', ]
    num_remarks = 333
    # only first 5 remarks for comparison
    nmax_remarks = 5
    remarks = [
        '2',
        '2 RESOLUTION.    1.54 ANGSTROMS.',
        '3',
        '3 REFINEMENT.',
        '3   PROGRAM     : REFMAC 5.5.0110',
    ]


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
    n_atoms = 18360
    pos_atom1 = np.array([11.89985657, 48.4455719, 19.09719849],
                         dtype=np.float32)
    vel_atom1 = np.array([-5.667593, 7.91380978, -3.00779533],
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
    vel_atom1 = np.array([-5.667593, 7.91380978, -3.00779533],
                         dtype=np.float32)
    dimensions = np.array([60., 50., 30., 90., 90., 90.], dtype=np.float32)


class RefCHARMMtriclinicDCD(object):
    topology = PSF_TRICLINIC
    trajectory = DCD_TRICLINIC
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    # dcd starts at t = 1ps
    ref_dimensions = np.array([
        [1., 35.44604, 35.06156, 34.1585, 91.32802, 61.73521, 44.40703],
        [2., 34.65957, 34.22689, 33.09897, 90.56206, 61.79192, 44.14549],
        [3., 34.52772, 34.66422, 33.53881, 90.55859, 63.11228, 40.14044],
        [4., 34.43749, 33.38432, 34.02133, 88.82457, 64.98057, 36.77397],
        [5., 33.73129, 32.47752, 34.18961, 89.88102, 65.89032, 36.10921],
        [6., 33.78703, 31.90317, 34.98833, 90.03092, 66.12877, 35.07141],
        [7., 33.24708, 31.18271, 34.9654, 93.11122, 68.17743, 35.73643],
        [8., 32.92599, 30.31393, 34.99197, 93.89051, 69.3799, 33.48945],
        [9., 32.15295, 30.43056, 34.96157, 96.01416, 71.50115, 32.56111],
        [10., 31.99748, 30.21518, 35.24292, 95.85821, 71.08429, 31.85939]
    ])


class RefNAMDtriclinicDCD(object):
    topology = PSF_NAMD_TRICLINIC
    trajectory = DCD_NAMD_TRICLINIC
    # vmd topology trajectory
    # molinfo 0 get {a b c alpha beta gamma}
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    ref_dimensions = np.array([
        [1., 38.426594, 38.393101, 44.759800, 90.000000, 90.000000, 60.028915],
    ])

