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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import pytest

import MDAnalysis as mda


# grab all known writers
# sort so test order is predictable for parallel tests
writers = sorted(set(mda._MULTIFRAME_WRITERS.values()) |
                 set(mda._SINGLEFRAME_WRITERS.values()),
                 key=lambda x: x.__name__)
known_ts_haters = [
    mda.coordinates.MOL2.MOL2Writer,
    mda.coordinates.PDB.PDBWriter,
    mda.coordinates.PDB.MultiPDBWriter,
    mda.coordinates.PQR.PQRWriter,
    mda.coordinates.PDBQT.PDBQTWriter,
    mda.coordinates.LAMMPS.DATAWriter,
    mda.coordinates.CRD.CRDWriter,
]


@pytest.mark.parametrize('writer', writers)
def test_ts_error(writer, tmpdir):
    u = mda.Universe.empty(10, trajectory=True)

    if writer == mda.coordinates.chemfiles.ChemfilesWriter:
        # chemfiles Writer exists but doesn't work without chemfiles
        if not mda.coordinates.chemfiles.check_chemfiles_version():
            pytest.skip("Chemfiles not available")
        fn = str(tmpdir.join('out.xtc'))
    elif writer == mda.coordinates.LAMMPS.DATAWriter:
        pytest.skip("DATAWriter requires integer atom types")
    elif writer == mda.coordinates.H5MD.H5MDWriter and not mda.coordinates.H5MD.HAS_H5PY:
        pytest.skip("skipping H5MDWriter test because h5py is not installed")
    else:
        fn = str(tmpdir.join('out.traj'))

    if (writer == mda.coordinates.PDB.PDBWriter or
        writer == mda.coordinates.PDB.MultiPDBWriter):
        errmsg = "PDBWriter cannot write Timestep objects directly"
    else:
        errmsg = "neither an AtomGroup or Universe"

    with pytest.raises(TypeError, match=errmsg):
        with writer(fn, n_atoms=u.atoms.n_atoms) as w:
            w.write(u.trajectory.ts)


@pytest.mark.parametrize('writer', writers)
def test_write_with_atomgroup(writer, tmpdir):
    u = mda.Universe.empty(10, trajectory=True)

    if writer == mda.coordinates.chemfiles.ChemfilesWriter:
        # chemfiles Writer exists but doesn't work without chemfiles
        if not mda.coordinates.chemfiles.check_chemfiles_version():
            pytest.skip("Chemfiles not available")
        fn = str(tmpdir.join('out.xtc'))
    elif writer == mda.coordinates.MOL2.MOL2Writer:
        pytest.skip("MOL2 only writes MOL2 back out")
    elif writer == mda.coordinates.LAMMPS.DATAWriter:
        pytest.skip("DATAWriter requires integer atom types")
    elif writer == mda.coordinates.H5MD.H5MDWriter and not mda.coordinates.H5MD.HAS_H5PY:
        pytest.skip("skipping H5MDWriter test because h5py is not installed")
    else:
        fn = str(tmpdir.join('out.traj'))

    # H5MDWriter raises ValueError if the trajectory has no units and
    # convert_units is True
    convert_units = (writer != mda.coordinates.H5MD.H5MDWriter)

    with writer(fn, n_atoms=u.atoms.n_atoms, convert_units=convert_units) as w:
        w.write(u.atoms)


@pytest.mark.parametrize('writer', writers)
def test_write_with_universe(writer, tmpdir):
    u = mda.Universe.empty(10, trajectory=True)

    if writer == mda.coordinates.chemfiles.ChemfilesWriter:
        # chemfiles Writer exists but doesn't work without chemfiles
        if not mda.coordinates.chemfiles.check_chemfiles_version():
            pytest.skip("Chemfiles not available")
        fn = str(tmpdir.join('out.xtc'))
    elif writer == mda.coordinates.MOL2.MOL2Writer:
        pytest.skip("MOL2 only writes MOL2 back out")
    elif writer == mda.coordinates.LAMMPS.DATAWriter:
        pytest.skip("DATAWriter requires integer atom types")
    elif writer == mda.coordinates.H5MD.H5MDWriter and not mda.coordinates.H5MD.HAS_H5PY:
        pytest.skip("skipping H5MDWriter test because h5py is not installed")
    else:
        fn = str(tmpdir.join('out.traj'))

    # H5MDWriter raises ValueError if the trajectory has no units and
    # convert_units is True
    convert_units = (writer != mda.coordinates.H5MD.H5MDWriter)

    with writer(fn, n_atoms=u.atoms.n_atoms, convert_units=convert_units) as w:
        w.write(u.atoms)
