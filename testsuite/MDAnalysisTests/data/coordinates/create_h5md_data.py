import pyh5md
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import COORDINATES_TRR, COORDINATES_TOPOLOGY

"""
This script converts the file test.trr to test.h5md, where
test.h5md will be used to as the reference trajectory for
tests using MDAnalysisTests.coordinates.base.BaseReference
"""

def create_test_trj(universe, filename):
    """uses pyh5md library to write h5md file"""
    # will eventually change this to use h5py instead

    with pyh5md.File(filename, 'w', creator='create_h5md_data.py') as f:
        # Add a trajectory group into particles group
        trajectory = f.particles_group('trajectory')

        # Add the positions, velocities, forces, n_atoms groups into the trajectory group
        trajectory_positions = pyh5md.element(trajectory,'positions', store='time',
                                              data=universe.trajectory.ts.positions, time=True)
        trajectory_velocities = pyh5md.element(trajectory, 'velocities', data=universe.trajectory.ts.velocities,
                                               step_from=trajectory_positions,store='time', time=True)
        trajectory_forces = pyh5md.element(trajectory, 'forces', data=universe.trajectory.ts.forces,
                                           step_from=trajectory_positions, store='time', time=True)
        trajectory_n_atoms = pyh5md.element(trajectory, 'n_atoms', store='fixed',
                                            data=universe.atoms.n_atoms)
        data_step = pyh5md.element(trajectory, 'data/step', store='time', data=universe.trajectory.ts.data['step'])
        data_lambda = pyh5md.element(trajectory, 'data/lambda', store='time', data=universe.trajectory.ts.data['lambda'])
        data_dt = pyh5md.element(trajectory, 'data/dt', store='time', data=universe.trajectory.ts.data['dt'])

        trajectory.create_box(dimension=3, boundary=['periodic', 'periodic', 'periodic'],
                              store='time', data=universe.trajectory.ts.triclinic_dimensions,
                              step_from=trajectory_positions)

        for ts in universe.trajectory:
            trajectory.box.edges.append(universe.trajectory.ts.triclinic_dimensions, ts.frame, time=ts.time)
            trajectory_positions.append(universe.trajectory.ts.positions, ts.frame, time=ts.time)
            trajectory_velocities.append(universe.trajectory.ts.velocities, ts.frame, time=ts.time)
            trajectory_forces.append(universe.trajectory.ts.forces, ts.frame, time=ts.time)
            data_step.append(universe.trajectory.ts.data['step'], ts.frame, time=ts.time)
            data_lambda.append(universe.trajectory.ts.data['lambda'], ts.frame, time=ts.time)
            data_dt.append(universe.trajectory.ts.data['dt'], ts.frame, time=ts.time)

def main():
    pdb = COORDINATES_TOPOLOGY
    trr = COORDINATES_TRR
    u = mda.Universe(pdb, trr)
    create_test_trj(u, 'test.h5md')

if __name__ == '__main__':
    main()
