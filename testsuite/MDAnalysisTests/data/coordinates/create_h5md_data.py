import MDAnalysis as mda
from MDAnalysis.tests.datafiles import COORDINATES_TRR, COORDINATES_TOPOLOGY
# To create the test trajectory, install pyh5md
# from https://github.com/pdebuyl/pyh5md
# with 'pip install pyh5md'
import pyh5md

"""
This script converts the file test.trr to test.h5md, where
test.h5md will be used to as the reference trajectory for
tests using MDAnalysisTests.coordinates.base.BaseReference
"""


def create_test_trj(uni, filename):
    """uses pyh5md library to write h5md file"""
    # will eventually change this to use h5py instead

    with pyh5md.File(filename, 'w', creator='create_h5md_data.py') as f:
        trajectory = f.particles_group('trajectory')
        obsv = f.require_group('observables')

        positions = pyh5md.element(trajectory, 'position', store='time',
                                   data=uni.trajectory.ts.positions,
                                   time=True)
        f['particles/trajectory/position/value'].attrs['unit'] = 'Angstrom'
        f['particles/trajectory/position/time'].attrs['unit'] = 'ps'

        velocities = pyh5md.element(trajectory, 'velocity',
                                    data=uni.trajectory.ts.velocities,
                                    step_from=positions,
                                    store='time', time=True)
        f['particles/trajectory/velocity/value'].attrs['unit'] = 'Angstrom ps-1'
        f['particles/trajectory/velocity/time'].attrs['unit'] = 'ps'

        forces = pyh5md.element(trajectory, 'force',
                                data=uni.trajectory.ts.forces,
                                step_from=positions,
                                store='time', time=True)
        f['particles/trajectory/force/value'].attrs['unit'] = 'kJ mol-1 Angstrom-1'
        f['particles/trajectory/force/time'].attrs['unit'] = 'ps'

        data_step = pyh5md.element(obsv, 'step',
                                   data=uni.trajectory.ts.data['step'],
                                   step_from=positions,
                                   store='time', time=True)
        data_lambda = pyh5md.element(obsv, 'lambda',
                                     data=uni.trajectory.ts.data['lambda'],
                                     step_from=positions,
                                     store='time', time=True)

        trajectory.create_box(dimension=3,
                              boundary=['periodic', 'periodic', 'periodic'],
                              store='time',
                              data=uni.trajectory.ts.triclinic_dimensions,
                              step_from=positions)
        f['particles/trajectory/box/edges/value'].attrs['unit'] = 'Angstrom'

        for ts in uni.trajectory:
            trajectory.box.edges.append(uni.trajectory.ts.triclinic_dimensions,
                                        ts.frame, time=ts.time)
            positions.append(uni.trajectory.ts.positions,
                             ts.frame, time=ts.time)
            velocities.append(uni.trajectory.ts.velocities,
                              ts.frame, time=ts.time)
            forces.append(uni.trajectory.ts.forces,
                          ts.frame, time=ts.time)
            data_step.append(uni.trajectory.ts.data['step'],
                             ts.frame, time=ts.time)
            data_lambda.append(uni.trajectory.ts.data['lambda'],
                               ts.frame, time=ts.time)

def main():
    pdb = COORDINATES_TOPOLOGY
    trr = COORDINATES_TRR
    u = mda.Universe(pdb, trr)
    create_test_trj(u, 'test.h5md')


if __name__ == '__main__':
    main()
