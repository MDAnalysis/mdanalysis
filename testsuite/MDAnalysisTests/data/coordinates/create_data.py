import MDAnalysis as mda
import numpy as np
from six.moves import range


def create_test_trj(uni, fname):
    n_atoms = uni.atoms.n_atoms
    pos = np.arange(3 * n_atoms).reshape(n_atoms, 3)
    uni.trajectory.ts.dt = 1
    orig_box = np.array([81.1, 82.2, 83.3, 75, 80, 85], dtype=np.float32)
    uni.trajectory.ts.dimensions = orig_box
    print(uni.trajectory)
    print(uni.trajectory.ts.__class__)
    with mda.Writer(fname, n_atoms) as w:
        for i in range(5):
            uni.atoms.positions = 2 ** i * pos
            uni.trajectory.ts.time = i
            uni.trajectory.ts.velocities = uni.atoms.positions / 10
            uni.trajectory.ts.forces = uni.atoms.positions / 100
            uni.trajectory.ts.frame = i
            uni.trajectory.ts.dimensions[:3] = orig_box[:3] + i
            uni.trajectory.ts.dimensions[3:] = orig_box[3:] + i * 0.1
            print(uni.trajectory.ts.dimensions)
            w.write(uni)


def main():
    pdb = 'test_topology.pdb'
    u = mda.Universe(pdb)

    create_test_trj(u, 'test.xyz')
    create_test_trj(u, 'test.xtc')
    create_test_trj(u, 'test.trr')

if __name__ == '__main__':
    main()
