import MDAnalysis as mda
import numpy as np
from six.moves import range


def create_test_trj(uni, fname, wrong_natoms=False):
    n_atoms = uni.atoms.n_atoms
    if wrong_natoms:
        n_atoms -= 2
    pos = np.arange(3 * n_atoms).reshape(n_atoms, 3)
    uni.trajectory.ts.dt = 1
    with mda.Writer(fname, n_atoms) as w:
        for i in range(5):
            uni.atoms.positions = 2 ** i * pos
            uni.trajectory.ts.time = i
            w.write(uni)


def main():
    pdb = 'test_topology.pdb'
    u = mda.Universe(pdb)

    create_test_trj(u, 'test.xyz')

if __name__ == '__main__':
    main()
