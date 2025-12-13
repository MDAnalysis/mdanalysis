# mypy: ignore-errors
import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from MDAnalysis.tests.datafiles import RANDOM_WALK, RANDOM_WALK_TOPO
from MDAnalysis.transformations import NoJump
import matplotlib.pyplot as plt

print("Using input files:")
print("Topology:", RANDOM_WALK_TOPO)
print("Trajectory:", RANDOM_WALK)

u = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)

# --- FIX: define box dimensions BEFORE NoJump runs ---
box = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]


def set_box(ts):
    ts.dimensions = box
    return ts


# Apply transformations in correct order
u.trajectory.add_transformations(
    set_box,  # must come first
    NoJump(u),  # requires PBC
)

# Compute MSD
msd = EinsteinMSD(u, select="all", msd_type="xyz", fft=False)
msd.run()

lagtimes = msd.results.delta_t_values
msd_values = msd.results.timeseries

plt.figure(figsize=(6, 4))
plt.plot(lagtimes, msd_values, label="Computed MSD")
plt.plot(lagtimes, 6 * lagtimes, "--", label="Theoretical 3D MSD (6τ)")
plt.xlabel("Lag time (τ)")
plt.ylabel("MSD")
plt.legend()
plt.title("Mean Squared Displacement with NoJump")

plt.savefig("example/msd_nojump.png", dpi=150, bbox_inches="tight")

plt.close()

print("Saved plot to msd_nojump.png")
