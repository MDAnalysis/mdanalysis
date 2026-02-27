import MDAnalysis as mda
from MDAnalysis.tests.datafiles import GSD
import numpy as np

print("Loading Universe...")
u = mda.Universe(GSD)
print(f"Number of frames: {len(u.trajectory)}")

# This should work (normal int)
print("Indexing with int(0):")
ts_int = u.trajectory[0]
print("frame:", ts_int.frame)

# This should fail (np.int64)
print("\nIndexing with np.int64(0):")
frame_np = np.int64(0)
ts_np = u.trajectory[frame_np] 