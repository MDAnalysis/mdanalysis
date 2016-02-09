
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.analysis.contacts import calculate_contacts
import numpy as np
import pandas as pd

ref = Universe("conf_protein.gro")
u = Universe("conf_protein.gro", "traj_protein_0.xtc")

x = len(ref.select_atoms("protein"))
selA = "not name H* and resid 72-95 and bynum {}:{}".format(1, x//2)
selB = "not name H* and resid 72-95 and bynum {}:{}".format(x//2, x)


df = calculate_contacts(ref, u, selA, selB)
print(df)
