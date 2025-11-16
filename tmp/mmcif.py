import MDAnalysis as mda
from io import StringIO
from pathlib import Path
from MDAnalysis.lib.util import NamedStream

u = mda.Universe(
    NamedStream(
        StringIO(open("/home/law/workspace/mdanalysis/tmp/1BD2.cif").read()),
        "tmp.cif",
    ),
    format="mmcif",
)

tmp_brk = u.select_atoms(f"resid 54 and segid D")

print(tmp_brk)

print(tmp_brk.select_atoms("name CA"))
