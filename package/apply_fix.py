import re, os, sys

files_and_replacements = {
    "MDAnalysis/__init__.py":                               '"MDAnalysis.__init__"',
    "MDAnalysis/analysis/align.py":                         '"MDAnalysis.analysis.align"',
    "MDAnalysis/analysis/atomicdistances.py":               '"MDAnalysis.analysis.atomicdistances"',
    "MDAnalysis/analysis/contacts.py":                      '"MDAnalysis.analysis.contacts"',
    "MDAnalysis/analysis/density.py":                       '"MDAnalysis.analysis.density"',
    "MDAnalysis/analysis/diffusionmap.py":                  '"MDAnalysis.analysis.diffusionmap"',
    "MDAnalysis/analysis/distances.py":                     '"MDAnalysis.analysis.distances"',
    "MDAnalysis/analysis/gnm.py":                           '"MDAnalysis.analysis.GNM"',
    "MDAnalysis/analysis/hydrogenbonds/wbridge_analysis.py":'"MDAnalysis.analysis.WaterBridgeAnalysis"',
    "MDAnalysis/analysis/legacy/x3dna.py":                  '"MDAnalysis.analysis.x3dna"',
    "MDAnalysis/analysis/msd.py":                           '"MDAnalysis.analysis.msd"',
    "MDAnalysis/analysis/rms.py":                           '"MDAnalysis.analysis.rmsd"',
    "MDAnalysis/converters/ParmEdParser.py":                '"MDAnalysis.converters.ParmEdParser"',
    "MDAnalysis/converters/RDKitParser.py":                 '"MDAnalysis.converters.RDKitParser"',
    "MDAnalysis/coordinates/IMD.py":                        '"MDAnalysis.coordinates.IMDReader"',
    "MDAnalysis/coordinates/PDB.py":                        '"MDAnalysis.coordinates.PBD"',
    "MDAnalysis/coordinates/TPR.py":                        '"MDAnalysis.coordinates.TPR"',
    "MDAnalysis/coordinates/TRC.py":                        '"MDAnalysis.coordinates.GROMOS11"',
    "MDAnalysis/coordinates/TRJ.py":                        '"MDAnalysis.coordinates.AMBER"',
    "MDAnalysis/coordinates/XYZ.py":                        '"MDAnalysis.coordinates.XYZ"',
    "MDAnalysis/core/universe.py":                          '"MDAnalysis.core.universe"',
    "MDAnalysis/guesser/base.py":                           '"MDAnalysis.guesser.base"',
    "MDAnalysis/topology/LAMMPSParser.py":                  '"MDAnalysis.topology.LAMMPS"',
    "MDAnalysis/topology/PDBParser.py":                     '"MDAnalysis.topology.PDBParser"',
    "MDAnalysis/topology/PSFParser.py":                     '"MDAnalysis.topology.PSF"',
    "MDAnalysis/topology/TOPParser.py":                     '"MDAnalysis.topology.TOPParser"',
    "MDAnalysis/topology/TPRParser.py":                     '"MDAnalysis.topology.TPRparser"',
}

# Run from the package/ directory
base = os.path.dirname(os.path.abspath(__file__))
# If run from repo root, prepend package/
if os.path.isdir(os.path.join(base, "package")):
    base = os.path.join(base, "package")

ok = 0
fail = 0
for rel_path, old_name in files_and_replacements.items():
    fpath = os.path.join(base, rel_path)
    if not os.path.exists(fpath):
        print(f"MISSING: {fpath}")
        fail += 1
        continue
    with open(fpath, "r", encoding="utf-8") as f:
        content = f.read()
    old = f"logging.getLogger({old_name})"
    new = "logging.getLogger(__name__)"
    if old in content:
        content = content.replace(old, new, 1)
        with open(fpath, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"FIXED : {rel_path}")
        ok += 1
    else:
        print(f"SKIP  : {rel_path}  (not found — already fixed?)")
        fail += 1

print(f"\nDone: {ok} fixed, {fail} skipped/missing")
