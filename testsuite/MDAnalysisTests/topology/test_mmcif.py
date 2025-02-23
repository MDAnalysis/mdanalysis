import MDAnalysis as mda
import pytest
from MDAnalysis.coordinates.MMCIF import HAS_GEMMI

from MDAnalysisTests.datafiles import MMCIF as MMCIF_FOLDER


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize(
    "mmcif_filename,n_chains",
    [
        (f"{MMCIF_FOLDER}/1YJP.cif", 1),
        (f"{MMCIF_FOLDER}/1YJP.cif.gz", 1),
        (f"{MMCIF_FOLDER}/7ETN.cif", 2),
        (f"{MMCIF_FOLDER}/7ETN.cif.gz", 2),
    ],
)
def test_chains(mmcif_filename, n_chains):
    u = mda.Universe(mmcif_filename)
    assert len(u.segments) == n_chains


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize(
    "mmcif_filename,sequence",
    [
        (
            f"{MMCIF_FOLDER}/1YJP.cif",
            ["GLY", "ASN", "ASN", "GLN", "GLN", "ASN", "TYR"],
        ),
        (
            f"{MMCIF_FOLDER}/1YJP.cif.gz",
            ["GLY", "ASN", "ASN", "GLN", "GLN", "ASN", "TYR"],
        ),
        (f"{MMCIF_FOLDER}/7ETN.cif", ["PRO", "PHE", "LEU", "ILE"]),
        (f"{MMCIF_FOLDER}/7ETN.cif.gz", ["PRO", "PHE", "LEU", "ILE"]),
    ],
)
def test_sequence(mmcif_filename, sequence):
    u = mda.Universe(mmcif_filename)
    in_structure = [
        str(res.resname)
        for res in u.select_atoms("protein and chainid A").residues
    ]
    assert in_structure == sequence, ":".join(in_structure)


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_wrong_format():
    with pytest.raises(ValueError):
        mda.Universe(f"{MMCIF_FOLDER}/1YJP_invalid.cif")


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_multimodel_warning_msg():
    with pytest.warns(
        UserWarning,
        match=r"MMCIF model .+ contains .+ different models, but only the first one will be used to assign the topology",
    ):
        mda.topology.MMCIFParser.MMCIFParser(
            f"{MMCIF_FOLDER}/multimodel_warning.cif"
        ).parse()
