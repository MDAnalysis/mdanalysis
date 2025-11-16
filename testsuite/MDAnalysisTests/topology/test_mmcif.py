import MDAnalysis as mda
import pytest
from MDAnalysis.coordinates.MMCIF import HAS_GEMMI

from MDAnalysisTests.datafiles import MMCIF as MMCIF_FOLDER


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize(
    "basename",
    [
        f"{MMCIF_FOLDER}/1BD2_short",
        f"{MMCIF_FOLDER}/1BD2",
        f"{MMCIF_FOLDER}/3PWP",
        f"{MMCIF_FOLDER}/3KPR",
    ],
)
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_legacy_pdb_vs_mmcif(basename):
    u_cif = mda.Universe(f"{basename}.cif.gz")
    u_pdb = mda.Universe(f"{basename}.pdb.gz")
    assert len(u_pdb.atoms) == len(u_cif.atoms)
    assert len(u_pdb.select_atoms("segid *")) == len(
        u_cif.select_atoms("segid *")
    )

    assert len(u_pdb.select_atoms("protein")) == len(
        u_cif.select_atoms("protein")
    )

    assert len(u_pdb.select_atoms("name CA and segid D")) == len(
        u_cif.select_atoms("name CA and segid D")
    )
    for segment in "ABCDE":
        for resid in [1, 10, 54, 72]:
            assert len(u_pdb.select_atoms(f"segid {segment}")) == len(
                u_cif.select_atoms(f"segid {segment}")
            )
            assert len(
                u_pdb.select_atoms(f"segid {segment} and resid {resid}")
            ) == len(u_cif.select_atoms(f"segid {segment} and resid {resid}"))
            assert len(
                u_pdb.select_atoms(
                    f"segid {segment} and resid {resid} and name CA"
                )
            ) == len(
                u_cif.select_atoms(
                    f"segid {segment} and resid {resid} and name CA"
                )
            )


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
