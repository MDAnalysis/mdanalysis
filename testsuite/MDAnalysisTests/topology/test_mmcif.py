import MDAnalysis as mda
import pytest
from pathlib import Path
from io import StringIO
import gzip
from numpy.testing import assert_allclose, assert_equal
from MDAnalysis.lib import util
from MDAnalysis.coordinates.MMCIF import HAS_GEMMI

from MDAnalysisTests.datafiles import MMCIF as MMCIF_FOLDER
from MDAnalysisTests.datafiles import CONECT, PDB_multiframe


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
def test_types_from_elements():
    # as in PDBParser, atom types come from the element column, and
    # elements are validated against known symbols
    u = mda.Universe(f"{MMCIF_FOLDER}/1YJP.cif")
    assert_equal(u.atoms.types, u.atoms.elements)
    assert set(u.atoms.elements) == {"C", "N", "O"}


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_unknown_elements():
    # 5PTI is a joint X-ray/neutron structure: 229 deuterium (D) atoms and
    # one unknown (X) atom, neither of which MDAnalysis knows as an element
    with pytest.warns(UserWarning, match="Unknown element") as record:
        u = mda.Universe(f"{MMCIF_FOLDER}/5PTI.cif.gz")
    assert {str(w.message).split()[2] for w in record} >= {"D", "X"}

    unknown = u.atoms[u.atoms.elements == ""]
    assert len(unknown) == 230
    # types keep the raw element column
    assert set(unknown.types) == {"D", "X"}
    # both unknowns intern to a single "" entry, so lookups find all of them
    elements = u._topology.elements
    assert len(set(elements.name_lookup)) == len(elements.name_lookup)
    assert (elements.nmidx == elements.namedict[""]).sum() == 230


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_altlocs():
    # gemmi stores absent altlocs as "\0", which must not leak into the
    # topology (it ends up in written PDB files otherwise)
    u = mda.Universe(f"{MMCIF_FOLDER}/3PWP.cif.gz")
    assert not any("\x00" in altloc for altloc in u.atoms.altLocs)
    with_altloc = u.atoms[[altloc != "" for altloc in u.atoms.altLocs]]
    assert len(with_altloc) == 12
    assert set(with_altloc.altLocs) == {"A", "B"}
    assert set(with_altloc.resnames) == {"GLN"}


def test_no_gemmi_raises(monkeypatch):
    monkeypatch.setattr(mda.topology.MMCIFParser, "HAS_GEMMI", False)
    with pytest.raises(ImportError, match="please install gemmi"):
        mda.topology.MMCIFParser.MMCIFParser(f"{MMCIF_FOLDER}/1YJP.cif")


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_wrong_format():
    with pytest.raises(ValueError):
        mda.Universe(f"{MMCIF_FOLDER}/1YJP_invalid.cif.gz")


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_unparseable_content():
    # content that fails both CIF and PDB parsing in gemmi must re-raise
    # the original CIF error
    garbage = util.NamedStream(
        StringIO("data_\n_nonsense.field ??? ][\n"), "garbage.cif"
    )
    with pytest.raises(ValueError, match="parse error"):
        mda.Universe(garbage)


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_blockless_content():
    # a comment-only file parses as a CIF document with no data blocks
    empty = util.NamedStream(StringIO("# no data blocks here\n"), "empty.cif")
    with pytest.raises(ValueError, match="No data blocks found"):
        mda.Universe(empty)


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_multimodel_warning_msg():
    with pytest.warns(
        UserWarning,
        match=(
            r"MMCIF model .+ contains .+ different models, but only the "
            r"first one will be used to assign the topology"
        ),
    ):
        mda.topology.MMCIFParser.MMCIFParser(
            f"{MMCIF_FOLDER}/multimodel_warning.cif.gz"
        ).parse()


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize(
    "filename,fmt",
    [
        (f"{MMCIF_FOLDER}/1BD2_short.cif.gz", None),
        (Path(f"{MMCIF_FOLDER}/1BD2_short.cif.gz"), None),
        (
            StringIO(util.anyopen(f"{MMCIF_FOLDER}/1BD2_short.cif.gz").read()),
            "CIF",
        ),
        (gzip.open(f"{MMCIF_FOLDER}/1BD2_short.cif.gz"), "CIF"),
        (
            util.NamedStream(
                StringIO(
                    util.anyopen(f"{MMCIF_FOLDER}/1BD2_short.cif.gz").read()
                ),
                "some_name.cif",
            ),
            "CIF",
        ),
        (f"{MMCIF_FOLDER}/1BD2_short.pdb.gz", None),
        (Path(f"{MMCIF_FOLDER}/1BD2_short.pdb.gz"), None),
        (
            StringIO(util.anyopen(f"{MMCIF_FOLDER}/1BD2_short.pdb.gz").read()),
            "CIF",
        ),
        (
            util.anyopen(f"{MMCIF_FOLDER}/1BD2_short.pdb.gz"),
            "CIF",
        ),
        (
            util.NamedStream(
                StringIO(
                    util.anyopen(f"{MMCIF_FOLDER}/1BD2_short.pdb.gz").read()
                ),
                "some_name.pdb",
            ),
            "CIF",
        ),
    ],
)
def test_input_methods(filename, fmt):
    mda.Universe(filename, topology_format=fmt)


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize("basename", ["1BD2_short", "1BD2", "3PWP", "3KPR"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pdb_gemmi_matches_mmcif(basename):
    u_pdb = mda.Universe(
        f"{MMCIF_FOLDER}/{basename}.pdb.gz", format="pdb_gemmi"
    )
    u_cif = mda.Universe(f"{MMCIF_FOLDER}/{basename}.cif.gz")
    assert isinstance(u_pdb.trajectory, mda.coordinates.MMCIF.MMCIFReader)
    for attr in (
        "names",
        "resnames",
        "resids",
        "elements",
        "chainIDs",
        "altLocs",
        "record_types",
        "icodes",
    ):
        assert_equal(getattr(u_pdb.atoms, attr), getattr(u_cif.atoms, attr))
    assert_allclose(u_pdb.atoms.positions, u_cif.atoms.positions, atol=1e-3)


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pdb_gemmi_hybrid36_resids():
    # residue numbers above 9999 are written in hybrid-36: A000 is 10000
    record = (
        "HETATM{:5d}  O   HOH W{:>4}    {:8.3f}   0.000   0.000"
        "  1.00  0.00           O"
    )
    lines = [
        record.format(i, resid, i)
        for i, resid in enumerate(["9999", "A000", "A001"], start=1)
    ]
    pdb = util.NamedStream(StringIO("\n".join(lines) + "\nEND\n"), "water.pdb")
    u = mda.Universe(pdb, format="pdb_gemmi")
    assert_equal(u.residues.resids, [9999, 10000, 10001])


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.parametrize("filename", [CONECT, PDB_multiframe])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pdb_gemmi_conect_bonds(filename):
    u_gemmi = mda.Universe(filename, format="pdb_gemmi")
    u_pdb = mda.Universe(filename)
    # bonds are matched by atom serial since gemmi may reorder atoms
    gemmi_bonds = {tuple(sorted(b.atoms.ids)) for b in u_gemmi.bonds}
    pdb_bonds = {tuple(sorted(b.atoms.ids)) for b in u_pdb.bonds}
    assert len(gemmi_bonds) > 0
    assert gemmi_bonds == pdb_bonds


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_mmcif_no_bonds():
    u = mda.Universe(f"{MMCIF_FOLDER}/1BD2.cif.gz")
    assert not hasattr(u.atoms, "bonds")


@pytest.mark.skipif(not HAS_GEMMI, reason="gemmi not installed")
def test_pdb_gemmi_conect_duplicate_serials():
    record = (
        "ATOM  {:5d}  C   UNK A   1    {:8.3f}   0.000   0.000"
        "  1.00  0.00           C"
    )
    lines = [record.format(serial, x) for x, serial in enumerate([1, 1, 2])]
    lines.append("CONECT    1    2")
    pdb = util.NamedStream(StringIO("\n".join(lines) + "\nEND\n"), "dup.pdb")
    with pytest.warns(UserWarning, match="bonds will not be parsed"):
        u = mda.Universe(pdb, format="pdb_gemmi")
    assert not hasattr(u.atoms, "bonds")
