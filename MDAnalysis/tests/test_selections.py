import MDAnalysis
import MDAnalysis.core.Selection
from MDAnalysis.tests.datafiles import PSF,DCD

from numpy.testing import *
from numpy import array, float32
from nose.plugins.attrib import attr


class TestSelections(TestCase):
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = MDAnalysis.Universe(PSF, DCD)

    def tearDown(self):
        del self.universe

    def test_segid(self):
        sel = self.universe.selectAtoms('segid 4AKE')
        assert_equal(sel.numberOfAtoms(), 3341, "failed to select segment 4AKE")
        assert_equal(sel._atoms, self.universe.s4AKE._atoms, 
                     "selected segment 4AKE is not the same as auto-generated segment s4AKE")
    def test_protein(self):
        sel = self.universe.selectAtoms('protein')
        assert_equal(sel.numberOfAtoms(), 3341, "failed to select protein")
        assert_equal(sel._atoms, self.universe.s4AKE._atoms, 
                     "selected protein is not the same as auto-generated protein segment s4AKE")

    def test_backbone(self):
        sel = self.universe.selectAtoms('backbone')
        assert_equal(sel.numberOfAtoms(), 855)

    def test_resid_single(self):
        sel = self.universe.selectAtoms('resid 100')
        assert_equal(sel.numberOfAtoms(), 7)
        assert_equal(sel.resnames(), ['GLY'])        

    def test_resid_range(self):
        sel = self.universe.selectAtoms('resid 100:105')
        assert_equal(sel.numberOfAtoms(), 89)
        assert_equal(sel.resnames(),  ['GLY', 'ILE', 'ASN', 'VAL', 'ASP', 'TYR'])        

    def test_resname(self):
        sel = self.universe.selectAtoms('resname LEU')
        assert_equal(sel.numberOfAtoms(), 304, "Failed to find all 'resname LEU' atoms.")
        assert_equal(sel.numberOfResidues(), 16, "Failed to find all 'resname LEU' residues.")
        assert_equal(sel._atoms, self.universe.s4AKE.LEU._atoms,
                     "selected 'resname LEU' atoms are not the same as aut-generated s4AKE.LEU")

    def test_name(self):
        sel = self.universe.selectAtoms('name CA')
        assert_equal(sel.numberOfAtoms(), 214)

    def test_atom(self):
        sel = self.universe.selectAtoms('atom 4AKE 100 CA')
        assert_equal(len(sel), 1)
        assert_equal(sel.resnames(), ['GLY'])
        assert_array_almost_equal(sel.coordinates(), 
                                  array([[ 20.38685226,  -3.44224262,  -5.92158318]], dtype=float32))

    def test_and(self):
        sel = self.universe.selectAtoms('resname GLY and resid 100')
        assert_equal(len(sel), 7)

    def test_or(self):
        sel = self.universe.selectAtoms('resname LYS or resname ARG')
        assert_equal(sel.numberOfResidues(), 31)

    def test_not(self):
        sel = self.universe.selectAtoms('not backbone')
        assert_equal(len(sel), 2486)


    # TODO:
    # add more test cases for around, point, prop, byres, bynum
    # and also for selection keywords such as 'nucleic' 

    def test_empty_selection(self):
        """Test empty selection: raises Excption but might be subject to change (see Issue 12)"""
        assert_raises(Exception, self.universe.selectAtoms, 'resname TRP')  # no Trp in AdK

    def test_parenthesized_expression(self):
        sel = self.universe.selectAtoms('( name CA or name CB ) and resname LEU')
        assert_equal(len(sel), 32)

    def test_no_space_around_parentheses(self):
        """Test that no space is needed around parentheses (Issue 43)."""
        # note: will currently be ERROR because it throws a ParseError
        sel = self.universe.selectAtoms('(name CA or name CB) and resname LEU')
        assert_equal(len(sel), 32)


    # TODO:
    # test for checking ordering and multiple comma-separated selections

    def test_concatenated_selection(self):
        E151 = self.universe.s4AKE.r151
        # note that this is not quite phi... HN should be C of prec. residue
        phi151 = E151.selectAtoms('name HN', 'name N', 'name CA', 'name CB')
        assert_equal(len(phi151), 4)
        assert_equal(phi151[0].name, 'HN', "wrong ordering in selection, should be HN-N-CA-CB") 
