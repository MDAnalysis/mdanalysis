from numpy.testing import assert_, assert_raises

import MDAnalysis as mda
from MDAnalysis.coordinates.base import Writer



class TestWriterCreation(object):
    class MagicWriter(Writer):
        # this writer does the 'magic' format
        format = 'MAGIC'

        def __init__(self, filename, n_atoms=None):
            self.filename = filename
            self.n_atoms = n_atoms

    class MultiMagicWriter(MagicWriter):
        # this writer does the 'magic' and 'magic2' formats
        # but *only* supports multiframe writing.
        format = ['MAGIC', 'MAGIC2']
        multiframe = True
        singleframe = False

    def test_default_multiframe(self):
        assert_(isinstance(mda.Writer('this.magic'), self.MultiMagicWriter))

    def test_singleframe(self):
        # check that singleframe=False has been respected
        assert_(isinstance(mda.Writer('this.magic', multiframe=False), self.MagicWriter))

    def test_multiframe_magic2(self):
        # this will work as we go for multiframe
        assert_(isinstance(mda.Writer('that.magic2'), self.MultiMagicWriter))

    def test_singleframe_magic2(self):
        # this should fail, there isn't a singleframe writer for magic2
        assert_raises(TypeError, 
                      mda.Writer, 'that.magic2', multiframe=False)
