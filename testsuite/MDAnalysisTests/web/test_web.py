import MDAnalysis as mda
import filecmp

import pytest
import requests

working_PDB_ID = '1DPX' # egg white lysozyme 

@pytest.fixture(scope="function")
def shared_cache_directory(tmp_path_factory):
    """Create Cache Directory"""   
    return tmp_path_factory.mktemp("cache") 


class Test_PDB_Downloader_BaseFunctionality():
    """Test Public API of Pdb_Downloader()"""
    
    valid_file_formats = ("pdb.gz", "pdb")

    def test_base_functionality(self):
        """Test file_formats and keywords in convert_to_universe()"""

        universe_list = []
        for file_format in self.valid_file_formats:
            downloader =  mda.web.PDBDownloader(PDB_ID=working_PDB_ID, file_format=file_format).download()
            universe_list.append(downloader.convert_to_universe(in_memory=True))

        assert all(isinstance(universe, mda.Universe) for universe in universe_list)

    def test_timeout(self):
        """Checks requests timeout Exception for fetch_pdb and PdbDownloader"""
        with pytest.raises(requests.exceptions.ConnectTimeout):
            mda.web.PDBDownloader(PDB_ID=working_PDB_ID).download(timeout=0.0000001)

    def test_invalid_id(self):
        """Test invalid id for PdbDownloader"""

        with pytest.raises(mda.web.downloaders.FileDownloadPDBError):
            mda.web.PDBDownloader(PDB_ID='BananaBoat').download()       

class Test_Pdb_Downloader_Cache():
    ### Cache Test underneath
    file_format = "pdb.gz" # default format for PdbDownloader

    def test_file_format(self):
        assert self.file_format == "pdb.gz"

    def test_file_name_cache(self, shared_cache_directory):
        """Test that cache is saved as ID.file_format"""
        downloader =  mda.web.PDBDownloader(PDB_ID=working_PDB_ID,)
        
        downloader.download(cache_path=shared_cache_directory)

        assert shared_cache_directory / f"{working_PDB_ID}.{self.file_format}"

    def test_loading_cache(self, shared_cache_directory):
        """Test that PdbDownloader.download() is reading cache"""

        downloader1 =  mda.web.PDBDownloader(PDB_ID=working_PDB_ID,
                                            file_format=self.file_format)

        downloader1.download(cache_path=shared_cache_directory)
        f1  = shared_cache_directory / f"{working_PDB_ID}.{self.file_format}"

        downloader2 =  mda.web.PDBDownloader(PDB_ID=working_PDB_ID,
                                            file_format=self.file_format)
        downloader2.download(cache_path=shared_cache_directory)
  
        f2  = shared_cache_directory / f"{working_PDB_ID}.{self.file_format}"
        
        assert filecmp.cmp(f1, f2, shallow=False) == True

class Test_Fetch_Pdb():
    """Tests Public API of MDAnalysis.fetch_pdb()"""

    def test_base_functionality(self):
        """Check default parameter of mda.fetch_pdb"""
        assert isinstance(mda.fetch_pdb(PDB_ID=working_PDB_ID), mda.Universe)

    def test_timeout(self):
        """Checks requests timeout Exception for fetch_pdb and PdbDownloader"""
        with pytest.raises(requests.exceptions.ConnectTimeout):
            mda.fetch_pdb(PDB_ID=working_PDB_ID, timeout=0.00000000001)

    def test_universe_keywords(self):
        """Check if keywords could be passed to Universe constructors"""
        assert isinstance(mda.fetch_pdb(working_PDB_ID, download_path=None, timeout=None, in_memory=True),
                          mda.Universe)








