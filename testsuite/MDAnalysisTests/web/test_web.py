import MDAnalysis as mda
import filecmp


import pytest
import requests



working_PDB_ID = '1DPX' # egg white lysozyme 
file_format = "pdb"


@pytest.fixture(scope="function")
def shared_cache_directory(tmp_path_factory):
    """Create Cache Directory"""   
    return tmp_path_factory.mktemp("cache") 


class Test_Pdb_Downloader_BaseFunctionality():
    """Test Public API of Pdb_Downloader()"""
   
    def test_base_functionality(self):
        """Test keywords in convert_to_universe()"""
        downloader =  mda.web.PdbDownloader(PDB_ID=working_PDB_ID).download()
    
        assert(isinstance(downloader.convert_to_universe(in_memory=True), mda.Universe))

    def test_timeout(self):
        """Checks requests timeout Exception for fetch_pdb and PdbDownloader"""
        with pytest.raises(requests.exceptions.ConnectTimeout):
            mda.web.PdbDownloader(PDB_ID=working_PDB_ID).download(timeout=0.00000001)

    def test_invalid_id(self):
        """Test invalid id for PdbDownloader"""

        with pytest.raises(mda.web.downloaders.FileDownloadPDBError):
            mda.web.PdbDownloader(PDB_ID='BananaBoat').download()

class Test_Pdb_Downloader_Cache():
    ### Cache Test underneath
    def test_file_name_cache(self, shared_cache_directory):
        """Test that cache is saved as ID.file_format"""
        downloader =  mda.web.PdbDownloader(PDB_ID=working_PDB_ID)
        
        downloader.download(cache_path=shared_cache_directory)

        assert shared_cache_directory / f"{working_PDB_ID}.{file_format}"

    def test_loading_cache(self, shared_cache_directory):
        """Test that PdbDownloader.download() is reading cache"""

        downloader1 =  mda.web.PdbDownloader(PDB_ID=working_PDB_ID,
                                            file_format=file_format)

        downloader1.download(cache_path=shared_cache_directory)
        f1  = shared_cache_directory / f"{working_PDB_ID}.{file_format}"

        downloader2 =  mda.web.PdbDownloader(PDB_ID=working_PDB_ID,
                                            file_format=file_format)
        downloader2.download(cache_path=shared_cache_directory)
  
        f2  = shared_cache_directory / f"{working_PDB_ID}.{file_format}"
        
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








