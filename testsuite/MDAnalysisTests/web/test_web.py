import MDAnalysis as mda
import pytest
import requests



working_PDB_ID =  '1DPX' # egg white lysozyme 
file_format="pdb"



def test_fetch_pdb_base_functionality():
    """Check default parameter of mda.fetch_pdb"""
    assert isinstance(mda.fetch_pdb(PDB_ID=working_PDB_ID), mda.Universe)


def test_fetch_pdb_timeout():
    """Checks requests timeout Exception for fetch_pdb and PdbDownloader"""
    with pytest.raises(requests.exceptions.ConnectTimeout):
        mda.fetch_pdb(PDB_ID=working_PDB_ID, timeout=0.000001)

def test_fetch_pdb_keywords():
    """Check if keywords could be passed to Universe consturctors"""
    assert isinstance(mda.fetch_pdb(working_PDB_ID, download_path=None, timeout=None, in_memory=True),
                      mda.Universe)


def test_file_namecache(tmp_path):
    """Test that cache is saved as ID.file_format"""
    test_directory = tmp_path
    
    downloader =  mda.web.PdbDownloader(PDB_ID=working_PDB_ID,
                                        file_format=file_format)

    downloader.download(cache_path=test_directory)

    assert test_directory / f"{working_PDB_ID}.{file_format}"


def test_invalid_id():
    """Test invalid id for PdbDownloader and fetch_pdb"""

    with pytest.raises(mda.web.downloaders.FileDownloadPDBError):
        mda.web.PdbDownloader(PDB_ID='BananaBoat').download()

    






