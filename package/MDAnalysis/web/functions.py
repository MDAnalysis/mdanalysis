from .downloaders import PdbDownloader

def fetch_pdb(PDB_ID, file_format="pdb.gz", download_path=None, timeout=None, **kwargs):
    """Fetchs PDB from RCSB"""

    downloader = PdbDownloader(PDB_ID, file_format).download(download_path, timeout=timeout)

    return downloader.convert_to_universe(**kwargs)
