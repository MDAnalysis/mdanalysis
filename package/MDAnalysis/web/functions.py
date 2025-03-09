from .downloaders import PDBDownloader

def fetch_pdb(PDB_ID, file_format="pdb.gz", download_path=None, timeout=None, progress_bar=False, **kwargs):
    """Fetchs PDB from RCSB"""

    downloader = PDBDownloader(PDB_ID, file_format).download(download_path, timeout, progress_bar)

    return downloader.convert_to_universe(**kwargs)
