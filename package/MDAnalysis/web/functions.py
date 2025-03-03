from .downloaders import PdbDownloader


def fetch_pdb(PDB_ID, download_path=None, timeout=30, **kwargs):
    """Fetchs PDB from RCSB"""

    downloader = PdbDownloader(PDB_ID).download(download_path,
                                                timeout=timeout)

    return downloader.convert_to_universe(**kwargs)
