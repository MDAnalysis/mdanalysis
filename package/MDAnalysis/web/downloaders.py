from pathlib import Path
from abc import ABC, abstractmethod
import tempfile

import requests
from MDAnalysis.core.universe import Universe

class FileDownloadPDBError(Exception):
    """
    Exception raised for errors in the file download process.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, 
        message="There was an error downloading the file from the Protein Data Bank. PDB or format for PDB code may not be available.",
        ):

        super().__init__(message)


class BaseDownloader(ABC):
    """Abstract Base Class for all Downloaders. Not meant to be directly initalized!"""

    def __str__(self):
        return f"Metadata: id={self.id}, file_format={self.file_format}, "
    
    def __del__(self):
        if self._file is not None:
            self._file.close() # Ensure file-like object is closed in child classes

    def __init__(self, id, file_format):
        """Declares attributes meant to be manipulated by child classes"""

        # Attributes are meant to get updated by child instances 
        self.id = id
        self.file_format = file_format

        self._file = None

    @abstractmethod
    def download():
        """This should be be implelemented setting self._file to a file like object"""
        pass

    def convert_to_universe(self, **kwargs):
        """Converts Downloaded file to a Universe"""

        if self._file is None:
            raise RuntimeError("File not set. Run download() to set file before convert_to_universe()")

        try:
            return Universe(self._file.name, topology_format=self.file_format.upper(), **kwargs)
        finally:
            self._file.close() 

class PdbDownloader(BaseDownloader):
    """Class to handle download PDBs from the RCSB"""
     
    def __init__(self, PDB_ID, file_format="pdb"):
        super().__init__(PDB_ID, file_format)

        self._download = False

    def _open_file(self, cache_path):
        """Sets private file-like attribute _file to a physical file on disk"""

        # create temporary file to save pdb file
        if cache_path is None: 
            self._file = tempfile.NamedTemporaryFile(mode='wb')
            self._download = True

        # Create/Parse download cache
        else: 
            named_file_path = Path(cache_path) /  f"{self.id}.{self.file_format}"
            
            # Found Cache, so don't download anything and open existing file
            if named_file_path.exists() and named_file_path.is_file():
                self._file = open(named_file_path, 'r')     
                self._download = False                 
            
            else: # No cache found, so create Cache
                self._file = open(named_file_path, 'wb')
                self._download = True


    def download(self, cache_path=None, timeout=None):
        """Downloads files from the RCSB"""

        # Sets self._file correctly
        self._open_file(cache_path)

        # Downloading file into temporary file or cached it into a pernament file
        if self._download:
            try:
                r = requests.get(f"https://files.rcsb.org/download/{self.id}.{self.file_format}",
                                timeout=timeout)
                r.raise_for_status()
                self._file.write(r.content)
            except requests.HTTPError:
                # This also deletes the undownloaded file since write() hasn't been called yet
                raise FileDownloadPDBError
            finally:
                # Closes File safely if saving to cache
                if cache_path is not None:
                    self._file.close()
                    
        return self
