from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.log import ProgressBar

from pathlib import Path
from abc import ABC, abstractmethod
#from io import BytesIO

import tempfile
import requests

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

# Dict to convert file name arguments for Universe(topology_format="")
TOPOLOGY_FORMAT_CONVERTER = {
    "physical_file_extension": tuple([["pdb.gz", "pdb"]],
                            ) ,
    "topology_format": tuple(["PDB",],
                             )
    }

def _file_format_to_topology_string(file_extension):
    """Converts file names to a string usable by Universe(topology_format=)"""
    for valid_file_extensions, valid_topology_string in zip(*TOPOLOGY_FORMAT_CONVERTER.values(), strict=True):
        if file_extension in valid_file_extensions:
            return valid_topology_string

class BaseDownloader(ABC):
    """Abstract Base Class for all File-Based Downloaders. Not meant to be directly initalized!"""

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
            return Universe(self._file.name,
                            topology_format=_file_format_to_topology_string(self.file_format),
                            **kwargs)
        finally:
            self._file.close() 

def _requests_progress_bar(requests_response, file_name, file_writer, return_writer=False):
    """Puts a progress bar while writing a file_like object from the web"""
    chunk_size = 1 
    r = requests_response

    with ProgressBar(total=len(r.content), unit='B', unit_scale=True, desc=file_name) as pb:
        for i in r.iter_content(chunk_size=chunk_size):
            file_writer.write(i)
            pb.update(chunk_size)

    if return_writer:
        return file_writer

class PdbDownloader(BaseDownloader):
    """Class to handle download PDBs from the RCSB"""
     
    def __init__(self, PDB_ID, file_format="pdb.gz"):

        super().__init__(PDB_ID, file_format)

        self._download = False

    def _open_file(self, cache_path):
        """Sets private file-like attribute _file to a physical file on disk"""

        # create temporary file to save pdb file
        if cache_path is None: 
            self._file = tempfile.NamedTemporaryFile('wb')
            self._download = True

        # Create/Parse download cache
        else: 
            named_file_path = Path(cache_path) / f"{self.id}.{self.file_format}"
            
            # Found Cache, so don't download anything and open existing file
            if named_file_path.exists() and named_file_path.is_file():
                self._file = open(named_file_path, 'r')     
                self._download = False                 
            
            else: # No cache found, so create Cache
                self._file = open(named_file_path, 'wb')
                self._download = True


    def download(self, cache_path=None, timeout=None, progress_bar=False):
        """Downloads files from the RCSB"""

        # Sets self._file correctly
        self._open_file(cache_path)

        # Downloading file into temporary file or cached it into a pernament file
        if self._download:
            try:
                r = requests.get(f"https://files.rcsb.org/download/{self.id}.{self.file_format}",
                                timeout=timeout, stream=progress_bar)
                
                r.raise_for_status()
                if progress_bar:
                    _requests_progress_bar(r, f"{self.id}.{self.file_format}", self._file)
                else:
                    self._file.write(r.content)
            except requests.HTTPError:
                raise FileDownloadPDBError
            finally:
                # Closes File safely if saving to cache
                if cache_path is not None:
                    self._file.close()
                    
        return self

