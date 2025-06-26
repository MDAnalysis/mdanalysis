from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.log import ProgressBar

from pathlib import Path
from abc import ABC, abstractmethod

import io
import gzip

import requests


class FileDownloadPDBError(Exception):
    """
    Exception raised for errors in the file download process.

    Attributes:
        message -- explanation of the error
    """

    def __init__(
        self,
        message="There was an error downloading the file from the Protein Data Bank. PDB or format for PDB code may not be available.",
    ):

        super().__init__(message)


def _file_format_to_topology_string(file_extension):
    """Converts downloaded file name to a string usable by Universe(topology_format=)"""

    TOPOLOGY_FORMAT_CONVERTER = {
        "physical_file_extension": tuple(
            [["pdb.gz", "pdb"]],
        ),
        "topology_format": tuple(
            [
                "PDB",
            ],
        ),
    }

    for valid_file_extensions, valid_topology_string in zip(
        *TOPOLOGY_FORMAT_CONVERTER.values(), strict=True
    ):
        if file_extension in valid_file_extensions:
            return valid_topology_string


class BaseDownloader(ABC):
    """Abstract Base Class for all File-Based Downloaders. Not meant to be directly initalized!"""

    def __str__(self):
        file_set = True if self._file is not None else False
        return f"{self.__class__.__name__}: Filename: {self.file_name}, File set? {file_set}"

    def __del__(self):
        if self._file is not None:
            self._file.close()  # Ensure file-like object is closed in child classes

    def __init__(self, id, file_format):
        """Declares attributes meant to be manipulated by child classes"""

        # Attributes are meant to get updated by child instances
        self.id = id
        self.file_format = file_format
        self.file_name = f"{self.id}.{self.file_format}"

        self._file = None

    @abstractmethod
    def download(self):
        """This should be be implemented setting self._file to a file like object"""
        pass

    def convert_to_universe(self, **kwargs):
        """Converts Downloaded file to a Universe"""

        if self._file is None:
            raise RuntimeError(
                "File not set. Run download() to set file before convert_to_universe()"
            )

        try:
            return Universe(
                topology=self._file,
                topology_format=_file_format_to_topology_string(
                    self.file_format
                ),
                **kwargs,
            )
        finally:
            self._file.close()


class PDBDownloader(BaseDownloader):
    """Class to handle download PDBs from the Protein Data Bank"""

    def __init__(self, PDB_ID, file_format="pdb.gz"):

        super().__init__(PDB_ID, file_format)

        self._download = False

    def _open_file(self, cache_path):
        """This method either load/create cache or reserve a spot in memory to store topologt"""

        if cache_path is None:
            self._file = io.BytesIO()
            self._download = True

        else:
            cache_file_path = Path(cache_path) / self.file_name

            # Found Cache, so don't download anything and open existing file
            # Note this doesn't check the content of the file!
            if cache_file_path.exists() and cache_file_path.is_file():
                self._file = open(cache_file_path, "r")
                self._download = False

            else:  # No cache found, so create Cache
                self._file = open(cache_file_path, "wb")
                self._download = True

    def _requests_progress_bar(self, requests_response):
        """Puts a progress bar when writing content with a request object"""
        chunk_size = (
            1  # Files are so small that you can read them one byte at a time
        )

        with ProgressBar(
            total=len(requests_response.content),
            unit="B",
            unit_scale=True,
            desc=self.file_name,
        ) as pb:
            for byte in requests_response.iter_content(chunk_size=chunk_size):
                self._file.write(byte)
                pb.update(chunk_size)

    def download(self, cache_path=None, timeout=None, progress_bar=False):
        """Downloads files from the Protein Data Bank"""

        self._open_file(cache_path)

        # Downloads file into temporary file or cached it into a pernament file
        if self._download:
            try:
                r = requests.get(
                    f"https://files.rcsb.org/download/{self.id}.{self.file_format}",
                    timeout=timeout,
                    stream=progress_bar,
                )

                # Moves to except block if invalid PDB code was declared!
                r.raise_for_status()

                if progress_bar:
                    self._requests_progress_bar(r)
                else:
                    self._file.write(r.content)

                # Download Cleanup!
                # Universe doesn't accept binary stream only a text stream
                # so need to convert BinaryIO to StringIO first!
                self._file.seek(0)
                if isinstance(self._file, io.BytesIO):
                    if self.file_format == "pdb.gz":
                        opened_gzip = gzip.open(self._file).read()
                        self._file = io.StringIO(opened_gzip.decode("UTF-8"))
                    elif self.file_format == "pdb":
                        self._file = io.StringIO(
                            self._file.read().decode("UTF-8")
                        )

                # Cleanly saves files to disk for cache
                elif cache_path is not None:
                    self._file.close()

            except requests.HTTPError:
                raise FileDownloadPDBError

        return self
