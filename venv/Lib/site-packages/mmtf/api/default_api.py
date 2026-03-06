import gzip

from mmtf.utils.constants import BASE_URL, BASE_URL_REDUCED
try:
    import urllib2
    from StringIO import StringIO
except:
    import urllib.request as urllib2
    from io import BytesIO as StringIO
import msgpack
from .mmtf_reader import MMTFDecoder
from .mmtf_writer import MMTFEncoder,TemplateEncoder


def _internet_on(address):
    """
    Check to see if the internet is on by pinging a set address.
    :param address: the IP or address to hit
    :return: a boolean - true if can be reached, false if not.
    """
    try:
        urllib2.urlopen(address, timeout=1)
        return True
    except urllib2.URLError as err:
        return False

def write_mmtf(file_path, input_data, input_function):
    """API function to write data as MMTF to a file

    :param file_path the path of the file to write
    :param input_data the input data in any user format
    :param input_function a function to converte input_data to an output format. Must contain all methods in TemplateEncoder
    """
    mmtf_encoder = MMTFEncoder()
    pass_data_on(input_data, input_function, mmtf_encoder)
    mmtf_encoder.write_file(file_path)

def pass_data_on(input_data, input_function, output_data):
    """Helper to pass data from one data structure to another.

    :param input_data the input data in any user format
    :param input_function a function to converte input_data to an output format. Must contain all methods in TemplateEncoder.
    :param output_data a data holder for the data to be put into."""
    input_function(input_data, output_data)
    return output_data


def get_raw_data_from_url(pdb_id, reduced=False):
    """" Get the msgpack unpacked data given a PDB id.

    :param pdb_id: the input PDB id
    :return the unpacked data (a dict) """
    url = get_url(pdb_id,reduced)
    request = urllib2.Request(url)
    request.add_header('Accept-encoding', 'gzip')
    response = urllib2.urlopen(request)
    if response.info().get('Content-Encoding') == 'gzip':
        data = ungzip_data(response.read())
    else:
        data = response.read()
    return _unpack(data)

def get_url(pdb_id,reduced=False):
    """Get the URL for the data for a given PDB id.

    :param pdb_id: the input PDB id
    :return the URL for this PDB id"""
    if reduced:
        return BASE_URL_REDUCED + pdb_id
    else:
        return BASE_URL + pdb_id

def _unpack(data):
    out_data = msgpack.unpackb(data.read(), raw=False)
    return out_data


def fetch(pdb_id):
    """Return a decoded API to the data from a PDB id.

    :param pdb_id: the input PDB id
    :return an API to decoded data """
    decoder = MMTFDecoder()
    decoder.decode_data(get_raw_data_from_url(pdb_id))
    return decoder

def parse(file_path):
    """Return a decoded API to the data from a file path.

    :param file_path: the input file path. Data is not entropy compressed (e.g. gzip)
    :return an API to decoded data """
    newDecoder = MMTFDecoder()
    with open(file_path, "rb") as fh:
        newDecoder.decode_data(_unpack(fh))
    return newDecoder


def parse_gzip(file_path):
    """Return a decoded API to the data from a file path. File is gzip compressed.
    :param file_path: the input file path. Data is gzip compressed.
    :return an API to decoded data"""
    newDecoder = MMTFDecoder()
    newDecoder.decode_data(_unpack(gzip.open(file_path, "rb")))
    return newDecoder


def ungzip_data(input_data):
    """Return a string of data after gzip decoding

    :param the input gziped data
    :return  the gzip decoded data"""
    buf = StringIO(input_data)
    f = gzip.GzipFile(fileobj=buf)
    return f
