

""""""# start delvewheel patch
def _delvewheel_init_patch_1_3_3():
    import os
    import sys
    libs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'netCDF4.libs'))
    is_pyinstaller = getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')
    if not is_pyinstaller or os.path.isdir(libs_dir):
        os.add_dll_directory(libs_dir)


_delvewheel_init_patch_1_3_3()
del _delvewheel_init_patch_1_3_3
# end delvewheel patch

# init for netCDF4. package
# Docstring comes from extension module _netCDF4.
from ._netCDF4 import *
# Need explicit imports for names beginning with underscores
from ._netCDF4 import __doc__
from ._netCDF4 import (__version__, __netcdf4libversion__, __hdf5libversion__,
                       __has_rename_grp__, __has_nc_inq_path__,
                       __has_nc_inq_format_extended__, __has_nc_open_mem__,
                       __has_nc_create_mem__, __has_cdf5_format__,
                       __has_parallel4_support__, __has_pnetcdf_support__,
                       __has_quantization_support__, __has_zstandard_support__,
                       __has_bzip2_support__, __has_blosc_support__, __has_szip_support__,
                       __has_set_alignment__)
import os
__all__ =\
['Dataset','Variable','Dimension','Group','MFDataset','MFTime','CompoundType','VLType','date2num','num2date','date2index','stringtochar','chartostring','stringtoarr','getlibversion','EnumType','get_chunk_cache','set_chunk_cache','set_alignment','get_alignment']
# if HDF5_PLUGIN_PATH not set, point to package path if plugins live there
pluginpath = os.path.join(__path__[0],'plugins')
if 'HDF5_PLUGIN_PATH' not in os.environ and\
    (os.path.exists(os.path.join(pluginpath,'lib__nczhdf5filters.so')) or\
     os.path.exists(os.path.join(pluginpath,'lib__nczhdf5filters.dylib'))):
    os.environ['HDF5_PLUGIN_PATH']=pluginpath
