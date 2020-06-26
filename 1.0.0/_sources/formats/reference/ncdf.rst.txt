.. -*- coding: utf-8 -*-
.. _NCDF-format:

===========================================
NCDF, NC (AMBER NetCDF trajectory)
===========================================

.. include:: classes/NCDF.txt

`AMBER binary trajectories`_ are automatically recognised by the file extension ".ncdf". The NCDF module uses :mod:`scipy.io.netcdf` and therefore :mod:`scipy` must be installed. 

.. _`AMBER binary trajectories`: https://ambermd.org/FileFormats.php#netcdf

Reading in
==========

Units are assumed to be the following default AMBER units:

    * length: Angstrom
    * time: ps

Currently, if other units are detected, MDAnalysis will raise a :exc:`NotImplementedError`. 


Writing out
===========

NCDF files are always written out in ångström and picoseconds. 

Although scale_factors can be read from NCDF files, they are not kept or used when writing NCDF files out.


**Writing with the netCDF4 module and potential issues**

    Although :mod:`scipy.io.netcdf` is very fast at reading NetCDF files, it is slow at writing them out. The netCDF4_ package is fast at writing (but slow at reading). This requires the compiled netcdf library to be installed. MDAnalysis tries to use :mod:`netCDF4` for writing if it is available, but will fall back to :mod:`scipy.io.netcdf` if it is not.

    **AMBER users** might have a hard time getting netCDF4 to work with a
    conda-based installation (as discussed in `Issue #506`_) because of the way
    that AMBER itself handles netcdf. In this scenario, MDAnalysis will simply switch to the scipy package. If you encounter this error and wish to use the the faster :mod:`netCDF4` writer, the only solution is to unload the AMBER environment.


    .. _netCDF4: https://unidata.github.io/netcdf4-python/
    .. _`Issue #506`:
       https://github.com/MDAnalysis/mdanalysis/issues/506#issuecomment-225081416

