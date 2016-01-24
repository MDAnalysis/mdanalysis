# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

# TPR parser and tpr support module
# Copyright (c) 2011 Zhuyi Xue
# Released under the  GNU Public Licence, v2


"""
Gromacs portable run input TPR format parser
============================================

The :mod:`~MDAnalysis.topology.TPRParser` module allows reading of a
Gromacs_ portable run input file (a `TPR file`_). Because
the file format of the TPR file is changing rapidly, not all versions
are currently supported. The known working versions and the
approximate Gromacs release numbers are listed in the table
:ref:`TPR format versions <TPR-format-table>`.

.. _`TPR-format-table`:

.. table:: TPR format versions and generations read by :func:`MDAnalysis.topology.TPRParser.parse`.

   ========== ============== ==================== =====
   TPX format TPX generation Gromacs release      read
   ========== ============== ==================== =====
   ??         ??             3.3, 3.3.1           no

   58         17             4.0, 4.0.2, 4.0.3,   yes
                             4.0.4, 4.0.5, 4.0.6,
                             4.0.7

   73         23             4.5.0, 4.5.1, 4.5.2, yes
                             4.5.3, 4.5.4, 4.5.5

   83         24             4.6, 4.6.1           yes

   100        26             5.0, 5.0.1, 5.0.2,   yes
                             5.0.3,5.0.4, 5.0.5

   103        26             5.1                  yes
   ========== ============== ==================== =====

For further discussion and notes see `Issue 2`_. Also add a comment to
`Issue 2`_ if a new or different TPR file format version should be
supported.

Bonded interactions available in Gromacs are described in table 5.5 of the
`Gromacs manual`_. The following ones are used to build the topology (see
`Issue 463`_):

* bonds: regular bonds (type 1), G96 bonds (type 2), Morse (type 3),
  cubic bonds (type 4), connections (type 5), harmonic potentials (type 6),
  FENE bonds (type 7), restraint potentials (type 10),
  tabulated potential with exclusion/connection (type 8),
  tabulated potential without exclusion/connection (type 9), constraints with
  exclusion/connection (type 1), constraints without exclusion/connection (type
  2)
* angles: regular angles (type 1), G96 angles (type 2), cross bond-bond
  (type3), cross-bond-angle (type 4), Urey-Bradley (type 5), quartic angles
  (type 6), restricted bending potential (type 10), tabulated angles (type 8)
* dihedrals: proper dihedrals (type 1 and type 9), Ryckaert-Bellemans dihedrals
  (type 3), Fourier dihedrals (type 5), restricted dihedrals (type 10),
  combined bending-torsion potentials (type 11), tabulated dihedral (type 8)
* impropers: improper dihedrals (type 2), periodic improper dihedrals (type 4)


Classes
---------

.. autoclass:: TPRParser
   :members:
   :inherited-members:

.. SeeAlso:: :mod:`MDAnalysis.topology.tpr`

Development notes
-----------------

The TPR reader is a pure-python implementation of a basic TPR
parser. Currently the following sections of the topology are parsed:

* Atoms: number, name, type, resname, resid, segid, mass, charge,
  [residue, segment, radius, bfactor, resnum]
* Bonds
* Angels
* Dihedrals
* Impropers

This tpr parser is written according to the following files

- :file:`{gromacs_dir}/src/kernel/gmxdump.c`
- :file:`{gromacs_dir}/src/gmxlib/tpxio.c` (the most important one)
- :file:`{gromacs_dir}/src/gmxlib/gmxfio_rw.c`
- :file:`{gromacs_dir}/src/gmxlib/gmxfio_xdr.c`
- :file:`{gromacs_dir}/include/gmxfiofio.h`

or their equivalent in more recent versions of Gromacs.

The function :func:`read_tpxheader` is based on the
`TPRReaderDevelopment`_ notes.  Functions with names starting with
``read_`` or ``do_`` are trying to be similar to those in
:file:`gmxdump.c` or :file:`tpxio.c`, those with ``extract_`` are new.

Wherever ``fver_err(fver)`` is used, it means the tpx version problem
has not been solved. Versions prior to Gromacs 4.0.x are not supported.

.. Links
.. _Gromacs: http://www.gromacs.org
.. _`Gromacs manual`: http://manual.gromacs.org/documentation/5.1/manual-5.1.pdf
.. _TPR file: http://manual.gromacs.org/current/online/tpr.html
.. _`Issue 2`: https://github.com/MDAnalysis/mdanalysis/issues/2
.. _`Issue 463`: https://github.com/MDAnalysis/mdanalysis/pull/463
.. _TPRReaderDevelopment: https://github.com/MDAnalysis/mdanalysis/wiki/TPRReaderDevelopment
"""
from __future__ import absolute_import
__author__ = "Zhuyi Xue"
__copyright__ = "GNU Public Licence, v2"

import xdrlib

from ..lib.util import anyopen
from .tpr import utils as U
from .base import TopologyReader

import logging
logger = logging.getLogger("MDAnalysis.topology.TPRparser")


class TPRParser(TopologyReader):
    """Read topology information from a Gromacs_ TPR_ file.

    .. SeeAlso:: :mod:`MDAnalysis.topology.TPR`

    .. _Gromacs: http://www.gromacs.org
    .. _TPR file: http://manual.gromacs.org/current/online/tpr.html
    """
    format = 'TPR'

    def parse(self):
        """Parse a Gromacs TPR file into a MDAnalysis internal topology structure.

        :Returns: ``structure`` dict
        """
        #ndo_int = U.ndo_int
        ndo_real = U.ndo_real
        #ndo_rvec = U.ndo_rvec
        #ndo_ivec = U.ndo_ivec

        tprf = anyopen(self.filename, mode='rb').read()
        data = xdrlib.Unpacker(tprf)
        try:
            th = U.read_tpxheader(data)                    # tpxheader
        except EOFError:
            msg = "{0}: Invalid tpr file or cannot be recognized".format(self.filename)
            logger.critical(msg)
            raise IOError(msg)

        self._log_header(th)

        V = th.fver                                    # since it's used very often

        state_ngtc = th.ngtc         # done init_state() in src/gmxlib/tpxio.c
        if th.bBox:
            U.extract_box_info(data, V)

        if state_ngtc > 0 and V >= 28:
            if V < 69:                      # redundancy due to  different versions
                ndo_real(data, state_ngtc)
            ndo_real(data, state_ngtc)        # relevant to Berendsen tcoupl_lambda

        if V < 26:
            U.fver_err(V)

        if th.bTop:
            tpr_top = U.do_mtop(data, V, self._u)
            structure = {
                'atoms': tpr_top.atoms,
                'bonds': tpr_top.bonds,
                'angles': tpr_top.angles,
                'dihedrals': tpr_top.dihe,
                'impropers': tpr_top.impr
            }
        else:
            msg = "{0}: No topology found in tpr file".format(self.filename)
            logger.critical(msg)
            raise IOError(msg)

        return structure

    # THE FOLLOWING CODE IS WORKING FOR TPX VERSION 58, BUT SINCE THESE INFO IS
    # NOT INTERESTED, SO IT'S NOT COVERED IN ALL VERSIONS. PARSING STOPS HERE.

    # if th.bX:
    #     ndo_rvec(data, th.natoms)

    # if th.bV:
    #     ndo_rvec(data, th.natoms)

    # if th.bF:
    #     ndo_rvec(data, th.natoms)

    # not useful at the moment
    # ePBC = -1;
    # bPeriodicMols = False
    # if th.bIr:
    #     # update
    #     data.unpack_int()                                # ePBC
    #     data.unpack_bool()                               # bPeriodicMols
    #     # 17 < 23. and ir (ir is from the c code, seems not apply here
    #     if th.fgen < setting.tpx_generation:
    #         # a crazily long (670 lines) function in c, slightly better here
    #         # (240 lines), so put it in setting.py
    #         utils.do_inputrec(data)

    def _log_header(self, th):
        logger.info("Gromacs version   : {0}".format(th.ver_str))
        logger.info("tpx version       : {0}".format(th.fver))
        logger.info("tpx generation    : {0}".format(th.fgen))
        logger.info("tpx number        : {0}".format(th.number))
        logger.info("tpx precision     : {0}".format(th.precision))
        logger.info("tpx file_tag      : {0}".format(th.file_tag))
        logger.info("tpx natoms        : {0}".format(th.natoms))
        logger.info("tpx ngtc          : {0}".format(th.ngtc))
        logger.info("tpx fep_state     : {0}".format(th.fep_state))
        logger.info("tpx lambda        : {0}".format(th.lamb))
        logger.debug("tpx bIr (input record): {0}".format(th.bIr))
        logger.debug("tpx bTop         : {0}".format(th.bTop))
        logger.debug("tpx bX           : {0}".format(th.bX))
        logger.debug("tpx bV           : {0}".format(th.bV))
        logger.debug("tpx bF           : {0}".format(th.bF))
        logger.debug("tpx bBox         : {0}".format(th.bBox))
