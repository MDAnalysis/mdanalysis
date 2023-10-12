# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
HOLE Analysis --- :mod:`MDAnalysis.analysis.hole2.helper`
=========================================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0

Helper functions used in :mod:`MDAnalysis.analysis.hole2.hole`
"""
import logging
import tempfile
import subprocess
import os
import numpy as np
import errno

from ...lib import util
from ...exceptions import ApplicationError
from .templates import (SIMPLE2_RAD, IGNORE_RESIDUES, hole_input,
                        hole_lines, exe_err)

logger = logging.getLogger(__name__)


def write_simplerad2(filename="simple2.rad"):
    """Write the built-in radii in :data:`SIMPLE2_RAD` to `filename`.

    Does nothing if `filename` already exists.

    Parameters
    ----------
    filename : str, optional
       output file name; the default is "simple2.rad"

    Returns
    -------
    filename : str
       returns the name of the data file
    """

    if not os.path.exists(filename):
        with open(filename, "w") as rad:
            rad.write(SIMPLE2_RAD + "\n")
        logger.debug("Created simple radii file {}".format(filename))
    return filename


def check_and_fix_long_filename(filename, tmpdir=os.path.curdir,
                                max_length=70,
                                make_symlink=True):
    """Return a file name suitable for HOLE.

    HOLE is limited to filenames <= ``max_length``. This method

    1. returns `filename` if HOLE can process it
    2. returns a relative path (see :func:`os.path.relpath`) if that shortens the
        path sufficiently
    3. creates a symlink to `filename` (:func:`os.symlink`) in a safe temporary
        directory and returns the path of the symlink.

    Parameters
    ----------
    filename : str
        file name to be processed
    tmpdir : str, optional
        By default the temporary directory is created inside the current
        directory in order to keep that path name short. This can be changed
        with the `tmpdir` keyword (e.g. one can use "/tmp"). The default is
        the current directory :data:`os.path.curdir`.

    Returns
    -------
    str
        path to the file that has a length less than
        ``max_length``

    Raises
    ------
    RuntimeError
        If none of the tricks for filename shortening worked. In this case,
        manually rename the file or recompile your version of HOLE.
    """

    if len(filename) <= max_length:
        return filename

    msg = ('HOLE will not read {} '
           'because it has more than {} characters.')
    logger.debug(msg.format(filename, max_length))

    # try a relative path
    new_name = os.path.relpath(filename)
    if len(new_name) <= max_length:
        msg = 'Using relative path: {} -> {}'
        logger.debug(msg.format(filename, new_name))
        return new_name

    if make_symlink:
        # shorten path by creating a symlink inside a safe temp dir
        _, ext = os.path.splitext(filename)
        dirname = os.path.relpath(tempfile.mkdtemp(dir=tmpdir))
        newname = os.path.join(dirname, os.path.basename(filename))
        if len(newname) > max_length:
            fd, newname = tempfile.mkstemp(suffix=ext, dir=dirname)
            os.close(fd)
            os.unlink(newname)

        if len(newname) > max_length:
            newname = os.path.relpath(newname)

        if len(newname) <= max_length:
            os.symlink(filename, newname)
            msg = 'Using symlink: {} -> {}'
            logger.debug(msg.format(filename, newname))
            return newname

    msg = 'Failed to shorten filename {}'
    raise RuntimeError(msg.format(filename))


def set_up_hole_input(pdbfile,
                      infile_text=None,
                      infile=None,
                      sphpdb_file='hole.sph',
                      vdwradii_file=None,
                      tmpdir=os.path.curdir,
                      sample=0.2,
                      end_radius=22,
                      cpoint=None,
                      cvect=None,
                      random_seed=None,
                      ignore_residues=IGNORE_RESIDUES,
                      output_level=0,
                      dcd=None,
                      dcd_iniskip=0,
                      dcd_step=1):
    """
    Generate HOLE input for the parameters.

    Parameters
    ----------
    pdbfile : str
        The `filename` is used as input for HOLE in the "COORD" card of the
        input file.  It specifies the name of a PDB coordinate file to be
        used. This must be in Brookhaven protein databank format or
        something closely approximating this. Both ATOM and HETATM records
        are read.

    infile_text: str, optional
        HOLE input text or template. If set to ``None``, the function will
        create the input text from the other parameters.
        Default: ``None``

    infile: str, optional
        File to write the HOLE input text for later inspection. If set to
        ``None``, the input text is not written out.
        Default: ``None``

    sphpdb_file : str, optional
        path to the HOLE sph file, a PDB-like file containing the
        coordinates of the pore centers.
        The coordinates are set to the sphere centres and the occupancies
        are the sphere radii. All centres are assigned the atom name QSS and
        residue name SPH and the residue number is set to the storage
        number of the centre. In VMD, sph
        objects are best displayed as "Points". Displaying .sph objects
        rather than rendered or dot surfaces can be useful to analyze the
        distance of particular atoms from the sphere-centre line.
        .sph files can be used to produce molecular graphical
        output from a hole run, by using the
        :program:`sph_process` program to read the .sph file.
        Default: 'hole.sph'

    vdwradii_file: str, optional
        path to the file specifying van der Waals radii for each atom. If
        set to ``None``, then a set of default radii,
        :data:`SIMPLE2_RAD`, is used (an extension of ``simple.rad`` from
        the HOLE distribution). Default: ``None``

    tmpdir: str, optional
        The temporary directory that files can be symlinked to, to shorten
        the path name. HOLE can only read filenames up to a certain length.
        Default: current working directory

    sample : float, optional
        distance of sample points in Å.
        Specifies the distance between the planes used in the HOLE
        procedure. The default value should be reasonable for most
        purposes. However, if you wish to visualize a very tight
        constriction then specify a smaller value.
        This value determines how many points in the pore profile are
        calculated. Default: 0.2

    end_radius : float, optional
        Radius in Å, which is considered to be the end of the pore. This
        keyword can be used to specify the radius above which the
        program regards a result as indicating that the end of the pore
        has been reached. This may need to be increased for large channels,
        or reduced for small channels. Default: 22.0

    cpoint : array_like, 'center_of_geometry' or None, optional
        coordinates of a point inside the pore, e.g. ``[12.3, 0.7,
        18.55]``. If set to ``None`` (the default) then HOLE's own search
        algorithm is used.
        ``cpoint`` specifies a point which lies within the channel. For
        simple channels (e.g. gramicidin), results do not show great
        sensitivity to the exact point taken. An easy way to produce an
        initial point is to use molecular graphics to find two atoms which
        lie either side of the pore and to average their coordinates. Or
        if the channel structure contains water molecules or counter ions
        then take the coordinates of one of these (and use the
        `ignore_residues` keyword to ignore them in the pore radius
        calculation).
        If this card is not specified, then HOLE (from version 2.2)
        attempts to guess where the channel will be. The procedure
        assumes the channel is reasonably symmetric. The initial guess on
        cpoint will be the centroid of all alpha carbon atoms (name 'CA'
        in pdb file). This is then refined by a crude grid search up to 5
        Å from the original position. This procedure works most of the
        time but is far from infallible — results should be
        carefully checked (with molecular graphics) if it is used.
        Default: None

    cvect : array_like, optional
        Search direction, should be parallel to the pore axis,
        e.g. ``[0,0,1]`` for the z-axis.
        If this keyword is ``None`` (the default), then HOLE attempts to guess
        where the channel will be. The procedure assumes that the channel is
        reasonably symmetric. The guess will be either along the X axis
        (1,0,0), Y axis (0,1,0) or Z axis (0,0,1). If the structure is not
        aligned on one of these axis the results will clearly be
        approximate. If a guess is used then results should be carefully
        checked. Default: None

    random_seed : int, optional
        integer number to start the random number generator.
        By default,
        :program:`hole` will use the time of the day.
        For reproducible runs (e.g., for testing) set ``random_seed``
        to an integer. Default: ``None``

    ignore_residues : array_like, optional
        sequence of three-letter residues that are not taken into
        account during the calculation; wildcards are *not*
        supported. Note that all residues must have 3 letters. Pad
        with space on the right-hand side if necessary.
        Default: {}.

    output_level : int, optional
        Determines the output of output in the ``outfile``.
        For automated processing, this must be < 3.
        0: Full text output,
        1: All text output given except "run in progress" (i.e.,
        detailed contemporary description of what HOLE is doing).
        2: Ditto plus no graph type output - only leaving minimum
        radius and conductance calculations.
        3: All text output other than input card mirroring and error messages
        turned off.
        Default: 0

    dcd : str, optional
        File name of DCD trajectory (must be supplied together with a
        matching PDB file `filename`) and then HOLE runs its analysis on
        each frame.
        It does multiple HOLE runs on positions taken from a CHARMM binary
        dynamics format DCD trajectory file. The ``dcd`` file must have
        exactly the same number of atoms in exactly the same order as the
        pdb file specified by ``pdbfile``. Note that if this option is used
        the pdb file is used as a template only - the coordinates are
        ignored. Note that structural parameters determined for each
        individual structure are written in a tagged format so that it is
        possible to extract the information from the text output file using
        a :program:`grep` command. The reading of the file can be
        controlled by the ``dcd_step`` keyword and/or setting
        ``dcd_iniskip`` to the number of frames to be skipped
        initially.

    dcd_step : int, optional
        step size for going through the trajectory (skips ``dcd_step-1``
        frames). Default: 1

    Returns
    -------
    str
        input text to run HOLE


    .. versionadded:: 1.0

    """.format(IGNORE_RESIDUES)

    short_filename = check_and_fix_long_filename(pdbfile, tmpdir=tmpdir)
    if vdwradii_file is not None:
        vdwradii_file = check_and_fix_long_filename(vdwradii_file,
                                                    tmpdir=tmpdir)
    else:
        vdwradii_file = write_simplerad2()

    if dcd is not None:
        dcd = check_and_fix_long_filename(dcd, tmpdir=tmpdir)

    if infile_text is None:
        infile_text = hole_input

    residues = ' '.join(ignore_residues)

    infile_text = infile_text.format(filename=pdbfile,
                                     coordinates=short_filename,
                                     radius=vdwradii_file,
                                     sphpdb=sphpdb_file,
                                     sample=sample,
                                     end_radius=end_radius,
                                     ignore=residues,
                                     output_level=output_level)

    if random_seed is not None:
        random_seed = int(random_seed)
        infile_text += hole_lines['random_seed'].format(random_seed)
        logger.info("Fixed random number seed {} for reproducible "
                    "runs.".format(random_seed))

    if cpoint is not None:
        if isinstance(cpoint, str):
            infile_text += 'CPOINT ' + cpoint + '\n'
        else:
            infile_text += hole_lines['cpoint'].format(*cpoint)
    else:
        logger.info("HOLE will guess CPOINT")

    if cvect is not None:
        infile_text += hole_lines['cvect'].format(*cvect)
    else:
        logger.info("HOLE will guess CVECT")

    if dcd is not None:
        infile_text += hole_lines['dcd'].format(dcd=dcd,
                                                iniskip=dcd_iniskip,
                                                step=dcd_step)

    if infile is not None:
        with open(infile, 'w') as f:
            f.write(infile_text)
        msg = 'Wrote HOLE input file {} for inspection'
        logger.debug(msg.format(infile))

    return infile_text


def run_hole(outfile, infile_text, executable):
    """Run the HOLE program.

    Parameters
    ----------
    outfile: str
        Output file name
    infile_text: str
        HOLE input text
        (typically generated by :func:`set_up_hole_input`)
    executable: str
        HOLE executable


    Returns
    -------
    str
        Output file name
    """
    with open(outfile, 'w') as output:
        proc = subprocess.Popen(executable, stdin=subprocess.PIPE,
                                stdout=output)
        stdout, stderr = proc.communicate(infile_text.encode('utf-8'))

    # check output in case of errors
    with open(outfile, 'r') as output:
        for line in output:
            if line.strip().startswith(('*** ERROR ***', 'ERROR')):
                proc.returncode = 255
                break

    # die in case of error
    if proc.returncode != 0:
        msg = 'HOLE failure ({}). Check output {}'
        logger.fatal(msg.format(proc.returncode, outfile))
        if stderr is not None:
            logger.fatal(stderr)
        raise ApplicationError(proc.returncode,
                               msg.format(executable, outfile))

    logger.info('HOLE finished. Output: {}'.format(outfile))
    return outfile


def collect_hole(outfile='hole.out'):
    """Collect data from HOLE output

    Parameters
    ----------
    outfile: str, optional
        HOLE output file to read. Default: 'hole.out'


    Returns
    -------
    dict
        Dictionary of HOLE profiles as record arrays
    """
    fmt = util.FORTRANReader('3F12')
    recarrays = {}

    with open(outfile, 'r') as output:
        toggle_read = False
        profile = 0
        records = []
        for line in output:
            line = line.rstrip()  # preserve columns in FORTRAN output
            # multiple frames from dcd in?
            if line.startswith(" Starting calculation for position number"):
                fields = line.split()
                profile = int(fields[5])
                records = []
                logger.debug('Started reading profile {}'.format(profile))
                continue

            # found data
            if line.startswith(' cenxyz.cvec'):
                toggle_read = True
                logger.debug('Started reading data')
                continue

            if toggle_read:
                if len(line.strip()) != 0:
                    try:
                        rxncoord, radius, cenlineD = fmt.read(line)
                    except:
                        msg = 'Problem parsing line: {}. Check output file {}'
                        logger.exception(msg.format(line, outfile))
                        raise
                    records.append((rxncoord, radius, cenlineD))
                    continue
                # end of data
                else:
                    toggle_read = False
                    names = ['rxn_coord', 'radius', 'cen_line_D']
                    recarr = np.rec.fromrecords(records,
                                                formats='f8, f8, f8',
                                                names=names)
                    recarrays[profile] = recarr

        return recarrays


def create_vmd_surface(sphpdb='hole.sph',
                       filename=None,
                       sph_process='sph_process',
                       sos_triangle='sos_triangle',
                       dot_density=15):
    """Create VMD surface file from sphpdb file.

    Parameters
    ----------
    sphpdb: str, optional
        sphpdb file to read. Default: 'hole.sph'
    filename: str, optional
        output VMD surface file. If ``None``, a temporary file
        is generated. Default: ``None``
    sph_process: str, optional
        Executable for ``sph_process`` program. Default: 'sph_process'
    sos_triangle: str, optional
        Executable for ``sos_triangle`` program. Default: 'sos_triangle'
    dot_density: int, optional
        density of facets for generating a 3D pore representation.
        The number controls the density of dots that will be used.
        A sphere of dots is placed on each centre determined in the
        Monte Carlo procedure. The actual number of dots written is
        controlled by ``dot_density`` and the ``sample`` level of the
        original analysis. ``dot_density`` should be set between 5
        (few dots per sphere) and 35 (many dots per sphere).
        Default: 15


    Returns
    -------
    str
        the output filename for the VMD surface

    """
    fd, tmp_sos = tempfile.mkstemp(suffix=".sos", text=True)
    os.close(fd)

    sph_process_path = util.which(sph_process)
    if sph_process_path is None:
        raise OSError(errno.ENOENT, exe_err.format(name=sph_process,
                                                   kw='sph_process'))
    base_path = os.path.dirname(sph_process_path)

    sos_triangle_path = util.which(sos_triangle)
    if sos_triangle_path is None:
        path = os.path.join(base_path, sos_triangle)
        sos_triangle_path = util.which(path)
    if sos_triangle_path is None:
        raise OSError(errno.ENOENT, exe_err.format(name=sos_triangle,
                                                   kw='sos_triangle'))
    try:
        output = subprocess.check_output([sph_process_path, "-sos", "-dotden",
                                          str(dot_density), "-color", sphpdb,
                                          tmp_sos], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as err:
        os.unlink(tmp_sos)
        logger.fatal("sph_process failed ({0})".format(err.returncode))
        raise OSError(err.returncode, "sph_process failed") from None
    except:
        os.unlink(tmp_sos)
        raise

    if filename is None:
        fd, filename = tempfile.mkstemp(suffix=".sos", text=True)
        os.close(fd)
    try:
        # Could check: os.devnull if subprocess.DEVNULL not available (>3.3)
        # Suppress stderr messages of sos_triangle
        with open(tmp_sos) as sos, open(filename, "w") as triangles, \
                open(os.devnull, 'w') as FNULL:
            subprocess.check_call(
                [sos_triangle_path, "-s"], stdin=sos, stdout=triangles,
                stderr=FNULL)
    except subprocess.CalledProcessError as err:
        logger.fatal("sos_triangle failed ({0})".format(err.returncode))
        raise OSError(err.returncode, "sos_triangle failed") from None
    finally:
        os.unlink(tmp_sos)

    return filename
