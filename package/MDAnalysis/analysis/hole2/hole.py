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
import os
import errno
import shutil
import tempfile
import textwrap
import logging
import itertools
import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import OrderedDict

from MDAnalysis.lib.util import deprecate
from ...exceptions import ApplicationError
from ..base import AnalysisBase
from ...lib import util
from .utils import (check_and_fix_long_filename, write_simplerad2,
                    set_up_hole_input, run_hole, collect_hole,
                    create_vmd_surface)
from .templates import (hole_input, hole_lines, vmd_script_array,
                        vmd_script_function, exe_err,
                        IGNORE_RESIDUES)

logger = logging.getLogger(__name__)


@deprecate(release="2.6.0", remove="3.0.0",
           message=("This method has been moved to the MDAKit hole2-mdakit: "
                    "https://github.com/MDAnalysis/hole2-mdakit"))
def hole(pdbfile,
         infile_text=None,
         infile=None,
         outfile='hole.out',
         sphpdb_file='hole.sph',
         vdwradii_file=None,
         executable='hole',
         tmpdir=os.path.curdir,
         sample=0.2,
         end_radius=22.0,
         cpoint=None,
         cvect=None,
         random_seed=None,
         ignore_residues=IGNORE_RESIDUES,
         output_level=0,
         dcd=None,
         dcd_iniskip=0,
         dcd_step=1,
         keep_files=True):
    r"""Run :program:`hole` on a single frame or a DCD trajectory.

    :program:`hole` is part of the HOLE_ suite of programs. It is used to
    analyze channels and cavities in proteins, especially ion channels.

    Only a subset of all `HOLE control parameters <http://www.holeprogram.org/doc/old/hole_d03.html>`_
    is supported and can be set with keyword arguments.

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
    infile: str, optional
        File to write the HOLE input text for later inspection. If set to
        ``None``, the input text is not written out.
    outfile : str, optional
        file name of the file collecting HOLE's output (which can be
        parsed using :meth:`collect_hole(outfile)`.
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
    vdwradii_file: str, optional
        path to the file specifying van der Waals radii for each atom. If
        set to ``None``, then a set of default radii,
        :data:`SIMPLE2_RAD`, is used (an extension of ``simple.rad`` from
        the HOLE distribution).
    executable: str, optional
        Path to the :program:`hole` executable.
        (e.g. ``~/hole2/exe/hole``). If
        :program:`hole` is found on the :envvar:`PATH`, then the bare
        executable name is sufficient.
    tmpdir: str, optional
        The temporary directory that files can be symlinked to, to shorten
        the path name. HOLE can only read filenames up to a certain length.
    sample : float, optional
        distance of sample points in Å.
        Specifies the distance between the planes used in the HOLE
        procedure. The default value should be reasonable for most
        purposes. However, if you wish to visualize a very tight
        constriction then specify a smaller value.
        This value determines how many points in the pore profile are
        calculated.
    end_radius : float, optional
        Radius in Å, which is considered to be the end of the pore. This
        keyword can be used to specify the radius above which the
        program regards a result as indicating that the end of the pore
        has been reached. This may need to be increased for large channels,
        or reduced for small channels.
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
        ``ignore_residues`` keyword to ignore them in the pore radius
        calculation).
        If this card is not specified, then HOLE (from version 2.2)
        attempts to guess where the channel will be. The procedure
        assumes the channel is reasonably symmetric. The initial guess on
        cpoint will be the centroid of all alpha carbon atoms (name 'CA'
        in pdb file). This is then refined by a crude grid search up to 5
        Å from the original position. This procedure works most of the
        time but is far from infallible — results should be
        carefully checked (with molecular graphics) if it is used.
    cvect : array_like, optional
        Search direction, should be parallel to the pore axis,
        e.g. ``[0,0,1]`` for the z-axis.
        If this keyword is ``None`` (the default), then HOLE attempts to guess
        where the channel will be. The procedure assumes that the channel is
        reasonably symmetric. The guess will be either along the X axis
        (1,0,0), Y axis (0,1,0) or Z axis (0,0,1). If the structure is not
        aligned on one of these axis the results will clearly be
        approximate. If a guess is used then results should be carefully
        checked.
    random_seed : int, optional
        integer number to start the random number generator.
        By default,
        :program:`hole` will use the time of the day.
        For reproducible runs (e.g., for testing) set ``random_seed``
        to an integer.
    ignore_residues : array_like, optional
        sequence of three-letter residues that are not taken into
        account during the calculation; wildcards are *not*
        supported. Note that all residues must have 3 letters. Pad
        with space on the right-hand side if necessary.
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
    dcd : str, optional
        File name of CHARMM-style DCD trajectory (must be supplied together with a
        matching PDB file `filename`) and then HOLE runs its analysis on
        each frame. HOLE can *not* read DCD trajectories written by MDAnalysis,
        which are NAMD-style (see Notes). Note that structural parameters
        determined for each individual structure are written in a tagged
        format so that it is possible to extract the information from the text
        output file using a :program:`grep` command. The reading of the file
        can be controlled by the ``dcd_step`` keyword and/or setting
        ``dcd_iniskip`` to the number of frames to be skipped
        initially.
    dcd_step : int, optional
        step size for going through the trajectory (skips ``dcd_step-1``
        frames).
    keep_files : bool, optional
        Whether to keep the HOLE output files and possible temporary
        symlinks after running the function.


    Returns
    -------

    dict
        A dictionary of :class:`numpy.recarray`\ s, indexed by frame.


    Notes
    -----
    - HOLE is very picky and does not read all DCD-like formats [#HOLEDCD]_.
      If in doubt, look into the `outfile` for error diagnostics.


    .. versionadded:: 1.0

    """

    if output_level > 3:
        msg = 'output_level ({}) needs to be < 3 in order to extract a HOLE profile!'
        warnings.warn(msg.format(output_level))

    # get executable
    exe = shutil.which(executable)
    if exe is None:
        raise OSError(errno.ENOENT, exe_err.format(name=executable,
                                                   kw='executable'))

    # get temp files
    tmp_files = [outfile, sphpdb_file]

    short_filename = check_and_fix_long_filename(pdbfile, tmpdir=tmpdir)
    if os.path.islink(short_filename):
        tmp_files.append(short_filename)

    if dcd is not None:
        dcd = check_and_fix_long_filename(dcd, tmpdir=tmpdir)
        if os.path.islink(dcd):
            tmp_files.append(dcd)

    if vdwradii_file is not None:
        vdwradii_file = check_and_fix_long_filename(vdwradii_file,
                                                    tmpdir=tmpdir)
    else:
        vdwradii_file = write_simplerad2()
        tmp_files.append(vdwradii_file)

    infile_text = set_up_hole_input(short_filename,
                                    infile_text=infile_text,
                                    infile=infile,
                                    sphpdb_file=sphpdb_file,
                                    vdwradii_file=vdwradii_file,
                                    tmpdir=tmpdir, sample=sample,
                                    end_radius=end_radius,
                                    cpoint=cpoint, cvect=cvect,
                                    random_seed=random_seed,
                                    ignore_residues=ignore_residues,
                                    output_level=output_level,
                                    dcd=dcd,
                                    dcd_iniskip=dcd_iniskip,
                                    dcd_step=dcd_step-1)

    run_hole(outfile=outfile, infile_text=infile_text, executable=exe)
    recarrays = collect_hole(outfile=outfile)

    if not keep_files:
        for file in tmp_files:
            try:
                os.unlink(file)
            except OSError:
                pass

    return recarrays


class HoleAnalysis(AnalysisBase):
    r"""
    Run :program:`hole` on a trajectory.

    :program:`hole` is part of the HOLE_ suite of programs. It is used to
    analyze channels and cavities in proteins, especially ion channels.

    Only a subset of all `HOLE control parameters <http://www.holeprogram.org/doc/old/hole_d03.html>`_
    is supported and can be set with keyword arguments.

    This class creates temporary PDB files for each frame and runs HOLE on
    the frame. It can be used normally, or as a context manager. If used as a
    context manager, the class will try to delete any temporary files created
    by HOLE, e.g. sphpdb files and logfiles. ::

        with hole2.HoleAnalysis(u, executable='~/hole2/exe/hole') as h2:
            h2.run()
            h2.create_vmd_surface()

    Parameters
    ----------

    universe : Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to.
    select : string, optional
        The selection string to create an atom selection that the HOLE
        analysis is applied to.
    vdwradii_file : str, optional
        path to the file specifying van der Waals radii for each atom. If
        set to ``None``, then a set of default radii,
        :data:`SIMPLE2_RAD`, is used (an extension of ``simple.rad`` from
        the HOLE distribution).
    executable : str, optional
        Path to the :program:`hole` executable.
        (e.g. ``~/hole2/exe/hole``). If
        :program:`hole` is found on the :envvar:`PATH`, then the bare
        executable name is sufficient.
    tmpdir : str, optional
        The temporary directory that files can be symlinked to, to shorten
        the path name. HOLE can only read filenames up to a certain length.
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
        ``ignore_residues`` keyword to ignore them in the pore radius
        calculation).
        If this card is not specified, then HOLE (from version 2.2)
        attempts to guess where the channel will be. The procedure
        assumes the channel is reasonably symmetric. The initial guess on
        cpoint will be the centroid of all alpha carbon atoms (name 'CA'
        in pdb file). This is then refined by a crude grid search up to 5
        Å from the original position. This procedure works most of the
        time but is far from infallible — results should be
        carefully checked (with molecular graphics) if it is used.
    cvect : array_like, optional
        Search direction, should be parallel to the pore axis,
        e.g. ``[0,0,1]`` for the z-axis.
        If this keyword is ``None`` (the default), then HOLE attempts to guess
        where the channel will be. The procedure assumes that the channel is
        reasonably symmetric. The guess will be either along the X axis
        (1,0,0), Y axis (0,1,0) or Z axis (0,0,1). If the structure is not
        aligned on one of these axis the results will clearly be
        approximate. If a guess is used then results should be carefully
        checked.
    sample : float, optional
        distance of sample points in Å.
        Specifies the distance between the planes used in the HOLE
        procedure. The default value should be reasonable for most
        purposes. However, if you wish to visualize a very tight
        constriction then specify a smaller value.
        This value determines how many points in the pore profile are
        calculated.
    end_radius : float, optional
        Radius in Å, which is considered to be the end of the pore. This
        keyword can be used to specify the radius above which the
        program regards a result as indicating that the end of the pore
        has been reached. This may need to be increased for large channels,
        or reduced for small channels.
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
    ignore_residues : array_like, optional
        sequence of three-letter residues that are not taken into
        account during the calculation; wildcards are *not*
        supported. Note that all residues must have 3 letters. Pad
        with space on the right-hand side if necessary.
    prefix : str, optional
        Prefix for HOLE output files.
    write_input_files : bool, optional
        Whether to write out the input HOLE text as files.
        Files are called `hole.inp`.


    Attributes
    ----------
    results.sphpdbs: numpy.ndarray
        Array of sphpdb filenames

        .. versionadded:: 2.0.0

    results.outfiles: numpy.ndarray
        Arrau of output filenames

        .. versionadded:: 2.0.0

    results.profiles: dict
        Profiles generated by HOLE2.
        A dictionary of :class:`numpy.recarray`\ s, indexed by frame.

        .. versionadded:: 2.0.0

    sphpdbs: numpy.ndarray
        Alias of :attr:`results.sphpdbs`

        .. deprecated:: 2.0.0
            This will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.sphpdbs` instead.

    outfiles: numpy.ndarray
        Alias of :attr:`results.outfiles`

        .. deprecated:: 2.0.0
            This will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.outfiles` instead.

    profiles: dict
        Alias of :attr:`results.profiles`

        .. deprecated:: 2.0.0
            This will be removed in MDAnalysis 3.0.0. Please use
            :attr:`results.profiles` instead.

    .. versionadded:: 1.0

    .. versionchanged:: 2.0.0
        :attr:`sphpdbs`, :attr:`outfiles` and :attr:`profiles `
        are now stored in a :class:`MDAnalysis.analysis.base.Results`
        instance.

    .. deprecated:: 2.6.0
       This class has been moved to the MDAKit
       `hole2-mdakit <https://github.com/MDAnalysis/hole2-mdakit>`_ and will
       be removed for the core MDAnalysis library in version 3.0

    """

    input_file = '{prefix}hole{i:03d}.inp'
    output_file = '{prefix}hole{i:03d}.out'
    sphpdb_file = '{prefix}hole{i:03d}.sph'

    input_file = '{prefix}hole{i:03d}.inp'
    output_file = '{prefix}hole{i:03d}.out'
    sphpdb_file = '{prefix}hole{i:03d}.sph'

    hole_header = textwrap.dedent("""
        ! Input file for Oliver Smart's HOLE program
        ! written by MDAnalysis.analysis.hole2.HoleAnalysis
        ! for a Universe
        ! u = mda.Universe({}
        !   )
        ! Frame {{i}}
        """)

    hole_body = textwrap.dedent("""
        COORD  {{coordinates}}
        RADIUS {radius}
        SPHPDB {{sphpdb}}
        SAMPLE {sample:f}
        ENDRAD {end_radius:f}
        IGNORE {ignore}
        SHORTO {output_level:d}
        """)

    _guess_cpoint = False

    def __init__(self, universe,
                 select='protein',
                 verbose=False,
                 ignore_residues=IGNORE_RESIDUES,
                 vdwradii_file=None,
                 executable='hole',
                 sos_triangle='sos_triangle',
                 sph_process='sph_process',
                 tmpdir=os.path.curdir,
                 cpoint=None,
                 cvect=None,
                 sample=0.2,
                 end_radius=22,
                 output_level=0,
                 prefix=None,
                 write_input_files=False):
        super(HoleAnalysis, self).__init__(universe.universe.trajectory,
                                           verbose=verbose)

        wmsg = ("This class has been moved to the MDAKit `hole2-mdakit` "
                "(https://github.com/MDAnalysis/hole2-mdakit) and will be "
                "removed in version 3.0.")
        warnings.warn(wmsg, DeprecationWarning)

        if output_level > 3:
            msg = 'output_level ({}) needs to be < 3 in order to extract a HOLE profile!'
            warnings.warn(msg.format(output_level))

        if prefix is None:
            prefix = ''

        if isinstance(cpoint, str):
            if 'geometry' in cpoint.lower():
                self._guess_cpoint = True
                self.cpoint = '{cpoint[0]:.10f} {cpoint[1]:.10f} {cpoint[2]:.10f}'
        else:
            self._guess_cpoint = False
            self.cpoint = cpoint

        self.prefix = prefix
        self.cvect = cvect
        self.sample = sample
        self.end_radius = end_radius
        self.output_level = output_level
        self.write_input_files = write_input_files
        self.select = select
        self.ag = universe.select_atoms(select, updating=True)
        self.universe = universe
        self.tmpdir = tmpdir
        self.ignore_residues = ignore_residues

        # --- finding executables ----
        hole = shutil.which(executable)
        if hole is None:
            raise OSError(errno.ENOENT, exe_err.format(name=executable,
                                                       kw='executable'))
        self.base_path = os.path.dirname(hole)

        sos_triangle_path = shutil.which(sos_triangle)
        if sos_triangle_path is None:
            path = os.path.join(self.base_path, sos_triangle)
            sos_triangle_path = shutil.which(path)
        if sos_triangle_path is None:
            raise OSError(errno.ENOENT, exe_err.format(name=sos_triangle,
                                                       kw='sos_triangle'))
        sph_process_path = shutil.which(sph_process)
        if sph_process_path is None:
            path = os.path.join(self.base_path, sph_process)
            sph_process_path = shutil.which(path)
        if sph_process_path is None:
            raise OSError(errno.ENOENT, exe_err.format(name=sph_process,
                                                       kw='sph_process'))

        self.exe = {
            'hole': hole,
            'sos_triangle': sos_triangle_path,
            'sph_process': sph_process_path
        }

        # --- setting up temp files ----
        self.tmp_files = []
        if vdwradii_file is not None:
            self.vdwradii_file = check_and_fix_long_filename(vdwradii_file,
                                                             tmpdir=self.tmpdir)
            if os.path.islink(self.vdwradii_file):
                self.tmp_files.append(self.vdwradii_file)
        else:
            self.vdwradii_file = write_simplerad2()
            self.tmp_files.append(self.vdwradii_file)

        # --- setting up input header ----
        filenames = [universe.filename]
        try:
            filenames.extend(universe.trajectory.filenames)
        except AttributeError:
            filenames.append(universe.trajectory.filename)
        filenames = [name for name in filenames if name is not None]
        hole_filenames = '\n!    '.join(filenames)
        self._input_header = self.hole_header.format(hole_filenames)

    def run(self, start=None, stop=None, step=None, verbose=None,
            random_seed=None):
        """
        Perform the calculation

        Parameters
        ----------
        start : int, optional
            start frame of analysis

        stop : int, optional
            stop frame of analysis

        step : int, optional
            number of frames to skip between each analysed frame

        verbose : bool, optional
            Turn on verbosity

        random_seed : int, optional
            integer number to start the random number generator.
            By default,
            :program:`hole` will use the time of the day.
            For reproducible runs (e.g., for testing) set ``random_seed``
            to an integer.
        """
        self.random_seed = random_seed
        return super(HoleAnalysis, self).run(start=start, stop=stop,
                                             step=step, verbose=verbose)

    @property
    def sphpdbs(self):
        wmsg = ("The `sphpdbs` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.sphpdbs` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.sphpdbs

    @property
    def outfiles(self):
        wmsg = ("The `outfiles` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.outfiles` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.outfiles

    @property
    def profiles(self):
        wmsg = ("The `profiles` attribute was deprecated in "
                "MDAnalysis 2.0.0 and will be removed in MDAnalysis 3.0.0. "
                "Please use `results.profiles` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.profiles

    def _prepare(self):
        """Set up containers and generate input file text"""
        # set up containers
        self.results.sphpdbs = np.zeros(self.n_frames, dtype=object)
        self.results.outfiles = np.zeros(self.n_frames, dtype=object)
        self.results.profiles = {}

        # generate input file
        body = set_up_hole_input('',
                                 infile_text=self.hole_body,
                                 infile=None,
                                 vdwradii_file=self.vdwradii_file,
                                 tmpdir=self.tmpdir,
                                 sample=self.sample,
                                 end_radius=self.end_radius,
                                 cpoint=self.cpoint,
                                 cvect=self.cvect,
                                 random_seed=self.random_seed,
                                 ignore_residues=self.ignore_residues,
                                 output_level=self.output_level,
                                 dcd=None)

        self.infile_text = self._input_header + body

    def guess_cpoint(self):
        """Guess a point inside the pore.

        This method simply uses the center of geometry of the selection as a
        guess.

        Returns
        -------
        float:
            center of geometry of selected AtomGroup
        """
        return self.ag.center_of_geometry()

    def _single_frame(self):
        """Run HOLE analysis and collect profiles"""
        # set up files
        frame = self._ts.frame
        i = self._frame_index
        outfile = self.output_file.format(prefix=self.prefix, i=frame)
        sphpdb = self.sphpdb_file.format(prefix=self.prefix, i=frame)
        self.results.sphpdbs[i] = sphpdb
        self.results.outfiles[i] = outfile
        if outfile not in self.tmp_files:
            self.tmp_files.append(outfile)
        if sphpdb not in self.tmp_files:
            self.tmp_files.append(sphpdb)
        else:
            self.tmp_files.append(sphpdb + '.old')

        # temp pdb
        logger.info('HOLE analysis frame {}'.format(frame))
        fd, pdbfile = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)  # close immediately (Issue 129)

        # get infile text
        fmt_kwargs = {'i': frame, 'coordinates': pdbfile, 'sphpdb': sphpdb}
        if self._guess_cpoint:
            fmt_kwargs['cpoint'] = self.guess_cpoint()
        infile_text = self.infile_text.format(**fmt_kwargs)

        if self.write_input_files:
            infile = self.input_file.format(prefix=self.prefix, i=frame)
            with open(infile, 'w') as f:
                f.write(infile_text)

        try:
            self.ag.write(pdbfile)
            run_hole(outfile=outfile, infile_text=infile_text,
                     executable=self.exe['hole'])
        finally:
            try:
                os.unlink(pdbfile)
            except OSError:
                pass

        recarrays = collect_hole(outfile=outfile)
        try:
            self.results.profiles[frame] = recarrays[0]
        except KeyError:
            msg = 'No profile found in HOLE output. Output level: {}'
            logger.info(msg.format(self.output_level))

    def create_vmd_surface(self, filename='hole.vmd', dot_density=15,
                           no_water_color='red', one_water_color='green',
                           double_water_color='blue'):
        """Process HOLE output to create a smooth pore surface suitable for VMD.

        Takes the ``sphpdb`` file for each frame and feeds it to `sph_process
        <http://www.holeprogram.org/doc/old/hole_d04.html#sph_process>`_ and
        `sos_triangle
        <http://www.holeprogram.org/doc/old/hole_d04.html#sos_triangle>`_ as
        described under `Visualization of HOLE results
        <http://www.holeprogram.org/doc/index.html>`_.

        Load the output file *filename* into VMD in :menuselection:`Extensions
        --> Tk Console` ::

           source hole.vmd

        The level of detail is determined by ``dot_density``.
        The surface will be colored by ``no_water_color``, ``one_water_color``, and
        ``double_water_color``. You can change these in the
        Tk Console::

            set no_water_color blue


        Parameters
        ----------
        filename: str, optional
            file to write the pore surfaces to.

        dot_density: int, optional
            density of facets for generating a 3D pore representation.
            The number controls the density of dots that will be used.
            A sphere of dots is placed on each centre determined in the
            Monte Carlo procedure. The actual number of dots written is
            controlled by ``dot_density`` and the ``sample`` level of the
            original analysis. ``dot_density`` should be set between 5
            (few dots per sphere) and 35 (many dots per sphere).

        no_water_color: str, optional
            Color of the surface where the pore radius is too tight for a
            water molecule.

        one_water_color: str, optional
            Color of the surface where the pore can fit one water molecule.

        double_water_color: str, optional
            Color of the surface where the radius is at least double the
            minimum radius for one water molecule.


        Returns
        -------
        str
            ``filename`` with the pore surfaces.

        """
        if not np.any(self.results.get("sphpdbs", [])):
            raise ValueError('No sphpdb files to read. Try calling run()')

        frames = []
        for i, frame in enumerate(self.frames):
            sphpdb = self.results.sphpdbs[i]
            tmp_tri = create_vmd_surface(sphpdb=sphpdb,
                                         sph_process=self.exe['sph_process'],
                                         sos_triangle=self.exe['sos_triangle'],
                                         dot_density=dot_density)

            shapes = [[], [], []]
            with open(tmp_tri) as f:
                for line in f:
                    if line.startswith('draw color'):
                        color = line.split()[-1].lower()
                        if color == 'red':
                            dest = shapes[0]
                        elif color == 'green':
                            dest = shapes[1]
                        elif color == 'blue':
                            dest = shapes[2]
                        else:
                            msg = 'Encountered unknown color {}'
                            raise ValueError(msg.format(color))

                    if line.startswith('draw trinorm'):
                        line = line.strip('draw trinorm').strip()
                        dest.append('{{ {} }}'.format(line))
            try:
                os.unlink(tmp_tri)
            except OSError:
                pass

            tri = '{ { ' + ' } { '.join(list(map(' '.join, shapes))) + ' } }'
            frames.append(f'set triangles({i}) ' + tri)

        trinorms = '\n'.join(frames)
        vmd_1 = vmd_script_array.format(no_water_color=no_water_color,
                                        one_water_color=one_water_color,
                                        double_water_color=double_water_color)
        vmd_text = vmd_1 + trinorms + vmd_script_function

        with open(filename, 'w') as f:
            f.write(vmd_text)

        return filename

    def min_radius(self):
        """Return the minimum radius over all profiles as a function of q"""
        profiles = self.results.get("profiles")
        if not profiles:
            raise ValueError('No profiles available. Try calling run()')
        return np.array([[q, p.radius.min()] for q, p in profiles.items()])

    def delete_temporary_files(self):
        """Delete temporary files"""
        for f in self.tmp_files:
            try:
                os.unlink(f)
            except OSError:
                pass
        self.tmp_files = []
        self.results.outfiles = []
        self.results.sphpdbs = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Delete temporary files on exit"""
        self.delete_temporary_files()

    def _process_plot_kwargs(self, frames=None,
                             color=None, cmap='viridis',
                             linestyle='-'):
        """Process the colors and linestyles for plotting

        Parameters
        ----------
        frames : array-like, optional
            Frames to plot. If ``None``, plots all of them.
        color : str or array-like, optional
            Color or colors for the plot. If ``None``, colors are
            drawn from ``cmap``.
        cmap : str, optional
            color map to make colors for the plot if ``color`` is
            not given. Names should be from the ``matplotlib.pyplot.cm``
            module.
        linestyle : str or array-like, optional
            Line style for the plot.


        Returns
        -------
        (array-like, array-like, array-like)
            frames, colors, linestyles
        """

        if frames is None:
            frames = self.frames
        else:
            frames = util.asiterable(frames)

        if color is None:
            colormap = matplotlib.colormaps.get_cmap(cmap)
            norm = matplotlib.colors.Normalize(vmin=min(frames),
                                               vmax=max(frames))
            colors = colormap(norm(frames))
        else:
            colors = itertools.cycle(util.asiterable(color))

        linestyles = itertools.cycle(util.asiterable(linestyle))

        return frames, colors, linestyles

    def plot(self, frames=None,
             color=None, cmap='viridis',
             linestyle='-', y_shift=0.0,
             label=True, ax=None,
             legend_loc='best', **kwargs):
        r"""Plot HOLE profiles :math:`R(\zeta)` in a 1D graph.

        Lines are colored according to the specified ``color`` or
        drawn from the color map ``cmap``. One line is
        plotted for each trajectory frame.

        Parameters
        ----------
        frames: array-like, optional
            Frames to plot. If ``None``, plots all of them.
        color: str or array-like, optional
            Color or colors for the plot. If ``None``, colors are
            drawn from ``cmap``.
        cmap: str, optional
            color map to make colors for the plot if ``color`` is
            not given. Names should be from the ``matplotlib.pyplot.cm``
            module.
        linestyle: str or array-like, optional
            Line style for the plot.
        y_shift : float, optional
            displace each :math:`R(\zeta)` profile by ``y_shift`` in the
            :math:`y`-direction for clearer visualization.
        label : bool or string, optional
            If ``False`` then no legend is
            displayed.
        ax : :class:`matplotlib.axes.Axes`
            If no `ax` is supplied or set to ``None`` then the plot will
            be added to the current active axes.
        legend_loc : str, optional
            Location of the legend.
        kwargs :  `**kwargs`
            All other `kwargs` are passed to :func:`matplotlib.pyplot.plot`.


        Returns
        -------
        ax : :class:`~matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """

        if not self.results.get("profiles"):
            raise ValueError('No profiles available. Try calling run()')

        if ax is None:
            fig, ax = plt.subplots()

        fcl = self._process_plot_kwargs(frames=frames, color=color,
                                        cmap=cmap, linestyle=linestyle)

        for i, (frame, c, ls) in enumerate(zip(*fcl)):
            profile = self.results.profiles[frame]
            dy = i*y_shift
            ax.plot(profile.rxn_coord, profile.radius+dy, color=c,
                    linestyle=ls, zorder=-frame, label=str(frame),
                    **kwargs)

        ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
        ax.set_ylabel(r"HOLE radius $R$ ($\AA$)")
        if label == True:
            ax.legend(loc=legend_loc)
        return ax

    def plot3D(self, frames=None,
               color=None, cmap='viridis',
               linestyle='-', ax=None, r_max=None,
               ylabel='Frames', **kwargs):
        r"""Stacked 3D graph of profiles :math:`R(\zeta)`.

        Lines are colored according to the specified ``color`` or
        drawn from the color map ``cmap``. One line is
        plotted for each trajectory frame.

        Parameters
        ----------
        frames : array-like, optional
            Frames to plot. If ``None``, plots all of them.
        color : str or array-like, optional
            Color or colors for the plot. If ``None``, colors are
            drawn from ``cmap``.
        cmap : str, optional
            color map to make colors for the plot if ``color`` is
            not given. Names should be from the ``matplotlib.pyplot.cm``
            module.
        linestyle : str or array-like, optional
            Line style for the plot.
        r_max : float, optional
            only display radii up to ``r_max``. If ``None``, all radii are
            plotted.
        ax : :class:`matplotlib.axes.Axes`
            If no `ax` is supplied or set to ``None`` then the plot will
            be added to the current active axes.
        ylabel : str, optional
            Y-axis label.
        **kwargs :
            All other `kwargs` are passed to :func:`matplotlib.pyplot.plot`.


        Returns
        -------
        ax : :class:`~mpl_toolkits.mplot3d.Axes3D`
             Axes with the plot, either `ax` or the current axes.

        """

        if not self.results.get("profiles"):
            raise ValueError('No profiles available. Try calling run()')

        from mpl_toolkits.mplot3d import Axes3D

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        fcl = self._process_plot_kwargs(frames=frames,
                                        color=color, cmap=cmap,
                                        linestyle=linestyle)

        for frame, c, ls in zip(*fcl):
            profile = self.results.profiles[frame]
            if r_max is None:
                radius = profile.radius
                rxn_coord = profile.rxn_coord
            else:
                # does not seem to work with masked arrays but with nan hack!
                # http://stackoverflow.com/questions/4913306/python-matplotlib-mplot3d-how-do-i-set-a-maximum-value-for-the-z-axis
                rxn_coord = profile.rxn_coord
                radius = profile.radius.copy()
                radius[radius > r_max] = np.nan
            ax.plot(rxn_coord, frame*np.ones_like(rxn_coord), radius,
                    color=c, linestyle=ls, zorder=-frame, label=str(frame),
                    **kwargs)

        ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
        ax.set_ylabel(ylabel)
        ax.set_zlabel(r"HOLE radius $R$ ($\AA$)")
        plt.tight_layout()

        return ax

    def over_order_parameters(self, order_parameters, frames=None):
        """Get HOLE profiles sorted over order parameters ``order_parameters``.

        Parameters
        ----------
        order_parameters : array-like or string
            Sequence or text file containing order parameters (float
            numbers) corresponding to the frames in the trajectory. Must
            be same length as trajectory.
        frames : array-like, optional
            Selected frames to return. If ``None``, returns all of them.

        Returns
        -------
        collections.OrderedDict
            sorted dictionary of {order_parameter:profile}

        """
        if not self.results.get("profiles"):
            raise ValueError('No profiles available. Try calling run()')
        if isinstance(order_parameters, str):
            try:
                order_parameters = np.loadtxt(order_parameters)
            except IOError:
                raise ValueError('Data file not found: {}'.format(order_parameters))
            except (ValueError, TypeError):
                msg = ('Could not parse given file: {}. '
                       '`order_parameters` must be array-like '
                       'or a filename with array data '
                       'that can be read by np.loadtxt')
                raise ValueError(msg.format(order_parameters))


        order_parameters = np.asarray(order_parameters)

        if len(order_parameters) != len(self._trajectory):
            msg = ('The number of order parameters ({}) must match the '
                   'length of the trajectory ({} frames)')
            raise ValueError(msg.format(len(order_parameters),
                                        len(self._trajectory)))

        if frames is None:
            frames = self.frames
        else:
            frames = np.asarray(util.asiterable(frames))

        idx = np.argsort(order_parameters[frames])
        sorted_frames = frames[idx]

        profiles = OrderedDict()
        for frame in sorted_frames:
            profiles[order_parameters[frame]] = self.results.profiles[frame]

        return profiles

    def plot_order_parameters(self, order_parameters,
                              aggregator=min,
                              frames=None,
                              color='blue',
                              linestyle='-', ax=None,
                              ylabel=r'Minimum HOLE pore radius $r$ ($\AA$)',
                              xlabel='Order parameter',
                              **kwargs):
        r"""Plot HOLE radii over order parameters. This function needs
        an ``aggregator`` function to reduce the ``radius`` array to a
        single value, e.g. ``min``, ``max``, or ``np.mean``.

        Parameters
        ----------
        order_parameters : array-like or string
            Sequence or text file containing order parameters (float
            numbers) corresponding to the frames in the trajectory. Must
            be same length as trajectory.
        aggregator : callable, optional
            Function applied to the radius array of each profile to
            reduce it to one representative value.
        frames : array-like, optional
            Frames to plot. If ``None``, plots all of them.
        color : str or array-like, optional
            Color for the plot.
        linestyle : str or array-like, optional
            Line style for the plot.
        ax : :class:`matplotlib.axes.Axes`
            If no `ax` is supplied or set to ``None`` then the plot will
            be added to the current active axes.
        xlabel : str, optional
            X-axis label.
        ylabel : str, optional
            Y-axis label.
        **kwargs :
            All other `kwargs` are passed to :func:`matplotlib.pyplot.plot`.


        Returns
        -------
        ax : :class:`~matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """

        op_profiles = self.over_order_parameters(order_parameters,
                                                 frames=frames)

        if ax is None:
            fig, ax = plt.subplots()

        data = [[x, aggregator(p.radius)] for x, p in op_profiles.items()]
        arr = np.array(data)
        ax.plot(arr[:, 0], arr[:, 1], color=color, linestyle=linestyle)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return ax

    def gather(self, frames=None, flat=False):
        """Gather the fields of each profile recarray together.

        Parameters
        ----------
        frames : int or iterable of ints, optional
            Profiles to include by frame. If ``None``, includes
            all frames.
        flat : bool, optional
            Whether to flatten the list of field arrays into a
            single array.

        Returns
        -------
        dict
            dictionary of fields
        """
        if frames is None:
            frames = self.frames
        frames = util.asiterable(frames)
        profiles = [self.results.profiles[k] for k in frames]

        rxncoords = [p.rxn_coord for p in profiles]
        radii = [p.radius for p in profiles]
        cen_line_Ds = [p.cen_line_D for p in profiles]

        if flat:
            rxncoords = np.concatenate(rxncoords)
            radii = np.concatenate(radii)
            cen_line_Ds = np.concatenate(cen_line_Ds)

        dct = {'rxn_coord': rxncoords,
               'radius': radii,
               'cen_line_D': cen_line_Ds}
        return dct

    def bin_radii(self, frames=None, bins=100, range=None):
        """Collects the pore radii into bins by reaction coordinate.

        Parameters
        ----------
        frames : int or iterable of ints, optional
            Profiles to include by frame. If ``None``, includes
            all frames.
        bins : int or iterable of edges, optional
            If bins is an int, it defines the number of equal-width bins in the given
            range. If bins is a sequence, it defines a monotonically increasing array of
            bin edges, including the rightmost edge, allowing for non-uniform bin widths.
        range : (float, float), optional
            The lower and upper range of the bins.
            If not provided, ``range`` is simply ``(a.min(), a.max())``,
            where ``a`` is the array of reaction coordinates.
            Values outside the range are ignored. The first element of the range must be
            less than or equal to the second.


        Returns
        -------
        list of arrays of floats
            List of radii present in each bin
        array of (float, float)
            Edges of each bin
        """
        agg = self.gather(frames=frames, flat=True)
        coords = agg['rxn_coord']

        if not util.iterable(bins):
            if range is None:
                range = (coords.min(), coords.max())
            xmin, xmax = range
            if xmin == xmax:
                xmin -= 0.5
                xmax += 0.5
            bins = np.linspace(xmin, xmax, bins+1, endpoint=True)
        else:
            bins = np.asarray(bins)
            bins = bins[np.argsort(bins)]

        idx = np.argsort(coords)
        coords = coords[idx]
        radii = agg['radius'][idx]
        # left: inserts at i where coords[:i] < edge
        # right: inserts at i where coords[:i] <= edge
        # r_ concatenates
        bin_idx = np.r_[coords.searchsorted(bins, side='right')]
        binned = [radii[i:j] for i, j in zip(bin_idx[:-1], bin_idx[1:])]
        return binned, bins

    def histogram_radii(self, aggregator=np.mean, frames=None,
                        bins=100, range=None):
        """Histograms the pore radii into bins by reaction coordinate,
        aggregate the radii with an `aggregator` function, and returns the
        aggregated radii and bin edges.

        Parameters
        ----------
        aggregator : callable, optional
            this function must take an iterable of floats and return a
            single value.
        frames : int or iterable of ints, optional
            Profiles to include by frame. If ``None``, includes
            all frames.
        bins : int or iterable of edges, optional
            If bins is an int, it defines the number of equal-width bins in the given
            range. If bins is a sequence, it defines a monotonically increasing array of
            bin edges, including the rightmost edge, allowing for non-uniform bin widths.
        range : (float, float), optional
            The lower and upper range of the bins.
            If not provided, ``range`` is simply ``(a.min(), a.max())``,
            where ``a`` is the array of reaction coordinates.
            Values outside the range are ignored. The first element of the range must be
            less than or equal to the second.


        Returns
        -------
        array of floats
            histogrammed, aggregate value of radii
        array of (float, float)
            Edges of each bin
        """
        binned, bins = self.bin_radii(frames=frames, bins=bins, range=range)
        return np.array(list(map(aggregator, binned))), bins

    def plot_mean_profile(self, bins=100, range=None,
                          frames=None, color='blue',
                          linestyle='-', ax=None,
                          xlabel='Frame', fill_alpha=0.3,
                          n_std=1, legend=True,
                          legend_loc='best',
                          **kwargs):
        """Collects the pore radii into bins by reaction coordinate.

        Parameters
        ----------
        frames : int or iterable of ints, optional
            Profiles to include by frame. If ``None``, includes
            all frames.
        bins : int or iterable of edges, optional
            If bins is an int, it defines the number of equal-width bins in the given
            range. If bins is a sequence, it defines a monotonically increasing array of
            bin edges, including the rightmost edge, allowing for non-uniform bin widths.
        range : (float, float), optional
            The lower and upper range of the bins.
            If not provided, ``range`` is simply ``(a.min(), a.max())``,
            where ``a`` is the array of reaction coordinates.
            Values outside the range are ignored. The first element of the range must be
            less than or equal to the second.
        color : str or array-like, optional
            Color for the plot.
        linestyle : str or array-like, optional
            Line style for the plot.
        ax : :class:`matplotlib.axes.Axes`
            If no `ax` is supplied or set to ``None`` then the plot will
            be added to the current active axes.
        xlabel : str, optional
            X-axis label.
        fill_alpha : float, optional
            Opacity of filled standard deviation area
        n_std : int, optional
            Number of standard deviations from the mean to fill between.
        legend : bool, optional
            Whether to plot a legend.
        legend_loc : str, optional
            Location of legend.
        **kwargs :
            All other `kwargs` are passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        ax : :class:`~matplotlib.axes.Axes`
             Axes with the plot, either `ax` or the current axes.

        """

        binned, bins = self.bin_radii(frames=frames, bins=bins, range=range)
        mean = np.array(list(map(np.mean, binned)))
        midpoints = 0.5 * (bins[1:] + bins[:-1])

        fig, ax = plt.subplots()
        if n_std:
            std = np.array(list(map(np.std, binned)))
            ax.fill_between(midpoints, mean-(n_std*std), mean+(n_std*std),
                            color=color, alpha=fill_alpha,
                            label='{} std'.format(n_std))
        ax.plot(midpoints, mean, color=color,
                linestyle=linestyle, label='mean', **kwargs)
        ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
        ax.set_ylabel(r"HOLE radius $R$ ($\AA$)")
        if legend:
            ax.legend(loc=legend_loc)
        return ax

    def plot3D_order_parameters(self, order_parameters,
                                frames=None,
                                color=None,
                                cmap='viridis',
                                linestyle='-', ax=None,
                                r_max=None,
                                ylabel=r'Order parameter',
                                **kwargs):
        r"""Plot HOLE radii over order parameters as a 3D graph.

        Lines are colored according to the specified ``color`` or
        drawn from the color map ``cmap``. One line is
        plotted for each trajectory frame.

        Parameters
        ----------
        order_parameters : array-like or string
            Sequence or text file containing order parameters(float
            numbers) corresponding to the frames in the trajectory. Must
            be same length as trajectory.
        frames : array-like, optional
            Frames to plot. If ``None``, plots all of them.
        color : str or array-like, optional
            Color or colors for the plot. If ``None``, colors are
            drawn from ``cmap``.
        cmap : str, optional
            color map to make colors for the plot if ``color`` is
            not given. Names should be from the ``matplotlib.pyplot.cm``
            module.
        linestyle : str or array-like, optional
            Line style for the plot.
        ax : : class: `matplotlib.axes.Axes`
            If no `ax` is supplied or set to ``None`` then the plot will
            be added to the current active axes.
        r_max : float, optional
            only display radii up to ``r_max``. If ``None``, all radii are
            plotted.
        ylabel : str, optional
            Y-axis label.
        **kwargs :
            All other `kwargs` are passed to: func: `matplotlib.pyplot.plot`.

        Returns
        -------
        ax: : class: `~mpl_toolkits.mplot3d.Axes3D`
             Axes with the plot, either `ax` or the current axes.

        """
        op_profiles = self.over_order_parameters(order_parameters,
                                                 frames=frames)

        from mpl_toolkits.mplot3d import Axes3D

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        ocl = self._process_plot_kwargs(frames=list(op_profiles.keys()),
                                        color=color, cmap=cmap,
                                        linestyle=linestyle)

        for op, c, ls in zip(*ocl):
            profile = op_profiles[op]
            if r_max is None:
                radius = profile.radius
                rxn_coord = profile.rxn_coord
            else:
                # does not seem to work with masked arrays but with nan hack!
                # http://stackoverflow.com/questions/4913306/python-matplotlib-mplot3d-how-do-i-set-a-maximum-value-for-the-z-axis
                rxn_coord = profile.rxn_coord
                radius = profile.radius.copy()
                radius[radius > r_max] = np.nan
            ax.plot(rxn_coord, op*np.ones_like(rxn_coord), radius,
                    color=c, linestyle=ls, zorder=int(-op), label=str(op),
                    **kwargs)

        ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
        ax.set_ylabel(ylabel)
        ax.set_zlabel(r"HOLE radius $R$ ($\AA$)")
        plt.tight_layout()
        return ax
