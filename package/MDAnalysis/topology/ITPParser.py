# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

r"""
ITP topology parser
===================

Reads a GROMACS ITP_ or TOP_ file to build the system. The topology will
contain atom IDs, segids, residue IDs, residue names, atom names, atom types,
charges, chargegroups, masses, moltypes, and molnums.
Any masses that are in the file will be read; any missing values will be guessed.
Bonds, angles, dihedrals and impropers are also read from the file.

If an ITP file is passed without a ``[ molecules ]`` directive, passing 
``infer_system=True`` (the default option) will create a Universe with 
1 molecule of each defined ``moleculetype``. 
If a ``[ molecules ]`` section is present, ``infer_system`` is ignored.

If files are included with the `#include` directive, they will also be read.
If they are not in the working directory, ITPParser will look for them in the
``include_dir`` directory. By default, this is set to
``include_dir="/usr/local/gromacs/share/gromacs/top/"``.
Variables can be defined with the `#define` directive in files, or by passing
in :ref:`keyword arguments <itp-define-kwargs>`.

Examples
--------

.. code-block:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import ITP_tip5p

    #  override charge of HW2 atom (defined in file as HW2_CHARGE)
    u = mda.Universe(ITP_tip5p, HW2_CHARGE=2, infer_system=True)

.. note::

    AMBER also uses topology files with the .top extension. To use ITPParser
    to read GROMACS top files, pass ``topology_format='ITP'``.

    ::

        import MDAnalysis as mda
        u = mda.Universe('topol.top', topology_format='ITP')


.. _itp-define-kwargs:

Preprocessor variables
----------------------

ITP files are often defined with lines that depend on 
whether a keyword flag is given. For example, this modified TIP5P water file:

.. code-block:: none

    [ moleculetype ]
    ; molname       nrexcl
    SOL             2

    #ifndef HW1_CHARGE
        #define HW1_CHARGE 0.241
    #endif

    #define HW2_CHARGE 0.241

    [ atoms ]
    ; id    at type res nr  residu name     at name         cg nr   charge
    1       opls_118     1       SOL              OW             1       0
    2       opls_119     1       SOL             HW1             1       HW1_CHARGE
    3       opls_119     1       SOL             HW2             1       HW2_CHARGE
    4       opls_120     1       SOL             LP1             1      -0.241
    5       opls_120     1       SOL             LP2             1      -0.241
    #ifdef EXTRA_ATOMS  ; added for keyword tests
    6       opls_120     1       SOL             LP3             1      -0.241
    7       opls_120     1       SOL             LP4             1       0.241
    #endif


Define these preprocessor variables by passing keyword arguments. Any arguments that you 
pass in *override* any variables defined in the file. For example, the universe below 
will have charges of 3 for the HW1 and HW2 atoms::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import ITP_tip5p

    u = mda.Universe(ITP_tip5p, EXTRA_ATOMS=True, HW1_CHARGE=3, HW2_CHARGE=3)

These keyword variables are **case-sensitive**. Note that if you set keywords to
``False`` or ``None``, they will be treated as if they are not defined in #ifdef conditions.

For example, the universe below will only have 5 atoms. ::

    u = mda.Universe(ITP_tip5p, EXTRA_ATOMS=False)


.. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file

.. _TOP: http://manual.gromacs.org/current/reference-manual/file-formats.html#top

Classes
-------

.. autoclass:: ITPParser
   :members:
   :inherited-members:

"""
from collections import defaultdict
import os

import logging
import numpy as np

import warnings
from ..lib.util import openany
from ..guesser.tables import SYMB2Z
from ..guesser.default_guesser import DefaultGuesser
from .base import TopologyReaderBase, change_squash, reduce_singular
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Moltypes,
    Molnums,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
    Bonds,
    Angles,
    Dihedrals,
    Impropers,
    AtomAttr,
    Elements,
)
from ..core.topology import Topology


class Chargegroups(AtomAttr):
    """The charge group for each Atom"""

    attrname = "chargegroups"
    singular = "chargegroup"


class GmxTopIterator:
    """
    Iterate over the lines of a TOP/ITP file and its included files

    The iterator strips comments, deals with #ifdef and #ifndef conditions,
    and substitutes defined variables.

    Defined variables passed into ``__init__`` *override* anything defined in the file.
    """

    def __init__(self, path, include_dir, defines):
        self._original_defines = defines
        self.defines = dict(**defines)  # copy
        self.include_dir = include_dir
        self.file_stack = [path]
        self.starting_file = path

    @property
    def current_file(self):
        return self.file_stack[-1]

    def __iter__(self):
        for line in self.iter_from_file(self.starting_file):
            yield line

    def iter_from_file(self, path):
        found_path = self.find_path(path)
        with openany(found_path) as infile:
            self.file_stack.append(infile)

            for line in self.clean_file_lines(infile):
                if line.startswith("#include"):
                    inc = line.split(None, 1)[1][1:-1]
                    for line in self.iter_from_file(inc):
                        yield line
                elif line.startswith("#define"):
                    self.define(line)
                elif line.startswith("#if"):
                    self.do_if(line, infile)
                elif line.startswith("#else"):
                    self.skip_until_endif(infile)
                elif line.startswith("#"):  # ignore #if and others
                    pass
                elif line:
                    line = self.substitute_defined(line)
                    yield line
            self.file_stack.pop()

    def define(self, line):
        try:
            _, variable, value = line.split(None, 2)
        except ValueError:
            _, variable = line.split()
            value = True

        # kwargs overrides files
        if variable not in self._original_defines:
            self.defines[variable] = value

    def substitute_defined(self, line):
        split = line.split()
        for k, v in self.defines.items():
            if k in split:
                split[split.index(k)] = str(v)
        line = " ".join(split)
        return line

    def clean_file_lines(self, infile):
        for line in infile:
            line = line.split(";")[0].strip()  # ; is for comments
            yield line

    def do_if(self, line, infile):
        ifdef, variable = line.split()
        if ifdef == "#ifdef":
            if self.defines.get(variable) in (False, None):
                self.skip_until_else(infile)
        elif ifdef == "#ifndef":
            if self.defines.get(variable) not in (False, None):
                self.skip_until_else(infile)

    def skip_until_else(self, infile):
        """Skip lines until #if condition ends"""
        for line in self.clean_file_lines(infile):
            if line.startswith("#if"):
                self.skip_until_endif(infile)
            elif line.startswith("#endif") or line.startswith("#else"):
                break
        else:
            raise IOError("Missing #endif in {}".format(self.current_file))

    def skip_until_endif(self, infile):
        """Skip lines until #endif"""
        for line in self.clean_file_lines(infile):
            if line.startswith("#if"):
                self.skip_until_endif(infile)
            elif line.startswith("#endif"):
                break
        else:
            raise IOError("Missing #endif in {}".format(self.current_file))

    def find_path(self, path):
        try:
            # in case of TextIOWrapper
            current_file = self.current_file.name
        except AttributeError:
            current_file = self.current_file

        try:
            path = os.path.abspath(path.name)
        except AttributeError:
            pass
        current_dir = os.path.dirname(current_file)
        dir_path = os.path.join(current_dir, path)
        if os.path.exists(dir_path):
            return dir_path
        include_path = os.path.join(self.include_dir, path)
        if os.path.exists(include_path):
            return include_path
        raise IOError("Could not find {}".format(path))


class Molecule:
    """Store moleculetype-specific attributes"""

    def __init__(self, name):
        self.name = name
        self.ids = []
        self.types = []
        self.resids = []
        self.resnames = []
        self.names = []
        self.chargegroups = []
        self.charges = []
        self.masses = []

        self.bonds = defaultdict(list)
        self.angles = defaultdict(list)
        self.dihedrals = defaultdict(list)
        self.impropers = defaultdict(list)

        self.parsers = {
            "atoms": self.parse_atoms,
            "bonds": self.parse_bonds,
            "angles": self.parse_angles,
            "dihedrals": self.parse_dihedrals,
            "constraints": self.parse_constraints,
            "settles": self.parse_settles,
        }

        self.resolved_residue_attrs = False

    @property
    def atom_order(self):
        return [
            self.ids,
            self.types,
            self.resids,
            self.resnames,
            self.names,
            self.chargegroups,
            self.charges,
            self.masses,
        ]

    @property
    def params(self):
        return [self.bonds, self.angles, self.dihedrals, self.impropers]

    def parse_atoms(self, line):
        values = line.split()
        for lst in self.atom_order:
            try:
                lst.append(values.pop(0))
            except IndexError:  # ran out of values
                lst.append("")

    def parse_bonds(self, line):
        self.add_param(
            line,
            self.bonds,
            n_funct=2,
            funct_values=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
        )

    def parse_angles(self, line):
        self.add_param(
            line,
            self.angles,
            n_funct=3,
            funct_values=(1, 2, 3, 4, 5, 6, 8, 10),
        )

    def parse_dihedrals(self, line):
        dih = self.add_param(
            line,
            self.dihedrals,
            n_funct=4,
            funct_values=(1, 3, 5, 8, 9, 10, 11),
        )
        if not dih:
            self.add_param(
                line, self.impropers, n_funct=4, funct_values=(2, 4)
            )

    def parse_constraints(self, line):
        self.add_param(line, self.bonds, n_funct=2, funct_values=(1, 2))

    def parse_settles(self, line):
        # [ settles ] is a triangular constraint for
        # water molecules.
        # In ITP files this is defined with only the
        # oxygen atom index. The next two atoms are
        # assumed to be hydrogens. Unlike TPRParser,
        # the manual only lists this format (as of 2019).
        # These are treated as 2 bonds.
        # No angle component is included to avoid discrepancies
        # with water molecules loaded from different MD engines.
        oxygen, funct, doh, dhh = line.split()
        try:
            base = self.index_ids([oxygen])[0]
        except ValueError:
            pass
        else:
            self.bonds[(base, base + 1)].append("settles")
            self.bonds[(base, base + 2)].append("settles")

    def resolve_residue_attrs(self):
        """Figure out residue borders and assign moltypes and molnums"""
        resids = np.array(self.resids, dtype=np.int32)
        resnames = np.array(self.resnames, dtype=object)
        self.residx, (self.resids, resnames) = change_squash(
            (resids,), (resids, resnames)
        )
        self.resnames = list(resnames)
        self.moltypes = [self.name] * len(self.resids)
        self.molnums = np.array([1] * len(self.resids))

        self.resolved_residue_attrs = True

    def shift_indices(
        self, atomid=0, resid=0, molnum=0, cgnr=0, n_res=0, n_atoms=0
    ):
        """
        Get attributes ready for adding onto a larger topology.

        Shifts atom indices, residue indices, molnums, and chargegroup numbers.

        Returns
        -------
        atom_attrs: list of lists
            attributes in the [ atoms ] section

        new_params: list of dicts
            Bonds, angles, dihedrals, impropers as dicts of shape {indices: parameters}

        molnums: list
        moltypes: list
        residx: list
        """
        if not self.resolved_residue_attrs:
            self.resolve_residue_attrs()

        resids = list(np.array(self.resids) + resid)
        residx = list(np.array(self.residx) + n_res)
        molnums = list(np.array(self.molnums) + molnum)
        ids = list(np.array(self.ids, dtype=int) + atomid)

        try:
            cg = np.array(self.chargegroups, dtype=int)
        except ValueError:
            cg = np.arange(1, len(self.chargegroups) + 1)
        chargegroups = list(cg + cgnr)

        atom_order = [
            ids,
            self.types,
            resids,
            self.resnames,
            self.names,
            chargegroups,
            self.charges,
            self.masses,
        ]

        new_params = []
        for p in self.params:
            new = {}
            for indices, values in p.items():
                new[tuple(np.array(indices) + n_atoms)] = values
            new_params.append(new)

        return atom_order, new_params, molnums, self.moltypes, residx

    def add_param(self, line, container, n_funct=2, funct_values=None):
        """Add defined GROMACS directive lines, only if the funct is in ``funct_values``"""
        if funct_values is None:
            funct_values = []
        values = line.split()
        funct = int(values[n_funct])
        if funct in funct_values:
            try:
                ids = self.index_ids(values[:n_funct])
                container[ids].append(funct)
            except ValueError:
                pass
            return True
        else:
            return False

    def index_ids(self, values):
        """
        Get indices of atom ids (list of strings)
        """
        return tuple(map(self.ids.index, values))


class ITPParser(TopologyReaderBase):
    """Read topology information from a GROMACS ITP_ or TOP_ file.

    Creates a Topology with the following Attributes:
    - ids
    - names
    - types
    - masses
    - charges
    - chargegroups
    - resids
    - resnames
    - segids
    - moltypes
    - molnums
    - bonds
    - angles
    - dihedrals
    - impropers

    .. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file
    .. _TOP: http://manual.gromacs.org/current/reference-manual/file-formats.html#top

    .. note::

        By default, atomtypes and masses will be guessed on Universe creation
        if they are not read from the input file.
        This may change in release 3.0.
        See :ref:`Guessers` for more information.

    .. versionchanged:: 2.2.0
      no longer adds angles for water molecules with SETTLE constraint
    .. versionchanged:: 2.8.0
      Removed mass guessing (attributes guessing takes place now
      through universe.guess_TopologyAttrs() API).
      Added guessed elements if all elements are valid to preserve partial
      mass guessing behavior

    """

    format = "ITP"

    def parse(
        self,
        include_dir="/usr/local/gromacs/share/gromacs/top/",
        infer_system=True,
        **kwargs,
    ):
        """Parse ITP file into Topology

        Parameters
        ----------
        include_dir: str, optional
            A directory in which to look for other files included
            from the original file, if the files are not first found
            in the current directory.
            Default: "/usr/local/gromacs/share/gromacs/top/"

        infer_system: bool, optional (default True)
            If a ``[ molecules ]`` directive is not found within the the
            topology file, create a Topology with one of every
            ``[ moleculetype ]`` defined. If a ``[ molecules ]`` directive is
            found, this keyword is ignored.

        Returns
        -------
        MDAnalysis *Topology* object
        """

        self.atomtypes = {}
        self.molecules = {}
        self._molecules = []  # for order
        self.current_mol = None
        self.parser = self._pass
        self.system_molecules = []

        # Open and check itp validity
        with openany(self.filename) as itpfile:
            self.lines = GmxTopIterator(itpfile, include_dir, kwargs)
            for line in self.lines:
                if "[" in line and "]" in line:
                    section = line.split("[")[1].split("]")[0].strip()

                    if section == "atomtypes":
                        self.parser = self.parse_atomtypes

                    elif section == "moleculetype":
                        self.parser = self.parse_moleculetype

                    elif section == "molecules":
                        self.parser = self.parse_molecules

                    elif self.current_mol:
                        self.parser = self.current_mol.parsers.get(
                            section, self._pass
                        )

                    else:
                        self.parser = self._pass

                else:
                    self.parser(line)

        if not self.system_molecules and infer_system:
            self.system_molecules = [x.name for x in self._molecules]

        self.build_system()

        self.types = np.array(self.types)
        self.charges = np.array(self.charges, dtype=object)
        self.masses = np.array(self.masses, dtype=object)

        if not all(self.charges):
            empty = self.charges == ""
            self.charges[empty] = [
                (
                    self.atomtypes.get(x)["charge"]
                    if x in self.atomtypes.keys()
                    else ""
                )
                for x in self.types[empty]
            ]

        if not all(self.masses):
            empty = self.masses == ""
            self.masses[empty] = [
                (
                    self.atomtypes.get(x)["mass"]
                    if x in self.atomtypes.keys()
                    else ""
                )
                for x in self.types[empty]
            ]

        attrs = []
        # atom stuff
        for vals, Attr, dtype in (
            (self.ids, Atomids, np.int32),
            (self.types, Atomtypes, object),
            (self.names, Atomnames, object),
            (self.chargegroups, Chargegroups, np.int32),
            (self.charges, Charges, np.float32),
        ):
            if all(vals):
                attrs.append(Attr(np.array(vals, dtype=dtype)))

        if not all(self.masses):
            empty = self.masses == ""
            self.masses[empty] = Masses.missing_value_label

        attrs.append(
            Masses(np.array(self.masses, dtype=np.float64), guessed=False)
        )

        self.elements = DefaultGuesser(None).guess_types(self.types)
        if all(e.capitalize() in SYMB2Z for e in self.elements):
            attrs.append(
                Elements(np.array(self.elements, dtype=object), guessed=True)
            )
            warnings.warn(
                "The elements attribute has been populated by guessing "
                "elements from atom types. This behaviour has been "
                "temporarily added to the ITPParser as we transition "
                "to the new guessing API. "
                "This behavior will be removed in release 3.0. "
                "Please see issue #4698 for more information. ",
                DeprecationWarning,
            )
        else:
            warnings.warn(
                "Element information is missing, elements attribute "
                "will not be populated. If needed these can be "
                "guessed using universe.guess_TopologyAttrs("
                "to_guess=['elements'])."
            )

        # residue stuff
        resids = np.array(self.resids, dtype=np.int32)
        resnames = np.array(self.resnames, dtype=object)
        molnums = np.array(self.molnums, dtype=np.int32)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))
        attrs.append(Moltypes(np.array(self.moltypes, dtype=object)))
        attrs.append(Molnums(molnums))

        n_atoms = len(self.ids)
        n_residues = len(self.resids)
        n_segments = len(self.system_molecules)
        attrs.append(Segids(np.array(self.system_molecules, dtype=object)))
        segidx = molnums - 1

        top = Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=self.residx,
            residue_segindex=segidx,
        )

        # connectivity stuff
        for dct, Attr, attrname in (
            (self.bonds, Bonds, "bonds"),
            (self.angles, Angles, "angles"),
            (self.dihedrals, Dihedrals, "dihedrals"),
            (self.impropers, Impropers, "impropers"),
        ):
            if dct:
                indices, types = zip(*list(dct.items()))
            else:
                indices, types = [], []

            types = [reduce_singular(t) for t in types]

            tattr = Attr(indices, types=types)
            top.add_TopologyAttr(tattr)

        return top

    def _pass(self, line):
        pass

    def parse_atomtypes(self, line):
        keys = ["type_bonded", "atomic_number", "mass", "charge", "p_type"]
        fields = line.split()
        if len(fields[5]) == 1 and fields[5].isalpha():
            values = fields[1:6]
        elif len(fields[3]) == 1 and fields[3].isalpha():
            values = "", "", fields[1], fields[2], fields[3]
        elif len(fields[4]) == 1 and fields[4].isalpha():
            if fields[1][0].isalpha():
                values = fields[1], "", fields[2], fields[3], fields[4]
            else:
                values = "", fields[1], fields[2], fields[3], fields[4]
        self.atomtypes[fields[0]] = dict(zip(keys, values))

    def parse_moleculetype(self, line):
        name = line.split()[0]
        self.current_mol = self.molecules[name] = Molecule(name)
        self._molecules.append(self.current_mol)

    def parse_molecules(self, line):
        name, n_mol = line.split()
        self.system_molecules.extend([name] * int(n_mol))

    def build_system(self):
        self.ids = []
        self.types = []
        self.resids = []
        self.resnames = []
        self.names = []
        self.chargegroups = []
        self.charges = []
        self.masses = []
        self.moltypes = []
        self.molnums = []
        self.residx = []

        self.atom_order = [
            self.ids,
            self.types,
            self.resids,
            self.resnames,
            self.names,
            self.chargegroups,
            self.charges,
            self.masses,
        ]

        self.bonds = defaultdict(list)
        self.angles = defaultdict(list)
        self.dihedrals = defaultdict(list)
        self.impropers = defaultdict(list)

        self.params = [self.bonds, self.angles, self.dihedrals, self.impropers]

        for i, moltype in enumerate(self.system_molecules):
            mol = self.molecules[moltype]

            atomid = self.ids[-1] if self.ids else 0
            resid = self.resids[-1] if self.resids else 0
            cgnr = self.chargegroups[-1] if self.chargegroups else 0
            n_res = len(self.resids)
            n_atoms = len(self.ids)

            shifted = mol.shift_indices(
                atomid=atomid,
                resid=resid,
                n_res=n_res,
                cgnr=cgnr,
                molnum=i,
                n_atoms=n_atoms,
            )
            atom_order, params, molnums, moltypes, residx = shifted

            for system_attr, mol_attr in zip(self.atom_order, atom_order):
                system_attr.extend(mol_attr)

            self.moltypes.extend(moltypes)
            self.molnums.extend(molnums)
            self.residx.extend(residx)

            for system_param, mol_param in zip(self.params, params):
                system_param.update(mol_param)
