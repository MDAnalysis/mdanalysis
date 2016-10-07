# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2014 Naveen Michaud-Agrawal,
# Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see AUTHORS for the full list)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Neutron reflectivity analysis --- :mod:`MDAnalysis.analysis.reflectivity`
==========================================================


"""
import MDAnalysis
from MDAnalysis.lib.log import ProgressMeter
from MDAnalysis.core.topologyattrs import Isotopes
from MDAnalysis.topology import tables, guessers
from MDAnalysis.topology.tables import kv2dict


TABLE_TYPE2ISOTOPE = """
CH3 CH3
CH2 CH2
CH1 CH1
CH  CH
CD3 CD3
CD2 CD2
CD  CD
C   C
CO  C
D   D
H   H
HC  HC
OC  O
O2  O
OE  O
OE2  O
OA2  O
OA   O
"""

type2isotope_dict = kv2dict(TABLE_TYPE2ISOTOPE)


#Varley F. Sears (1992): Neutron scattering lengths and cross sections, Neutron News, 3:3, 26-7
#mass (amu)
#abundance (fraction)
#nucleus spin
#scattering lengths (fm)
#
#element wt abundance   nuc.spin  b_c (real)  b_c (imag) b_i (real)  b_i (imag)
TABLE_ISOTOPEDATA = """
  H    1     0.99985     0.5      -3.7406       0.000    25.274       0.
  H    2     0.00015     1.0       6.6710       0.000    4.0400       0.
  H    3     0           0.5       4.7920       0.000   -1.0400       0.
  He   3     0.0000014   0.5       5.7400      -1.483   -2.5000       2.568
  He   4     0.9999986   0.0       3.2600       0.000    0.0000       0.
  Li   6     0.075       1.0       2.0000      -0.261   -1.8900       0.26
  Li   7     0.925      -1.5       2.2200       0.000   -2.4900       0.
  Be   9     1          -1.5       7.7900       0.000    0.1200       0.
  B    10    0.2         3.0      -0.1         -1.066   -4.7000       1.231
  B    11    0.8        -1.5       6.6500       0.000    1.3000       0.
  C    12    0.989       0.0       6.6511       0.000    0.0000       0.
  C    13    0.011      -0.5       6.1900       0.000   -0.5200       0.
  N    14    0.9963      1.0       9.3700       0.000    2.0000       0.
  N    15    0.0037     -0.5       6.4400       0.000   -0.0200       0.
  O    16    0.99762     0.0       5.8030       0.000    0.0000       0.
  O    17    0.00038     2.5       5.7800       0.000    0.1800       0.
  O    18    0.002       0.0       5.8400       0.000    0.0000       0.
  F    19    1           0.5       5.6540       0.000   -0.8200       0.
  Ne   20    0.9051      0.0       4.6310       0.000    0.0000       0.
  Ne   21    0.0027      1.5       6.6600       0.000    0.6000       0.
  Ne   22    0.0922      0.0       3.8700       0.000    0.0000       0.
  Na   23    1           1.5       3.6300       0.000    3.5900       0.
  Mg   24    0.7899      0.0       5.6600       0.000    0.0000       0.
  Mg   25    0.1000      2.5       3.6200       0.000    1.4800       0.
  Mg   26    0.1101      0.0       4.8900       0.000    0.0000       0.
  Al   27    1           2.5       3.4490       0.000    0.2560       0.
  Si   28    0.9223      0.0       4.1070       0.000    0.0000       0.
  Si   29    0.0467      0.5       4.7000       0.000    0.0900       0.
  Si   30    0.0310      0.0       4.5800       0.000    0.0000       0.
  P    31    1           0.5       5.1300       0.000    0.2000       0.
  S    32    0.9502      0.0       2.8040       0.000    0.0000       0.
  S    33    0.0075      1.5       4.7400       0.000    1.5000       0.
  S    34    0.0421      0.0       3.4800       0.000    0.0000       0.
  S    36    0.0002      0.0       3.0000       0.000    0.0000       0.
  Cl   35    0.7577      1.5       11.650       0.000    6.1000       0.
  Cl   37    0.2423      1.5       3.0800       0.000    0.1000       0.
  Ar   36    0.00337     0.0       24.900       0.000    0.0000       0.
  Ar   38    0.00063     0.0       3.5000       0.000    0.0000       0.
  Ar   40    0.99600     0.0       1.8300       0.000    0.0000       0.
  K    39    0.93258     1.5       3.7400       0.000    1.4000       0.
  K    40    0.00012    -4.0       3.0000       0.000    0.0000       0.
  K    41    0.06730     1.5       2.6900       0.000    1.5000       0.
  Ca   40    0.96941     1.5       4.8000       0.000    0.0000       0.
  Ca   42    0.00647     0.0       3.3600       0.000    0.0000       0.
  Ca   43    0.00135     3.5      -1.5600       0.000    0.0000       0.
  Ca   44    0.02086     0.0       1.4200       0.000    0.0000       0.
  Ca   46    0.00004     0.0       3.6000       0.000    0.0000       0.
  Ca   48    0.00187     0.0       0.3900       0.000    0.0000       0.
"""


def guess_isotopes(u):
    """Add isotope topology attribute to the atoms in a universe and make
    guesses about what the isotopic composition of each atom should be.

    Parameters
    ----------
    u : MDAnalysis.Universe
        A universe instance to populate with isotope attributes

    Remarks
    -------
    After guessing, it is highly recommended to check the guesses and fix
    any errors. For specific isotopes of hydrogen, the following formula
    symbols should be used:

    - Pr = protium
    - D = deuterium
    - T = tritium
    
    Where H is used, the natural isotopic composition will be assumed.
    Example:
        CH4 = methane
        CH3D = methane with natural isotopic hydrogen except for one
        which is deuterium
        CPr4 = methane with no deuterium (all protium)

    """

    if 'isotopes' not in dir(u.atoms):
        u.add_TopologyAttr(Isotopes(np.empty(len(u.atoms), dtype='S10')))

    for atom in u.atoms:
        try:
            atom.isotope = type2isotope_dict[atom.type]
        except KeyError:
            try:
                atom.isotope = guessers.guess_atom_element(atom.type)
            except AttributeError:
                atom.isotope = guessers.guess_atom_element(atom.name)


def parse_isotope(formula):
    """Convert a formula in the style C2H6O into a dictionary whose keys
    are elements and values are the number of each element.
    """
    letters = []
    digits = []
    element_num_dict = {}
    formula = formula.strip()+' '
    for char in formula:
        if char.isdigit():
            digits.append(char)
        else:
            if char != char.upper():
                letters.append(char)
            else:
                element = ''.join(letters)
                if element:
                    if digits:
                        num = int(''.join(digits))
                        digits = []
                    else:
                        num = 1
                    try:
                        element_num_dict[element] += num
                    except:
                        element_num_dict[element] = num
                letters = [char]
    return element_num_dict


def bc_atom(element=None, isotope_dict=None): 
    # calculate total
    # calculate coherent scattering length from an element and isotope dict
    lines = TABLE_ISOTOPEDATA.splitlines()
    for line in lines:
        words = line.split()
        if words and words[0] == element:
            if isotope_dict == None: abundance = float(words[2])
            else: abundance = isotope_dict[int(words[1])]
            bc_r = float(words[4])
            bc_i = float(words[5])
            try:
                scat_length += abundance*(bc_r + 1j*bc_i) 
            except NameError:
                scat_length = abundance*(bc_r + 1j*bc_i) 
    try:
        return scat_length
    except NameError:
        raise ValueError("Couldn't calculate bc for element{}".format(element))


def element_dict_to_scat_len(element_num_dict):
    sl = 0+0j
    for element in element_num_dict:
        element_name = element
        isotope_dict = None
        if element == 'Pr': 
            isotope_dict = {1:1,2:0,3:0}
            element_name = 'H'
        if element == 'D': 
            isotope_dict = {1:0,2:1,3:0}
            element_name = 'H'
        if element == 'T': 
            isotope_dict = {1:0,2:0,3:1}
            element_name = 'H'
        sl += element_num_dict[element]*bc_atom(element_name,isotope_dict)
    return sl


def isotope_to_scat_len(isotope):
    element_num_dict = parse_isotope(isotope)
    sl = element_dict_to_scat_len(element_num_dict)
    return sl


def linear_scattering_length_density(selection, bins=10, range=None,
                                     start=None, stop=None, step=None,
                                     interval=1, quiet=False):
    """
    selection : AtomGroup
        Selection to perform the analysis on
    bins
        `bins` parameter passed to numpy.histogram
    range
        `range` parameter passed to numpy.histogram
    start, stop, step
        Slice the trajectory as ``trajectory[start"stop:step]``; default
        is to read the whole trajectory.
    quiet
        Print status update to the screen for every *interval* frame? [``False``]
        - ``True``: no status updates when a new frame is processed
        - ``False``: status update every frame (including number of atoms
          processed, which is interesting with ``update_selection=True``)
    interval
        Show status update every *interval* frame [1]

    Remarks
    -------
        If the water should be null-reflecting, simply omit it from the
        atom selection provided.
    """
    if type(selection) == MDAnalysis.Universe:
        u = selection
        selection = u.atoms
    else:
        u = selection.universe
    isotope_to_scatlen_dict = {}
    for i in set(selection.isotopes):
        isotope_to_scatlen_dict[i] = isotope_to_scat_len(i)
    scatlens = np.array([isotope_to_scatlen_dict[i] for i in selection.isotopes])

    pm = ProgressMeter(u.trajectory.n_frames, interval=interval, quiet=quiet,
                       format="Computing SLD profile for "
                       "%(n_atoms)6d atoms in frame "
                       "%(step)5d/%(numsteps)d  [%(percentage)5.1f%%]\r")
    if range is None:
        range = (0., u.dimensions[2])
    Lz = range[1] - range[0]
    n_frames = 0.
    for ts in u.trajectory[start:stop:step]:
        pm.echo(ts.frame, n_atoms=len(selection))
        V = u.dimensions[0]*u.dimensions[1]*Lz
        # TODO PBC wrap? or let user do that manually?
        sld_real, edges = np.histogram(selection.positions[:, 2], bins=bins,
                                        range=range,
                                        weights=np.real(scatlens))
        sld_imag, edges = np.histogram(selection.positions[:, 2], bins=bins,
                                        range=range,
                                        weights=np.imag(scatlens))
        sld_complex = sld_real+1j*sld_imag
        try:
            sld_sum += sld_complex/V*len(sld_complex)
        except NameError:
            sld_sum = sld_complex/V*len(sld_complex)
        n_frames += 1.
    return sld_sum*1e-5/n_frames # fm to angstroms



