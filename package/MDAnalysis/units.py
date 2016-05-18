# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

r"""
Constants and unit conversion --- :mod:`MDAnalysis.units`
===============================================================

The base units of MDAnalysis trajectories are the **Å** (**ångström**) for
**length** and **ps** (**pico second**) for **time**. By default, all positions
are in Å and all times are in ps, regardless of how the MD code stored
trajectory data. By default, MDAnalysis converts automatically to the
MDAnalysis units when reading trajectories and converts back when writing. This
makes it possible to write scripts that can be agnostic of the specifics of how
a particular MD code stores trajectory data. Other base units are listed in the
table on :ref:`table-baseunits`.

.. _table-baseunits:

.. Table:: Base units in MDAnalysis

   =========== ============== ===============================================
   quantity    unit            SI units
   =========== ============== ===============================================
   length       Å              :math:`10^{-10}` m
   time         ps             :math:`10^{-12}` s
   energy       kJ/mol         :math:`1.66053892103219 \times 10^{-21}` J
   charge       :math:`e`      :math:`1.602176565 \times 10^{-19}` As
   force        kJ/(mol·Å)     :math:`1.66053892103219 \times 10^{-11}` J/m
   speed        Å/ps           :math:`100` m/s
   =========== ============== ===============================================

Implementation notes
--------------------

All conversions with :func:`convert` are carried out in a simple fashion: the
conversion factor :math:`f_{b,b'}` from the base unit :math:`b` to another unit
:math:`b'` is precomputed and stored (see :ref:`Data`). The numerical value of
a quantity in unit :math:`b` is :math:`X/b` (e.g. for :math:`X =
1.23\,\mathrm{ps}` the numerical value is :math:`X/\mathrm{ps} =
1.23`). [#funits]_

The new numerical value :math:`X'/b'` of the quantity (in units of :math:`b'`)
is then

.. math::

   X'/b' = f_{b,b'} X/b

The function :func:`get_conversion_factor` returns the appropriate factor
:math:`f_{b,b'}`.

Conversion between different units is always carried out via the base unit as
an intermediate step::

    x is in u1: from u1 to b:  x'  = x  / factor[u1]
                from b  to u2: x'' = x' * factor[u2]
    so f[u1,u2] = factor[u2]/factor[u1]


Conversions
~~~~~~~~~~~

Examples for how to calculate some of the conversion factors that are
hard-coded in :mod:`~MDAnalysis.units` (see :ref:`Data`).

density:
  Base unit is :math:`\mathrm{Å}^{-3}`::

     n/x = n/A**3 * densityUnit_factor[x]

  Example for how to calculate the conversion factor
  :math:`f_{\mathrm{Å}^{-3},\mathrm{nm}^{-3}}` from :math:`\mathrm{Å^-3}` to
  :math:`\mathrm{nm}^{-3}`:

  .. math::

     f_{\mathrm{Å}^{-3},\mathrm{nm}^{-3}}
           = \frac{1\,\mathrm{nm}^{-3}}{1\,\mathrm{Å}^{-3}}
           = \frac{(10\,\mathrm{Å})^{-3}}{1\,\mathrm{Å}^{-3}}
           = 10^{-3}

concentration:
  Example for how to convert the conversion factor to Molar (mol/l)::

     factor = 1 A**-3 / (N_Avogadro * (10**-9 dm)**-3)

  relative to a density rho0 in g/cm^3::

    M(H2O) = 18 g/mol   Molar mass of water

    factor = 1/(1e-24 * N_Avogadro / M(H2O))

  from `rho/rho0 = n/(N_A * M**-1) / rho0`  where `[n] = 1/Volume`, `[rho] = mass/Volume`


.. SeeAlso:: Maybe we should simply use Quantities_ or :mod:`scipy.constants`?

.. _Quantities: http://packages.python.org/quantities/

Functions
---------

.. autofunction:: get_conversion_factor
.. autofunction:: convert

.. _Data:

Data
----

.. autodata:: constants
.. autodata:: lengthUnit_factor
.. autodata:: N_Avogadro
.. autodata:: water
.. autodata:: densityUnit_factor
.. autodata:: timeUnit_factor
.. autodata:: speedUnit_factor
.. autodata:: forceUnit_factor
.. autodata:: chargeUnit_factor
.. autodata:: conversion_factor
.. autodata:: unit_types


References and footnotes
------------------------

.. [Jorgensen1998]  W. Jorgensen, C. Jenson, J Comp Chem 19 (1998), 1179-1186

.. _AKMA: http://www.charmm.org/documentation/c37b1/usage.html#%20AKMA
.. _electron charge: http://physics.nist.gov/cgi-bin/cuu/Value?e
.. _`Avogadro's constant`: http://physics.nist.gov/cgi-bin/cuu/Value?na

.. Rubric:: Footnotes

.. [#funits] One can also consider the conversion factor to carry
   units :math:`b'/b`, in which case the conversion formula would
   become

   .. math::

      X' = f_{b,b'} X

"""

from __future__ import unicode_literals

#: `Avogadro's constant`_ in mol**-1.
#:
#: .. deprecated:: 0.9.0
#:    Use *N_Avogadro* in dict :data:`constants` instead.
#:
#: .. versionchanged:: 0.9.0
#:    Use CODATA 2010 value, which differs from the previously used value
#:    6.02214179e+23 mol**-1 by -5.00000000e+16 mol**-1.
N_Avogadro = 6.02214129e+23  # mol**-1

#: Values of physical constants are taken from `CODATA 2010 at NIST`_. The
#: thermochemical calorie is defined in the `ISO 80000-5:2007`_ standard
#: and is also listed in the `NIST Guide to SI: Appendix B.8: Factors for Units`_.
#:
#: .. _`CODATA 2010 at NIST`:
#:    http://physics.nist.gov/cuu/Constants/
#: .. _`ISO 80000-5:2007`:
#:    http://www.iso.org/iso/catalogue_detail?csnumber=31890
#: .. _`NIST Guide to SI: Appendix B.8: Factors for Units`:
#:    http://physics.nist.gov/Pubs/SP811/appenB8.html#C
#:
#: .. versionadded:: 0.9.0
#
# NOTE: Whenever a constant is added to the constants dict, you also
#       MUST add an appropriate entry to
#       test_units:TestConstants.constants_reference !
constants = {
    'N_Avogadro': 6.02214129e+23,          # mol**-1
    'elementary_charge': 1.602176565e-19,  # As
    'calorie': 4.184,                      # J
}

#: The basic unit of *length* in MDAnalysis is the Angstrom.
#: Conversion factors between the base unit and other lengthUnits *x* are stored.
#: Conversions follow `L/x = L/Angstrom * lengthUnit_factor[x]`.
#: *x* can be *nm*/*nanometer* or *fm*.
lengthUnit_factor = {
    'Angstrom': 1.0, 'A': 1.0, 'angstrom': 1.0,
    u'\u212b': 1.0,   # Unicode and UTF-8 encoded symbol for angstroms
    'nm': 1.0 / 10, 'nanometer': 1.0 / 10,
    'pm': 1e2, 'picometer': 1e2,
    'fm': 1e5, 'femtometer': 1e5,
}


#: water density values at T=298K, P=1atm [Jorgensen1998]_
#:  ======== =========
#:  model    g cm**-3
#:  ======== =========
#:    SPC     0.985(1)
#:    TIP3P   1.002(1)
#:    TIP4P   1.001(1)
#:    exp     0.997
#:  ======== =========
#:
#: and molar mass 18.016 g mol**-1.
water = {
    'exp': 0.997, 'SPC': 0.985, 'TIP3P': 1.002, 'TIP4P': 1.001,  # in g cm**-3
    'MolarMass': 18.016,  # in g mol**-1
}

#: The basic unit for *densities* is Angstroem**(-3), i.e.
#: the volume per molecule in A**3. Especially for water
#: it can be convenient to measure the density relative to bulk, and
#: hence a number of values are pre-stored in :data:`water`.
densityUnit_factor = {
    'Angstrom^{-3}': 1 / 1.0, 'A^{-3}': 1 / 1.0,
    '\u212b^{-3}': 1 / 1.0,
    'nm^{-3}': 1 / 1e-3, 'nanometer^{-3}': 1 / 1e-3,
    'Molar': 1 / (1e-27 * constants['N_Avogadro']),
    'SPC': 1 / (1e-24 * constants['N_Avogadro'] * water['SPC'] / water['MolarMass']),
    'TIP3P': 1 / (1e-24 * constants['N_Avogadro'] * water['TIP3P'] / water['MolarMass']),
    'TIP4P': 1 / (1e-24 * constants['N_Avogadro'] * water['TIP4P'] / water['MolarMass']),
    'water': 1 / (1e-24 * constants['N_Avogadro'] * water['exp'] / water['MolarMass']),
}


#: For *time*, the basic unit is ps; in particular CHARMM's
#: 1 AKMA_ time unit = 4.888821E-14 sec is supported.
timeUnit_factor = {
    'ps': 1.0, 'picosecond': 1.0,  # 1/1.0
    'fs': 1e3, 'femtosecond': 1e3,  # 1/1e-3,
    'ns': 1e-3, 'nanosecond': 1e-3,  # 1/1e3,
    'second': 1e-12, 'sec': 1e-12, 's': 1e-12,  # 1/1e12,
    'AKMA': 1 / 4.888821e-2,
}
# getting the factor f:  1200ps * f = 1.2 ns  ==> f = 1/1000 ns/ps

#: For *speed*, the basic unit is Angstrom/ps.
speedUnit_factor = {
    'Angstrom/ps': 1.0, 'A/ps': 1.0, '\u212b/ps': 1.0,
    'Angstrom/picosecond': 1.0,
    'angstrom/picosecond': 1.0,  # 1
    'Angstrom/AKMA': 4.888821e-2,
    'nm/ps': 0.1, 'nanometer/ps': 0.1, 'nanometer/picosecond': 0.1,  # 1/10
    'nm/ns': 0.1 / 1e-3,
    'pm/ps': 1e2,
    'm/s': 1e-10 / 1e-12,
}
# (TODO: build this combinatorically from lengthUnit and timeUnit)

#: *Energy* is measured in kJ/mol.
energyUnit_factor = {
    'kJ/mol': 1.0,
    'kcal/mol': 1/constants['calorie'],
    'J': 1e3/constants['N_Avogadro'],
    'eV': 1e3/(constants['N_Avogadro'] * constants['elementary_charge']),
    }

#: For *force* the basic unit is kJ/(mol*Angstrom).
forceUnit_factor = {
    'kJ/(mol*Angstrom)': 1.0, 'kJ/(mol*A)': 1.0,
    'kJ/(mol*\u212b)': 1.0,
    'kJ/(mol*nm)': 10.0,
    'Newton': 1e13/constants['N_Avogadro'],
    'N': 1e13/constants['N_Avogadro'],
    'J/m': 1e13/constants['N_Avogadro'],
    'kcal/(mol*Angstrom)': 1/constants['calorie'],
}
# (TODO: build this combinatorically from lengthUnit and energyUnit)

#: *Charge* is measured in multiples of the `electron charge`_ *e*, with the value
#: *elementary_charge* in :data:`constants`.
#: The `conversion factor to Amber charge units`_ is 18.2223.
#:
#: .. _`conversion factor to Amber charge units`: http://ambermd.org/formats.html#parm
#:
#: .. versionchanged:: 0.9.0
#:    Use CODATA 2010 value for *elementary charge*, which differs from the previously used value
#:    *e* =  1.602176487 x 10**(-19) C by 7.8000000e-27 C.
chargeUnit_factor = {
    'e': 1.0,
    'Amber': 18.2223,  # http://ambermd.org/formats.html#parm
    'C': constants['elementary_charge'], 'As': constants['elementary_charge'],
}

#: :data:`conversion_factor` is used by :func:`get_conversion_factor`:
#: Note: any observable with a unit (i.e. one with an entry in
#: the :attr:`unit` attribute) needs an entry in :data:`conversion_factor`
conversion_factor = {
    'length': lengthUnit_factor,
    'density': densityUnit_factor,
    'time': timeUnit_factor,
    'charge': chargeUnit_factor,
    'speed': speedUnit_factor,
    'force': forceUnit_factor,
    'energy': energyUnit_factor,
}

#: Generated lookup table (dict): returns the type of unit for a known input unit.
#: Note: Any unit must be *unique* because this dict is used to guess the
#: unit type.
unit_types = {}
for utype, ufactor in conversion_factor.items():
    for unit in ufactor.keys():
        assert not unit in unit_types  # see comment!
        unit_types[unit] = utype


def get_conversion_factor(unit_type, u1, u2):
    """generate the conversion factor u1 -> u2 by using the base unit as an intermediate

    f[u1 -> u2] = factor[u2]/factor[u1]

    Conversion of X (in u1) to X' (in u2):

        X' = conversion_factor * X
    """
    # x is in u1: from u1 to b:  x'  = x  / factor[u1]
    #             from b  to u2: x'' = x' * factor[u2]
    # so f[u1,u2] = factor[u2]/factor[u1]
    return conversion_factor[unit_type][u2] / conversion_factor[unit_type][u1]


def convert(x, u1, u2):
    """Convert value *x* in unit *u1* to new value in *u2*.

    :Returns: Converted value.

    :Raises: :Exc:`ValueError` if the units are not known or if
             one attempts to convert between incompatible
             units.
    """
    try:
        ut1 = unit_types[u1]
    except KeyError:
        raise ValueError("unit '{0}' not recognized.\n"
                          "It must be one of {1}.".format(u1, ", ".join(unit_types)))
    try:
        ut2 = unit_types[u2]
    except KeyError:
        raise ValueError("unit '{0}' not recognized.\n"
                          "It must be one of {1}.".format(u2, ", ".join(unit_types)))
    if ut1 != ut2:
        raise ValueError("Cannot convert between unit types "
                         "{0} --> {1}".format(u1, u2))
    return x * get_conversion_factor(ut1, u1, u2)
