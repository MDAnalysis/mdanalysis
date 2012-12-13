# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
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
Constants and unit conversion --- :mod:`MDAnalysis.core.units`
===============================================================

The base units of MDAnalysis are *Angstrom* for *length* (1 Angstrom =
0.1 nm = 10^-10 m) and *ps* (pico second) for *time* (1 ps = 10^-12
sec). For *force* we adopted kJ/(mol*Angstrom).

All conversions: the conversion factor f to a unit b' for a quantity X
(whose numeric value relative to the base unit b is stored in the
program) is a quantity with unit `b'/b`. In the dictionaries below only
the numeric value `f(b->b')` is stored::

  X/b' = f(b->b') * X/b

:func:`get_conversion_factor` returns the appropriate factor f(b->b').

Conversion is done via the base units::

    x is in u1: from u1 to b:  x'  = x  / factor[u1]
                from b  to u2: x'' = x' * factor[u2]
    so f[u1,u2] = factor[u2]/factor[u1]


Conversions
-----------

density conversion factor. Base unit is A**-3::

   n/x = n/A**3 * densityUnit_factor[x]

nm::

   f = 1 A^-3/1 nm^-3 = 1/(10A)^-3 = 1/1000

Molar::

   factor = 1 A**-3 / (N_Avogadro * (10**-9 dm)**-3)

relative to a density rho0 in g/cm^3::

    M(H2O) = 18 g/mol   Molar mass of water

    factor = 1/(1e-24 * N_Avogadro / M(H2O))

from `rho/rho0 = n/(N_A * M**-1) / rho0`  where `[n] = 1/Volume`, `[rho] = mass/Volume`


.. SeeAlso:: Maybe we should simply use Quantities_?

.. _Quantities: http://packages.python.org/quantities/

Functions
---------

.. autofunction:: get_conversion_factor
.. autofunction:: convert

Data
----

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


References
----------

.. [Jorgensen1998]  W. Jorgensen, C. Jenson, J Comp Chem 19 (1998), 1179-1186

.. _AKMA: http://www.charmm.org/html/documentation/c36b1/usage.html#%20AKMA
.. _electron charge: http://physics.nist.gov/cgi-bin/cuu/Value?e
.. _`Avogadro's constant`: http://physics.nist.gov/cgi-bin/cuu/Value?na

"""

#: The basic unit of *length* in MDAnalysis is the Angstrom.
#: Conversion factors between the base unit and other lengthUnits *x* are stored.
#: Conversions follow `L/x = L/Angstrom * lengthUnit_factor[x]`.
#: *x* can be *nm*/*nanometer* or *fm*.
lengthUnit_factor = {'Angstrom': 1.0, 'A': 1.0, 'Å': 1.0, 'angstrom': 1.0,
                     'nm': 1.0/10, 'nanometer': 1.0/10,
                     'pm': 1e2, 'picometer': 1e2,
                     'fm': 1e5, 'femtometer': 1e5,
                     }


#: `Avogadro's constant`_ in mol**-1,
N_Avogadro = 6.02214179e+23  # mol**-1

#: water density values ay 1179: T=298K, P=1atm [Jorgensen1998]_
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
water = {'exp':0.997, 'SPC':0.985, 'TIP3P':1.002, 'TIP4P':1.001,  # in g cm**-3
         'MolarMass': 18.016,                                     # in g mol**-1
         }

#: The basic unit for *densities* is Angstroem**(-3), i.e.
#: the volume per molecule in A**3. Especially for water
#: it can be convenient to measure the density relative to bulk, and
#: hence a number of values are pre-stored in :data:`water`.
densityUnit_factor = {
    'Angstrom^{-3}': 1/1.0, 'A^{-3}': 1/1.0, 'Å^{-3}': 1/1.0,
    'nm^{-3}': 1/1e-3, 'nanometer^{-3}': 1/1e-3,
    'Molar': 1/(1e-27*N_Avogadro),
    'SPC':   1/(1e-24*N_Avogadro*water['SPC']  / water['MolarMass']),
    'TIP3P': 1/(1e-24*N_Avogadro*water['TIP3P']/ water['MolarMass']),
    'TIP4P': 1/(1e-24*N_Avogadro*water['TIP4P']/ water['MolarMass']),
    'water': 1/(1e-24*N_Avogadro*water['exp']  / water['MolarMass']),
    }


#: For *time*, the basic unit is ps; in particular CHARMM's
#: 1 AKMA_ time unit = 4.888821E-14 sec is supported.
timeUnit_factor = {'ps': 1.0, 'picosecond': 1.0,    # 1/1.0
                   'fs': 1e3, 'femtosecond': 1e3,   # 1/1e-3,
                   'ns': 1e-3, 'nanosecond': 1e-3,  # 1/1e3,
                   'second': 1e-12, 'sec':  1e-12, 's':  1e-12, # 1/1e12,
                   'AKMA': 1/4.888821e-2,
                   }
# getting the factor f:  1200ps * f = 1.2 ns  ==> f = 1/1000 ns/ps

#: For *speed*, the basic unit is Angstrom/ps.
speedUnit_factor = {'Angstrom/ps': 1.0, 'A/ps': 1.0, 'Å/ps': 1.0, 'Angstrom/picosecond': 1.0,  'angstrom/picosecond': 1.0, # 1
                    'Angstrom/AKMA': 4.888821e-2,
                    'nm/ps': 0.1, 'nanometer/ps': 0.1, 'nanometer/picosecond': 0.1,       # 1/10
                    'nm/ns': 0.1/1e-3,
                    'pm/ps': 1e2,
                    'm/s': 1e-10/1e-12,
                    }
# (TODO: build this combinatorically from lengthUnit and timeUnit)

#: For *force* the basic unit is kJ/(mol*Angstrom).
forceUnit_factor = {'kJ/(mol*Angstrom)': 1.0, 'kJ/(mol*A)': 1.0, 'kJ/(mol*Å)': 1.0,
                    'kJ/(mol*nm)': 10.0,
                    }
# (TODO: build this combinatorically from lengthUnit and ... a new energyUnit)

#: *Charge* is measured in multiples of the `electron charge`_
#: *e* =  1.602176487 x 10**(-19) C.
chargeUnit_factor = {'e': 1.0,
                     'Amber': 18.2223,  # http://ambermd.org/formats.html#parm
                     'C': 1.602176487e-19, 'As': 1.602176487e-19,
                     }

#: :data:`conversion_factor` is used by :func:`get_conversion_factor`:
#: Note: any observable with a unit (i.e. one with an entry in
#: the :attr:`unit` attribute) needs an entry in :data:`conversion_factor`
conversion_factor = {'length': lengthUnit_factor,
                     'density': densityUnit_factor,
                     'time': timeUnit_factor,
                     'charge': chargeUnit_factor,
                     'speed': speedUnit_factor,
                     'force': forceUnit_factor,
                     }

#: Generated lookup table (dict): returns the type of unit for a known input unit.
#: Note: Any unit must be *unique* because this dict is used to guess the
#: unit type.
unit_types = {}
for utype,ufactor in conversion_factor.items():
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
    """Convert value in unit *u1* to *u2*."""
    try:
        ut1 = unit_types[u1]
        ut2 = unit_types[u2]
    except KeyError:
        raise ValueError("units must be one of %r, not %r or %r" %
                         (unit_types.keys(), u1, u2))
    if ut1 != ut2:
        raise ValueError("Cannot convert between unit types %(ut1)s --> %(ut2)s" %
                         vars())
    return x * get_conversion_factor(ut1, u1, u2)


