# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
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
#     Molecular Dynamics Simulations. J. Comput. Chem. (2011),
#     doi:10.1002/jcc.21787
#

# MDAnalysis -- nucleic acid analysis
# Copyright (c) 2011 Elizabeth Denning <denniej0@gmail.com>

"""
Generating nucleic acid based information from PDB file or trajectories --- :mod:`MDAnalysis.analysis.NuclInfo`
=============================================================================

:Author: Elizabeth Denning
:Year: 2011
:Copyright: GNU Public License v3

The module provides functions to represent nuclic acid data, in
particular backbone dihedrals, chi dihedrals, AS or CP phase angles,
     WC N1-N3 distances, C2-O2 distances, N6-O4 distances, O6-N4 distances

"""

from MDAnalysis import *
import numpy
from math import *

def wc_pair(i,bp):
        if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DC","DT","U","C","T","CYT","THY","URA"]:
		a1,a2 = "N3","N1"
	if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DG","DA","A","G","ADE","GUA"]:
		a1,a2 = "N1","N3"
	wc_dist = universe.selectAtoms(" (resid %s and name %s) or (resid %s and name %s) "%(i,a1,bp,a2))
	wc = norm(wc_dist[0].pos - wc_dist[1].pos)
	return wc

def minor_pair(i,bp):
        if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DC","DT","U","C","T","CYT","THY","URA"]:
	        a1,a2 = "O2","C2"
	if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DG","DA","A","G","ADE","GUA"]:
		a1,a2 = "C2","O2"
	c2o2_dist = universe.selectAtoms(" (resid %s and name %s) or (resid %s and name %s) "%(i,a1,bp,a2))
	c2o2 = norm(c2o2_dist[0].pos - c2o2_dist[1].pos)
	return c2o2

def major_pair(i,bp):
        if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DC","DG","C","G","CYT","GUA"]:
	        a1,a2 = "O6","N4"
	if universe.selectAtoms(" resid %s "%(i)).resnames()[0] in ["DT","DA","A","T","U","ADE","THY","URA"]: 
		a1,a2 = "N6","O4"               
	no_dist = universe.selectAtoms(" (resid %s or resid %s) and (name %s or name %s) "%(i,bp,a1,a2))
	major = norm(no_dist[0].pos - no_dist[1].pos)
	return major


def phase_cp(seg,i):
        atom1 = universe.selectAtoms(" atom %s %s O4\' "%(seg,i))
	atom2 = universe.selectAtoms(" atom %s %s C1\' "%(seg,i))
	atom3 = universe.selectAtoms(" atom %s %s C2\' "%(seg,i))
	atom4 = universe.selectAtoms(" atom %s %s C3\' "%(seg,i))
	atom5 = universe.selectAtoms(" atom %s %s C4\' "%(seg,i))

	data1 = atom1.coordinates()
	data2 = atom2.coordinates()
	data3 = atom3.coordinates()
	data4 = atom4.coordinates()
	data5 = atom5.coordinates()


	r0 = ( data1 + data2 + data3 + data4 + data5 ) * (1.0/ 5.0)
	r1 = data1 - r0
	r2 = data2 - r0
	r3 = data3 - r0
	r4 = data4 - r0
	r5 = data5 - r0

	R1 = ( ( r1 * sin(2 * 3.14159265 * 0.0 / 5.0) )\
		+ ( r2 * sin(2 * 3.14159265 * 1.0 / 5.0) )\
		+ ( r3 * sin(2 * 3.14159265 * 2.0 / 5.0) )\
		+ ( r4 * sin(2 * 3.14159265 * 3.0 / 5.0) )\
		+ ( r5 * sin(2 * 3.14159265 * 4.0 / 5.0) ) )

	R2 = ( ( r1 * cos(2 * 3.14159265 * 0.0 / 5.0) )\
		+ ( r2 * cos(2 * 3.14159265 * 1.0 / 5.0) )\
		+ ( r3 * cos(2 * 3.14159265 * 2.0 / 5.0) )\
		+ ( r4 * cos(2 * 3.14159265 * 3.0 / 5.0) )\
		+ ( r5 * cos(2 * 3.14159265 * 4.0 / 5.0) ) )

	x = numpy.cross(R1[0],R2[0])

	n = x / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )

	r1_d = numpy.dot(r1,n)
	r2_d = numpy.dot(r2,n)
	r3_d = numpy.dot(r3,n)
	r4_d = numpy.dot(r4,n)
	r5_d = numpy.dot(r5,n)

	D = ( (r1_d * sin(4 * 3.14159265 * 0.0 / 5.0))\
	        + (r2_d * sin(4 * 3.14159265 * 1.0 / 5.0))\
		+ (r3_d * sin(4 * 3.14159265 * 2.0 / 5.0))\
		+ (r4_d * sin(4 * 3.14159265 * 3.0 / 5.0))\
		+ (r5_d * sin(4 * 3.14159265 * 4.0 / 5.0)) ) * -1 * sqrt(2.0/5.0)

	C = ( (r1_d * cos(4 * 3.14159265 * 0.0 / 5.0))\
		+ (r2_d * cos(4 * 3.14159265 * 1.0 / 5.0))\
		+ (r3_d * cos(4 * 3.14159265 * 2.0 / 5.0))\
		+ (r4_d * cos(4 * 3.14159265 * 3.0 / 5.0))\
		+ (r5_d * cos(4 * 3.14159265 * 4.0 / 5.0)) ) * sqrt(2.0/5.0)
	
	phase_ang = ( atan2(D,C) + ( 3.14159265 / 2. ) ) * 180. / pi
	if phase_ang < 0:
		phase_ang = phase_ang+360
	else:
		phase_ang
	return phase_ang


def phase_as(seg,i):
        angle1 = universe.selectAtoms(" atom %s %s C1\' "%(seg,i)," atom %s %s C2\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s C4\' "%(seg,i))
	angle2 = universe.selectAtoms(" atom %s %s C2\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s O4\' "%(seg,i))
	angle3 = universe.selectAtoms(" atom %s %s C3\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i))
	angle4 = universe.selectAtoms(" atom %s %s C4\' "%(seg,i)," atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s C2\' "%(seg,i))
	angle5 = universe.selectAtoms(" atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s C2\' "%(seg,i)," atom %s %s C3\' "%(seg,i))

	data1 = angle1.dihedral()
	data2 = angle2.dihedral()
	data3 = angle3.dihedral()
	data4 = angle4.dihedral()
	data5 = angle5.dihedral()

	B = ( (data1 * sin(2 * 2 * 3.1416 * (1-1.)/5.))\
		+ (data2 * sin(2 * 2 * 3.1416 * (2-1.)/5.))\
		+ (data3 * sin(2 * 2 * 3.1416 * (3-1.)/5.))\
		+ (data4 * sin(2 * 2 * 3.1416 * (4-1.)/5.))\
		+ (data5 * sin(2 * 2 * 3.1416 * (5-1.)/5.)) ) * -2. / 5.

	A = ( (data1 * cos(2 * 2 * 3.1416 * (1-1.)/5.))\
		+ (data2 * cos(2 * 2 * 3.1416 * (2-1.)/5.))\
		+ (data3 * cos(2 * 2 * 3.1416 * (3-1.)/5.))\
		+ (data4 * cos(2 * 2 * 3.1416 * (4-1.)/5.))\
		+ (data5 * cos(2 * 2 * 3.1416 * (5-1.)/5.)) ) * 2. / 5.

	phase_ang = atan2(B,A) * 180. / pi
	if phase_ang < 0:
		phase_ang = phase_ang+360
	else:
		phase_ang
	return phase_ang

def tors(seg,i):	
	a = universe.selectAtoms(" atom %s %s O3\' "%(seg,i-1)," atom %s %s P  "%(seg,i)," atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i))
	b = universe.selectAtoms(" atom %s %s P    "%(seg,i)," atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i))
	g = universe.selectAtoms(" atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i))
	d = universe.selectAtoms(" atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i))
	e = universe.selectAtoms(" atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i)," atom %s %s P    "%(seg,i+1))
	z = universe.selectAtoms(" atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i)," atom %s %s P    "%(seg,i+1)," atom %s %s O5\' "%(seg,i+1))
	try:
		c = universe.selectAtoms(" atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s N1 "%(seg,i)," atom %s %s C2  "%(seg,i))
	except:
		c = universe.selectAtoms(" atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s N9 "%(seg,i)," atom %s %s C4  "%(seg,i))

	alpha = a.dihedral()
	beta = b.dihedral()
	gamma = g.dihedral()
	delta = d.dihedral()
	epsilon = e.dihedral()
	zeta = z.dihedral()
	chi = c.dihedral()

	if alpha < 0 :
		alpha = alpha+360
	if beta < 0 :
		beta = beta+360
	if gamma < 0 :
		gamma = gamma+360
	if epsilon < 0:
		epsilon = epsilon+360
	if zeta < 0 :
		zeta = zeta+360
	if chi < 0 :
		chi = chi+360
	return alpha, beta, gamma, delta, epsilon, zeta, chi

def tors_alpha(seg,i):
	a = universe.selectAtoms(" atom %s %s O3\' "%(seg,i-1)," atom %s %s P  "%(seg,i)," atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i))
	alpha = a.dihedral()
	if alpha < 0 :
	        alpha = alpha+360
	return alpha

def tors_beta(seg,i):
        b = universe.selectAtoms(" atom %s %s P    "%(seg,i)," atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i))
	beta = b.dihedral()
	if beta < 0 :
		beta = beta+360
	return beta

def tors_gamma(seg,i):
	g = universe.selectAtoms(" atom %s %s O5\' "%(seg,i)," atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i))
	gamma = g.dihedral()
	if gamma < 0:
		gamma = gamma + 360
	return gamma

def tors_delta(seg,i):
	d = universe.selectAtoms(" atom %s %s C5\' "%(seg,i)," atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i))
	delta = d.dihedral()
	if delta < 0:
		delta = delta + 360
	return delta

def tors_eps(seg,i):
	e = universe.selectAtoms(" atom %s %s C4\' "%(seg,i)," atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i)," atom %s %s P    "%(seg,i+1))
	epsilon = e.dihedral()
	if epsilon < 0:
		epsilon = epsilon + 360
	return epsilon

def tors_zeta(seg,i):
	z = universe.selectAtoms(" atom %s %s C3\' "%(seg,i)," atom %s %s O3\' "%(seg,i)," atom %s %s P    "%(seg,i+1)," atom %s %s O5\' "%(seg,i+1))
	zeta = z.dihedral()
	if zeta < 0:
		zeta = zeta + 360
	return zeta

def tors_chi(seg,i):
	try:
		c = universe.selectAtoms(" atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s N1 "%(seg,i)," atom %s %s C2  "%(seg,i))
	except:
		c = universe.selectAtoms(" atom %s %s O4\' "%(seg,i)," atom %s %s C1\' "%(seg,i)," atom %s %s N9 "%(seg,i)," atom %s %s C4  "%(seg,i))
	chi = c.dihedral()
	if chi < 360:
		chi = chi + 360
	return chi

def hydroxyl(seg,i):
        h = universe.selectAtoms(" atom %s %s C1\' "%(seg,i)," atom %s %s C2\' "%(seg,i)," atom %s %s O2\' "%(seg,i)," atom %s %s H2\'\' "%(seg,i))
        hydr = h.dihedral()
        if hydr < 360:
                hydr = hydr + 360
        return hydr  
