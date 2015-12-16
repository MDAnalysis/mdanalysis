#!/usr/bin/env python
# Example script, part of MDAnalysis
"""
:Author: Jan Domanski
:Year: 2010
:Copyright: GNU Public License v3

MDAnalysis example: Leaflet indentification
===========================================
The script was used with MARTINI simulations to monitor the multimeric state
of several (12) WALP/KALP peptides of variable lenght.
Aggregative properties of an array of identical peptides can be traced with
time and plotted using matlibplot module in python.

The script ASSUMES that the pepties are defined on the beginning of conf file.

In house, the script will also be extended to calculate eg. angle between
peptides in dimeric arrangement.

Example of figure obtained, Fig 3 in: (there, a different script is used
but to the same end)
L.V. Schafer, D.H. de Jong, A. Holt, A.J. Rzepiela, A.H. de Vries, B. Poolman,
J.A. Killian, S.J. Marrink. Lipid packing drives the segregation of
transmembrane helices  into disordered lipid domains in model biomembranes.
PNAS, doi:10.1073/pnas.1009362108, open access
"""

import MDAnalysis
from pylab import *

# periodicity has to be turned on in the search, turn off the KDTree that is
# faster but ignores periodicity
MDAnalysis.core.flags['use_periodic_selections'] = True
MDAnalysis.core.flags['use_KDTree_routines'] = False

conf = "conf.gro"
traj = "traj.xtc"
universe = None

# peptide configuration
peptide_selection, peptide_dictionary = ({}, {})
peptide_conf = {"number": 12, "lenght": 23}

# geometry search cutoff (in Angstrom)
cutoff = 8


def __main__():
    universe = MDAnalysis.Universe(conf, traj)

    peptide_selection, peptide_dictionary = define_peptides(peptide_conf["number"], peptide_conf["lenght"])

    partners, clusters, multimers = analyze(skip=1000)

    plot(multimers)


def define_peptides(number_of_peptides, lenght_of_peptide):
    """
    Returns two dictionaries,
    selection - key is peptide id, values is selection string
    e.g. {0: "resid 1-31", 1: "resid 32-62"}

    lookup - key is residue id, value is peptide id. it is used when a protein
    bead is found but the peptide to which it belongs needs to be identified
    e.g. {1: 0, 2: 0, 3:0 ... 31:0, 32: 1, 33: 1}
    """
    selection = {}
    lookup = {}

    for i in range(number_of_peptides):
        selection[i] = "resid {0:d}-{1:d}".format(lenght_of_peptide * i + 1, lenght_of_peptide * (i + 1))
        index = lenght_of_peptide * i + 1
        while (index <= lenght_of_peptide * (i + 1)):
            lookup[index] = i
            index += 1

    return selection, lookup


def analyze(partners=None, clusters=None, multimers=None, skip=1000):
    # initialize multimers var, this var stores output of the analysis in a
    # format that is easy to plot. at least in the deafault implementation.
    if partners is None:
        partners = {}
    if clusters is None:
        clusters = {}
    if multimers is None:
        multimers = [{}, {}, {}, {}]
    for d in multimers:
        d["x"] = []
        d["y"] = []

    for ts in universe.trajectory:
        if not ts.frame % skip == 0 and ts.frame != 1:
            continue
        print "Stepping... Frame {0:d}, time {1:d} ns".format(ts.frame, ts.time / 1000)

        p = find_partners(peptide_selection, peptide_dictionary)
        partners[ts.frame] = p

        # the c(lusters) variable is probably the best starting point for
        # implementing your own processing of the data - contact analysis,
        # helix tilt in dimers present - whatever. see the function doc.
        c = find_clusters(p)
        clusters[ts.frame] = c

        # this gets messy, since the data has to be rearranged spec. to the
        # plotting library used
        t0, t1, t2, t3 = find_multimers(c)

        multimers[0]["x"].append(ts.time / 1000)
        multimers[0]["y"].append(t0)
        multimers[1]["x"].append(ts.time / 1000)
        multimers[1]["y"].append(t1)
        multimers[2]["x"].append(ts.time / 1000)
        multimers[2]["y"].append(t2)
        multimers[3]["x"].append(ts.time / 1000)
        multimers[3]["y"].append(t3)

    return partners, clusters, multimers


def find_partners(peptide_list, lookup):
    """
    Perform the cut-off restricted search in the proximity of all peptides
    present in the system and identitfies, for a given peptide, the set
    (see python doc on 'set' data structure) of interacting peptides.
    e.g. {0: [(3, 5)], 3: [(0,5)], 5: [(3,0)]}
    """
    ret = {}
    for id, selection in peptide_list.items():
        atom_list = universe.select_atoms("around {0:d} ({1!s}) and not resname W".format(cutoff, selection))
        ret[id] = set()
        for atom in atom_list:
            if atom.resname == "CHOL" or atom.resname == "DPPC" or atom.resname == "DUPC":
                continue
            if not atom.resid in lookup:
                continue
            walp_id = lookup[atom.resid]
            if walp_id in ret[id]:
                continue
            # print atom
            ret[id].add(walp_id)
    return ret


def find_clusters(partners_dictionary):
    """
    Function recieves as input the data generated by 'find_partners' in the
    format :
    {0: [(3, 5)], 3: [(0,5)], 5: [(3,0)]}
    Key is peptide id, the values is a set of its partners.
    As one can see, the contacts defined here are redundant and not in a useful
    form.
    What would be best is output in the form:
    [set(0,3,5), set(1), set(2), set(4)...]
    assuming that only peptides 0,3,5 form a trimer and rest is present as
    monomeres.

    TODO there is a way to elegantize this function, it's easy. i leave this as
    a riddle to the users.
    """
    # stores a list of set objects
    cluster_list = []

    cluster_dict = {}

    # history of assigned peptides
    history = set()

    for id1, partner_list in partners_dictionary.items():
        history.add(id1)
        if len(partner_list) == 0:
            continue
        for id2 in partner_list:
            if not id1 in cluster_dict and not id2 in cluster_dict:
                s = set()
                s.add(id1)
                s.add(id2)
                cluster_list.append(s)
                cluster_dict[id1] = s
                cluster_dict[id2] = s
                if id2 in history:
                    history.remove(id2)
                if id1 in history:
                    history.remove(id1)
            elif id1 in cluster_dict and not id2 in cluster_dict:
                s = cluster_dict[id1]
                s.add(id2)
                cluster_dict[id2] = s
                if id2 in history:
                    history.remove(id2)
            elif id2 in cluster_dict and not id1 in cluster_dict:
                s = cluster_dict[id2]
                s.add(id1)
                cluster_dict[id1] = s
                if id1 in history:
                    history.remove(id1)
    # add the monomers
    for peptide in history:
        if peptide in cluster_dict:
            continue
        monomer = set()
        monomer.add(peptide)
        cluster_list.append(monomer)

    return cluster_list


def find_multimers(clusters_list):
    multimers = [0, 0, 0, 0]  # mono, di, tr, multimers
    for cluster in clusters_list:
        i = len(cluster)
        if i == 1:
            multimers[0] += 1
        elif i == 2:
            multimers[1] += 1
        elif i == 3:
            multimers[2] += 1
        else:
            multimers[3] += 1
    return multimers[0], multimers[1], multimers[2], multimers[3]


def plot(multimers):
    """
    Default plotting function. Change according to your needs.
    """
    title("KALP23")

    subplot(411)
    step(multimers[0]["x"], multimers[0]["y"], color="black", linewidth=2)
    ylim(0, 14)
    xlim(0, 8000)
    text(6000, 10, 'monomeres', fontsize=14)

    subplot(412)
    step(multimers[1]["x"], multimers[1]["y"], color="black", linewidth=2)
    ylim(0, 14)
    xlim(0, 8000)
    text(6000, 10, 'dimers', fontsize=14)

    subplot(413)
    step(multimers[2]["x"], multimers[2]["y"], color="black", linewidth=2)
    ylim(0, 14)
    xlim(0, 8000)
    text(6000, 10, 'trimers', fontsize=14)

    subplot(414)
    step(multimers[3]["x"], multimers[3]["y"], color="black", linewidth=2)
    ylim(0, 14)
    xlim(0, 8000)
    text(6000, 10, 'more', fontsize=14)

    xlabel("Time (ns)")
    show()


if __name__ == "__main__":
    main()
