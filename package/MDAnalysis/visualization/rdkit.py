# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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

try:
    get_ipython()
except NameError:
    raise ImportError("You must be in an interactive python shell (IPython or "
                      "Notebook) to use this")
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw.rdMolDraw2D import (
        PrepareMolForDrawing, MolDrawOptions, MolDraw2DCairo, MolDraw2DSVG)
except ImportError:
    raise ImportError("RDKit is needed to use the RDKit drawing code, but it "
                      "doesn't appear to be installed")

from functools import partial
from io import BytesIO
from base64 import b64encode

from IPython.display import display, HTML, SVG
from PIL import Image

from ..core.groups import AtomGroup


class RDKitDrawer:
    shell = get_ipython()
    _repr_format = {}

    def __init__(self, size=(450, 250), max_atoms=200, removeHs=True,
                 kekulize=True, drawOptions=None, useSVG=False):
        self.size = size
        self.max_atoms = max_atoms
        self.removeHs = removeHs
        self.kekulize = kekulize
        self.drawOptions = drawOptions or MolDrawOptions()
        self.reset_all_repr()
        self.useSVG = useSVG
        if useSVG:
            self.add_repr(AtomGroup, 'image/svg+xml', self._repr_atomgroup)
        else:
            self.add_repr(AtomGroup, 'image/png', self._repr_atomgroup)

    def _repr_atomgroup(self, ag):
        if ag.n_atoms <= self.max_atoms:
            return self.atomgroup_to_image(ag)

    def _prepare_atomgroup_for_drawing(self, ag, traj=False):
        mol = ag.convert_to("RDKIT")
        if self.removeHs:
            mol = Chem.RemoveHs(mol)
        # remove 3D coordinates for a clearer image
        if not traj:
            mol.RemoveAllConformers()
        try:
            mol = PrepareMolForDrawing(mol, kekulize=self.kekulize)
        except ValueError:
            mol = PrepareMolForDrawing(mol, kekulize=False)
        return mol

    def atomgroup_to_image(self, ag, **kwargs):
        traj = kwargs.pop("traj", False)
        mol = self._prepare_atomgroup_for_drawing(ag, traj)
        size = kwargs.pop("size", self.size)
        useSVG = kwargs.pop("useSVG", self.useSVG)
        drawer = MolDraw2DSVG if useSVG else MolDraw2DCairo
        d2d = drawer(*size)
        opts = kwargs.pop("drawOptions", self.drawOptions)
        opts.prepareMolsBeforeDrawing = False
        d2d.SetDrawOptions(opts)
        d2d.DrawMolecule(mol, **kwargs)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def atomgroup_to_gif(self, ag, output=None, legend="Frame {}",
                         frame_duration=200, traj_slice=slice(0, None, None),
                         **kwargs):
        mol = self._prepare_atomgroup_for_drawing(ag, traj=True)
        if mol and hasattr(ag.universe, "trajectory"):
            pngs = []
            for ts in ag.universe.trajectory[traj_slice]:
                img = self.atomgroup_to_image(ag, traj=True, useSVG=False,
                                              legend=legend.format(ts.frame),
                                              **kwargs)
                pngs.append(img)
            img, *imgs = [Image.open(BytesIO(png)) for png in pngs]
            # write to file, or display if output is None
            buffer = partial(open, output, "wb") if output else BytesIO
            with buffer() as fp:
                img.save(fp=fp, format='GIF', append_images=imgs,
                        save_all=True, duration=frame_duration, loop=0)
                if isinstance(fp, BytesIO):
                    b64 = b64encode(fp.getvalue()).decode("ascii")
                    display(HTML(f"<img src='data:image/gif;base64,{b64}' />"))

    def reset_repr(self, obj=AtomGroup, fmt="image/png"):
        self.shell.display_formatter.formatters[fmt].pop(obj)
        self._repr_format.pop(obj)

    def reset_all_repr(self):
        for obj, fmt in self._repr_format.items():
            self.shell.display_formatter.formatters[fmt].pop(obj)
        self._repr_format.clear()

    def add_repr(self, obj, fmt, func):
        # adds a new formatter to obj
        formatters = self.shell.display_formatter.formatters
        formatters[fmt].for_type(obj, func)
        self._repr_format[obj] = fmt


# add default AtomGroup rich display on import
RDKitDrawer()
