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

"""RDKit Drawer --- :mod:`MDAnalysis.visualization.rdkit`
=========================================================

A module to wrap RDKit's drawing code and modify the representation of small
AtomGroups in interactive notebooks.
"""

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

from .base import FormatterBase
from ..core.groups import AtomGroup
from .. import _FORMATTERS


class RDKitDrawer(FormatterBase):
    """Wrapper class to RDKit's drawing code

    .. versionadded:: 2.0.0
    """
    format = "RDKIT"

    def __init__(self, size=(450, 250), max_atoms=200, removeHs=True,
                 kekulize=True, drawOptions=None, useSVG=False):
        """Set default parameters for the drawer

        Parameters
        ----------
        size : tuple
            default size for images
        max_atoms : int
            AtomGroups with more atoms that this limit won't be displayed as
            images and will use the default representation instead
        removeHs : bool
            Remove hydrogens from the image
        kekulize : bool
            Use the Kekule representation for aromatic systems
        drawOptions : None or rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
            RDKit MolDrawOption object passed when calling the drawing code.
            Use it to set the label size, highlight color...etc.
        useSVG : bool
            Use SVG images instead of PNG
        """
        self.size = size
        self.max_atoms = max_atoms
        self.removeHs = removeHs
        self.kekulize = kekulize
        self.drawOptions = drawOptions or MolDrawOptions()
        # remove any previous RDKIT formatters
        self.reset_all_repr()
        self.useSVG = useSVG
        # add the new formatter
        if useSVG:
            self.add_repr(AtomGroup, 'image/svg+xml', self._repr_atomgroup)
        else:
            self.add_repr(AtomGroup, 'image/png', self._repr_atomgroup)

    def _repr_atomgroup(self, ag):
        """Wrapper method to :meth:`~atomgroup_to_image`

        Returns
        -------
        repr : str or None
            Raw image data (str) if the AtomGroup has less atoms than
            ``max_atoms`` or else (None) the default AtomGroup representation
        """
        if ag.n_atoms <= self.max_atoms:
            return self.atomgroup_to_image(ag)

    def _prepare_atomgroup_for_drawing(self, ag, keep_3D=False):
        """Prepare the AtomGroup and resulting mol for the drawing code

        Parameters
        ----------
        ag : MDAnalysis.core.groups.AtomGroup
            The AtomGroup to prepare
        keep_3D : bool
            Keep or remove 3D coordinates to generate the image
        """
        mol = ag.convert_to("RDKIT")
        if self.removeHs:
            mol = Chem.RemoveHs(mol)
        # remove 3D coordinates for a clearer image
        if not keep_3D:
            mol.RemoveAllConformers()
        try:
            mol = PrepareMolForDrawing(mol, kekulize=self.kekulize)
        except ValueError:
            mol = PrepareMolForDrawing(mol, kekulize=False)
        return mol

    def atomgroup_to_image(self, ag, **kwargs):
        r"""Create an image from an AtomGroup

        Parameters
        ----------
        ag : MDAnalysis.core.groups.AtomGroup
            The AtomGroup to display
        keep_3D : bool
            Keep or remove 3D coordinates to generate the image
        size : tuple
            size of the output image
        useSVG : bool
            Output an SVG instead of a PNG
        drawOptions : None or rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
            Parameters passed to RDKit for drawing the molecule
        legend : str
            Legend displayed on the image
        **kwargs : object
            Other parameters passed to :meth:`~rdkit.Chem.Draw.rdMolDraw2D.MolDraw2D.DrawMolecule`

        Returns
        -------
        img : str
            Raw PNG or SVG code
        """
        keep_3D = kwargs.pop("keep_3D", False)
        mol = self._prepare_atomgroup_for_drawing(ag, keep_3D)
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
                         frame_duration=200, start=None, stop=None, step=None,
                         **kwargs):
        r"""Create a GIF from an AtomGroup

        Parameters
        ----------
        ag : MDAnalysis.core.groups.AtomGroup
            The AtomGroup to display
        output : None or str
            Either a path to save the gif (str), or ``None`` to display the GIF
            inline
        legend : str
            Format string used for the legend of the GIF. ``{}`` will be
            replaced by the frame number
        frame_duration : int or list
            Duration of each frame for the GIF
        start : None or int
            Start frame for the GIF
        stop : None or int
            End frame for the GIF
        step : None or int
            Skip 1 frame for every ``step`` frames in the trajectory
        **kwargs : object
            Other parameters used by :meth:`~atomgroup_to_image`
        """
        mol = self._prepare_atomgroup_for_drawing(ag, keep_3D=True)
        if mol and hasattr(ag.universe, "trajectory"):
            pngs = []
            for ts in ag.universe.trajectory[start:stop:step]:
                img = self.atomgroup_to_image(ag, keep_3D=True, useSVG=False,
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


# add default AtomGroup rich display on import
RDKitDrawer()
