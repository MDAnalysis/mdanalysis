"""MDAnalysis Sphinx theme."""

import logging
import os
import pathlib

import sass
from sass import SassColor, SassFunction
from sphinx.util import console

from ._version import get_versions

logger = logging.getLogger(__name__)

try:
    from .authors import __authors__
except ImportError:
    logger.info('Could not find authors.py, __authors__ will be empty.')
    __authors__ = []

__version__ = get_versions()["version"]
del get_versions


def setup(app):
    """Setup connects events to the sitemap builder"""
    app.connect("build-finished", compile_css)

    app.site_pages = []
    app.add_html_theme(
        "mdanalysis_sphinx_theme", html_theme_path()[0]
    )
    return {
        "version": __version__,
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }


def hex_to_rgb(hex):
    """Convert a hex color to RGB"""
    hex = hex.lstrip("#")
    return tuple(int(hex[i : i + 2], 16) for i in (0, 2, 4))


def compile_css(app, exception):
    """Compile Bulma SASS into CSS"""
    if exception is not None:
        return

    theme_path = pathlib.Path(html_theme_path()[0])
    src = theme_path / "sass/site.sass"
    dest = pathlib.Path(app.outdir) / "_static/site.css"

    if not dest.parent.exists():
        return

    config = app.config["html_theme_options"]
    COLORS = {
        "mdanalysis-orange": (255, 146, 0),
        "mdanalysis-code-orange": (202, 101, 0),
        "white": (255, 255, 255),
        "dark-gray": (52, 49, 49),
    }
    theme_defaults = {
        "color_accent": "mdanalysis-code-orange",
        "sidebar_logo_background": "white",
        "mobile_navbar_background": "dark-gray",
    }
    function_colors = {}
    custom_sass_functions = {}
    for option, default in theme_defaults.items():
        theme_option = config.get(option, default)
        if theme_option in COLORS:
            color = COLORS[theme_option]
        else:
            color = hex_to_rgb(theme_option)
        function_colors[option] = color

    if config.get("css_minify", False):
        output_style = "compressed"
        source_comments = False
    else:
        output_style = "expanded"
        source_comments = True

    custom_sass_functions["mobile_navbar_background"] = lambda: SassColor(*function_colors["mobile_navbar_background"], 1)
    custom_sass_functions["sidebar_logo_background"] = lambda: SassColor(*function_colors["sidebar_logo_background"], 1)
    custom_sass_functions["color_accent"] = lambda: SassColor(*function_colors["color_accent"], 1)
    custom_sass_functions["hyphenate"] = lambda: config.get(
        "html_hyphenate_and_justify", False
    )

    css = sass.compile(
        filename=str(src),
        output_style=output_style,
        custom_functions=custom_sass_functions,
    )

    print(f"Writing compiled SASS to {console.colorize('blue', str(dest))}")

    with open(dest, "w") as f:
        print(css, file=f)


def html_theme_path():
    return [os.path.dirname(os.path.abspath(__file__))]


from . import _version
__version__ = _version.get_versions()['version']
