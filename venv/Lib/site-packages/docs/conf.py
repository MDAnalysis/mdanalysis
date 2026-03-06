# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------


# -- Project information -----------------------------------------------------

project = "MDAnalysis Sphinx theme"
html_title = "MDAnalysis Sphinx theme"
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("."))
import mdanalysis_sphinx_theme  # noqa: E402

authors = ", ".join(mdanalysis_sphinx_theme.__authors__)
copyright = "2023, " + authors

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "mdanalysis_sphinx_theme",
    "sphinxcontrib.autodoc_pydantic",
    "sphinx_search.extension",
]

# Autodoc settings
autosummary_generate = True

autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "member-order": "bysource",
}
autodoc_preserve_defaults = True
autodoc_typehints_format = "short"
# Workaround for autodoc_typehints_format not working for attributes
# see https://github.com/sphinx-doc/sphinx/issues/10290#issuecomment-1079740009
python_use_unqualified_type_names = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_attr_annotations = True
napoleon_custom_sections = [("attributes", "params_style")]
napoleon_use_rtype = False
napoleon_use_param = True
napoleon_use_ivar = True
napoleon_preprocess_types = True

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "mdanalysis_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named 'default.css' will overwrite the builtin 'default.css'.
html_static_path = ["_static"]

# If HTML theme settings isn't lines 90-160, remember to change customization.rst

# -- HTML theme settings ------------------------------------------------
html_show_sourcelink = True
html_sidebars = {
    "**": ["globaltoc.html", "localtoc.html", "searchbox.html", "navigation.html"],
    "customization": ["globaltoc.html", "searchbox.html"],
    "subpage/second-subsubpage": ["globaltoc.html", "searchbox.html"],
}

html_theme = "mdanalysis_sphinx_theme"

extra_nav_links = {}
extra_nav_links["MDAnalysis"] = "http://mdanalysis.org"
extra_nav_links["MDAnalysis docs"] = "http://docs.mdanalysis.org"
extra_nav_links["@mdanalysis"] = "https://twitter.com/mdanalysis"


# material theme options (see theme.conf for more information)
html_theme_options = {
    # whether to apply official MDAnalysis styling
    # e.g. using the official MDAnalysis logo and favicon
    # and using the MDAnalysis privacy policy
    "mda_official": True,
    # Extra navigation links to show on the sidebar, before the table of contents
    "extra_nav_links": extra_nav_links,
    # The background colour of the logo area in the navigation bar
    "sidebar_logo_background": "#ffffff",
    # The background colour of the top navigation bar on mobile
    "mobile_navbar_background": "dark-gray",

    # other options inherited from sphinx_rtd_theme
    # Only display logo and not the project name on sidebar
    "logo_only": True,
    # Display the version number on the sidebar
    "display_version": True,
    # Where to display "next" and "previous" buttons
    "prev_next_buttons_location": "bottom",
    # Add an icon next to external links
    "style_external_links": False,
    # If enabled, navigation entries are not expandable
    "collapse_navigation": True,
    # If enabled, the navigation bar scrolls with the main page
    "sticky_navigation": True,
    # Maximum depth of the contents tree -- set to -1 for unlimited
    "navigation_depth": 4,
    # Whether to include hidden toctrees in the navigation bar
    "includehidden": True,
    # Whether to hide page subheadings from navigation
    "titles_only": False,
}

# If HTML theme settings isn't lines 90-160, remember to change customization.rst

language = "en"
html_last_updated_fmt = ""

todo_include_todos = True

html_use_index = True
html_domain_indices = True

nbsphinx_execute = "always"
nbsphinx_kernel_name = "python3"

extlinks = {
    "duref": (
        "http://docutils.sourceforge.net/docs/ref/rst/" "restructuredtext.html#%s",
        "%s",
    ),
    "durole": ("http://docutils.sourceforge.net/docs/ref/rst/" "roles.html#%s", "%s"),
    "dudir": (
        "http://docutils.sourceforge.net/docs/ref/rst/" "directives.html#%s",
        "%s",
    ),
}
