"""Sphinx configuration for the SPECTRA documentation."""

import os
import sys
from datetime import date


sys.path.insert(0, os.path.abspath("../src"))

project = "spectra-cobra"
author = "Pavan Kumar S"
copyright = f"{date.today().year}, Nirav Lab"

try:
    from spectra_cobra import __version__ as release
except ImportError:  # pragma: no cover - docs build without the package
    release = "0.1.0.dev0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.graphviz",
]

# Graphviz: SVG keeps the network diagrams sharp at any zoom.
graphviz_output_format = "svg"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = f"{project} {release}"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
}

# Numpydoc-style docstrings, which is what the package uses.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_member_order = "bysource"
autodoc_typehints = "description"
# Members are declared per directive, not globally: the flat API index needs
# to suppress them so the dedicated pages own the anchors.
autodoc_default_options = {
    "show-inheritance": True,
}
# cobra is heavy and not needed to render signatures on a docs builder.
autodoc_mock_imports = []

autosummary_generate = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "cobra": ("https://cobrapy.readthedocs.io/en/latest/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

nitpicky = False
