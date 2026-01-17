# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import mock
import os
import sys
import matplotlib
import importlib.util
import sphinx_gallery
# NOTE: avoid storing unpicklable objects (like classes/functions) in Sphinx
# config values, otherwise Sphinx emits config.cache warnings when pickling the
# environment.

sys.path.insert(
    0, os.path.abspath("../../simcoon-python-builder/include/simcoon/python_wrappers")
)

import subprocess

# -- Run Doxygen to generate XML for Breathe ---------------------------------
# Build the Doxygen XML documentation if it doesn't exist.
# Note: When using the Makefile, Doxygen is run as a dependency before Sphinx,
# so this block won't execute. This fallback is useful for:
# - ReadTheDocs builds
# - Direct sphinx-build invocations without make
doxygen_xml_dir = os.path.join(os.path.dirname(__file__), "_build", "doxygen", "xml")
if not os.path.exists(doxygen_xml_dir):
    print("Running Doxygen to generate XML for Breathe...")
    subprocess.call("doxygen Doxyfile", shell=True, cwd=os.path.dirname(__file__))

sphinx_gallery.EXAMPLES_DIR = os.path.abspath("../examples")
print("Sphinx-Gallery examples dir:", sphinx_gallery.EXAMPLES_DIR)

project = "Simcoon"
copyright = "2025, 3MAH development team"
author = "3MAH development team"

# Respect the release of simcoon
# Default values
release = version = "unknown"

# Path to the version file
version_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../python-setup/simcoon/__version__.py")
)

# Load the __version__ variable safely if the file exists
if os.path.exists(version_path):
    spec = importlib.util.spec_from_file_location("simcoon_version", version_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    release = getattr(mod, "__version__", "unknown")
    version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_gallery.gen_gallery",
    "matplotlib.sphinxext.plot_directive",
    "breathe",
]

# -- Breathe Configuration ---------------------------------------------------
# Path to the Doxygen XML output
breathe_projects = {
    "simcoon": os.path.join(os.path.dirname(__file__), "_build", "doxygen", "xml")
}
breathe_default_project = "simcoon"
breathe_default_members = ("members", "undoc-members")

matplotlib.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)

# include matplotlib plots
plot_include_source = True
plot_formats = [("png", 100), "pdf"]
plot_html_show_source_link = True

# Ensure index.rst is the master file instead of 'contents.rst'
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
# html_theme = 'sphinxdoc'
html_logo = "_static/simcoonLogos_ss_fond.png"


html_theme_options = {
    "canonical_url": "https://microgen.readthedocs.io/en/latest/",
    # 'logo_only': False,
    # 'display_version': True,
    # 'prev_next_buttons_location': 'bottom',
    "style_external_links": True,
    # 'vcs_pageview_mode': '',
    "style_nav_header_background": "#24445C",
    # Toc options
    "collapse_navigation": False,
    # 'sticky_navigation': True,
    # 'navigation_depth': 4,
    # 'includehidden': True,
    # 'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Sphinx Gallery Options
sphinx_gallery_conf = {
    # path to your examples scripts
    "examples_dirs": ["../examples/"],
    # path where to save gallery generated examples
    "gallery_dirs": ["examples"],
    # Pattern to search for example files
    "filename_pattern": r"\.py",
    # Sort gallery example by file name instead of number of lines (default)
    # Use a dotted-path string instead of the imported class to keep this
    # configuration picklable.
    "within_subsection_order": "sphinx_gallery.sorting.FileNameSortKey",
    "ignore_pattern": r"/(results|data)$",
    "copyfile_regex": r"^(data/.+\.(txt|csv|json)|results/.+\.(txt|csv))$",
    "nested_sections": False,  # disables the extra toctree line
}
