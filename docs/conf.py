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
import datetime
import mock
import os
import sys
import matplotlib
import importlib.util
import sphinx_gallery

# Imported before sphinx-gallery captures stdout: an editable install rebuilds
# on first import and its cmake subprocess needs a real file descriptor.
import simcoon  # noqa: F401
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
    os.makedirs(os.path.dirname(doxygen_xml_dir), exist_ok=True)
    subprocess.call("doxygen Doxyfile", shell=True, cwd=os.path.dirname(__file__))

sphinx_gallery.EXAMPLES_DIR = os.path.abspath("../examples")

sys.path.insert(0, os.path.dirname(__file__))
import equation_thumbnails

equation_thumbnails.build(os.path.join(os.path.dirname(__file__), "_static", "thumbnails"))
print("Sphinx-Gallery examples dir:", sphinx_gallery.EXAMPLES_DIR)

project = "Simcoon"
author = "3MAH development team"
copyright = f"2025-{datetime.date.today().year}, {author}"

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
    "sphinx_design",
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

# Suppress specific warnings
suppress_warnings = ["toc.not_included"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
# html_theme = 'sphinxdoc'
html_logo = "_static/simcoonLogos_ss_fond.png"
html_title = f"simcoon {release}"


html_theme_options = {
    "canonical_url": "https://3mah.github.io/simcoon-docs/",
    "logo_only": True,
    "style_external_links": True,
    "style_nav_header_background": "#24445C",
    "collapse_navigation": False,
    "navigation_depth": 3,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sphinx = False

# -- Sphinx Gallery Options
sphinx_gallery_conf = {
    # path to your examples scripts
    "examples_dirs": ["../examples/"],
    # path where to save gallery generated examples
    "gallery_dirs": ["examples"],
    # Pattern to search for example files
    "filename_pattern": r"\.py",
    # Sections in reading order, from the material point to the composite.
    "subsection_order": [
        f"../examples/{category}"
        for category in (
            "mechanical",
            "thermomechanical",
            "hyperelasticity",
            "continuum_mechanics",
            "analysis",
            "heterogeneous",
            "identification",
            "ml",
            "solver",
            "umats",
            "umat_tutorial",
            "data",
        )
    ],
    # Sort gallery example by file name instead of number of lines (default)
    # Use a dotted-path string instead of the imported class to keep this
    # configuration picklable.
    "within_subsection_order": "sphinx_gallery.sorting.FileNameSortKey",
    "ignore_pattern": r"/(results|data)$",
    "copyfile_regex": r"^(data/.+\.(txt|csv|json)|results/.+\.(txt|csv))$",
    "nested_sections": False,  # disables the extra toctree line
}
