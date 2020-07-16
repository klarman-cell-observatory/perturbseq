# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

#settings inspired from scanpy's documentation
import os
import sys
import warnings
from pathlib import Path
from datetime import datetime


HERE=Path(__file__).parent
sys.path.insert(0, str(HERE.parent))
on_rtd=os.environ.get('READTHEDOCS') == True
print(sys.path)

# -- Project information -----------------------------------------------------

project = 'perturbseq'
copyright = '2020, Oana Ursu'
author = 'Oana Ursu'

# The full version, including alpha/beta/rc tags
release = '0.0.0'
needs_sphinx = '2.0' 

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = ['sphinx.ext.autodoc','sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
#'sphinx_autodoc_typehints',
"sphinx_rtd_theme"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
source_suffix='.rst'
master_doc='index'
default_role='literal'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style='sphinx'

# -- Options for HTML output -------------------------------------------------


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

autosummary_generate = True
autodoc_member_order = 'bysource'
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False
api_dir = HERE / 'api'  # function_images

html_theme = 'sphinx_rtd_theme'
#html_theme_options = dict(navigation_depth=4, logo_only=True)  # Only show the logo
html_context = dict(
    display_github=True,  
    github_user='oursu',  
    github_repo='perturbseq',  
    github_version='master',  
    conf_py_path='/docs/',  # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_show_sphinx = False
