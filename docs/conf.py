# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'coremstools'
copyright = '2024, Christian Dewey & Rene Boiteau'
author = 'Christian Dewey & Rene Boiteau'
release = '0.0.1'

import os
import sys
sys.path.insert(0,os.path.abspath('..'))
sys.path.append(os.path.abspath('/Users/christiandewey/CoreMS'))
sys.path.append(os.path.abspath('/Users/christiandewey/CoreMS/ext_lib'))
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc','sphinx.ext.coverage','sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'Python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

