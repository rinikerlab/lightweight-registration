import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'lightweight-registration'
copyright = '2025, Jessica Braun and Greg Landrum'
author = 'Jessica Braun and Greg Landrum'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

import os
import sys
print("KOKO:", os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('..'))
print("KOKO:", os.path.abspath('../lwreg'))
sys.path.insert(0, os.path.abspath('../lwreg'))