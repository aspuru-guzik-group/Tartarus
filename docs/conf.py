# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from unittest.mock import Mock

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Tartarus'
copyright = '2023, Akshat Nigam'
author = 'Akshat Nigam'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_mock_imports = ["numpy", "torch", "sascorer", "openbabel", "pyscf", "morfeus", "polanyi", "scipy", "sklearn"]
autosummary_generate = True

MOCK_MODULES = ["rdkit", "rdkit.Chem", "rdkit.Chem.Crippen", "rdkit.Chem.Descriptors", "rdkit.Chem.Lipinski", "rdkit.Chem.rdmolops", "rdkit.Chem.rdMolDescriptors"]
for module_name in MOCK_MODULES:
    if module_name == "rdkit.Chem":
        Chem = Mock()
        Chem.RDConfig = Mock()
        Chem.RDConfig.RDContribDir = '/tmp'
        sys.modules[module_name] = Chem
    else:
        sys.modules[module_name] = Mock()

os.environ["XTBHOME"] = "/tmp"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
