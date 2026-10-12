# Sphinx configuration for the CPMD manual.
#
# Build the web pages with   sphinx-build -b html doc doc/_build/html
# and the PDF manual with    sphinx-build -M latexpdf doc doc/_build

import os
import sys
from datetime import date

sys.path.insert(0, os.path.abspath('_ext'))

with open(os.path.join(os.path.dirname(__file__), '..', 'VERSION')) as f:
    release = f.read().strip()
version = release

project = 'CPMD'
author = 'The CPMD consortium'
copyright = '1990-%d IBM Corp., 1994-2001 Max Planck Institute Stuttgart' % date.today().year

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'cpmd_keywords',
]

primary_domain = 'cpmd'
default_role = 'literal'
todo_include_todos = True
numfig = True
math_numfig = True
math_eqref_format = '{number}'
exclude_patterns = ['_build', 'README.md', 'cmake.md']
linkcheck_ignore = [r'mailto:.*']
linkcheck_timeout = 20

# -- HTML ------------------------------------------------------------------
html_theme = 'furo'
html_title = 'CPMD %s manual' % release
html_static_path = ['_static']
html_css_files = ['cpmd.css']
html_theme_options = {
    'source_repository': 'https://github.com/ivan-pi/CPMD/',
    'source_branch': 'main',
    'source_directory': 'doc/',
}

# -- LaTeX / PDF -------------------------------------------------------------
latex_engine = 'pdflatex'
latex_documents = [
    ('index', 'cpmd-manual.tex', 'CPMD %s Manual' % release, author, 'manual'),
]
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'preamble': r'''
\usepackage{amsmath,amssymb}
\setcounter{tocdepth}{1}
''',
    'printindex': r'\footnotesize\raggedright\printindex',
}
latex_show_urls = 'footnote'
