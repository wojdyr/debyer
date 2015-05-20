# -*- coding: utf-8 -*-
import os

needs_sphinx = '1.3'
extensions = ['sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
source_encoding = 'utf-8'
master_doc = 'index'
project = u'Debyer'
version = '0.4'
release = version
language = 'en'
exclude_patterns = ['_build']
default_role = 'math'
pygments_style = 'sphinx'


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'alabaster'
html_title = "Debyer documentation"
html_short_title = "Debyer docs"
html_logo = 'debye-logo.png'
#html_favicon = None
#html_last_updated_fmt = '%b %d, %Y' # ''
#html_use_smartypants = True
html_sidebars = { '**': ['localtoc.html'] }
#html_additional_pages = {}
html_domain_indices = False
html_use_index = False
html_show_sourcelink = True
#html_show_sphinx = True
copyright = '2015, Marcin Wojdyr'
html_show_copyright = False
#html_file_suffix = None
#html_static_path = ['_static']
#html_static_path = ['debyer.css']
#html_style = 'debyer.css'
html_theme_options = {
        'analytics_id': 'UA-17365358-2',
        }

# Apr 2015: the default value on RTD is http://... and doesn't work with https
#if os.getenv('READTHEDOCS'):
#    mathjax_path = 'https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

latex_elements = { 'papersize': 'a4paper', 'pointsize': '11pt',
                    #'preamble': '',
                 }

# (source start file, target name, title, author, documentclass [howto/manual])
latex_documents = [
  ('index', 'Debyer.tex', 'Debyer Documentation', '', 'howto'),
]

latex_logo = 'debye-logo.png'
#latex_use_parts = False
#latex_show_pagerefs = False
#latex_show_urls = False
latex_domain_indices = False

