# -*- coding: utf-8 -*-
import os

needs_sphinx = '1.1'
extensions = ['sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
source_encoding = 'utf-8'
master_doc = 'index'
project = u'Debyer'
version = 'rev86'
release = version
language = 'en'
exclude_patterns = ['_build']
default_role = 'math'
pygments_style = 'sphinx'


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'default'
html_title = "Debyer documentation"
html_short_title = "Debyer docs"
html_logo = 'debye-logo.png'
#html_favicon = None
html_static_path = ['_static']
#html_last_updated_fmt = '%b %d, %Y' # ''
#html_use_smartypants = True
html_sidebars = { '**': ['localtoc.html', 'relations.html'] }
#html_additional_pages = {}
html_domain_indices = False
html_use_index = False
html_show_sourcelink = True
#html_show_sphinx = True
html_show_copyright = False
#html_file_suffix = None
html_static_path = ['debyer.css']
html_style = 'debyer.css'

# Apr 2015: the default value on RTD is http://... and doesn't work with https
if os.getenv('READTHEDOCS'):
    mathjax_path = '//cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'

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


