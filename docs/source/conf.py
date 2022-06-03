# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------


project = 'MD DaVis'
copyright = '2019, Dibyajyoti Maity'
author = 'Dibyajyoti Maity'

# The short X.Y version.
version = release = '0.4.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'myst_parser',
    'sphinx_panels',
    'sphinx_click',
    'sphinx_tabs.tabs',
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary

pygments_style = 'trac'
extlinks = {
    'issue': ('https://github.com/djmaity/md_davis/issues/%s', '#'),
    'pr': ('https://github.com/djmaity/md_davis/pull/%s', 'PR #'),
}
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_extra_path = ['_extras']
html_use_smartypants = True
html_last_updated_fmt = '%b %d, %Y'
html_split_index = False
html_sidebars = {
   '**': ['search-field.html', 'globaltoc.html', 'sidebar-ethical-ads.html']
}
html_short_title = '%s-%s' % (project, version)
html_css_files = [
    'css/custom.css',
]
html_logo = "_static/MD_DaVis_Logo.png"
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/djmaity/md-davis",
            "icon": "fab fa-github-square",
        },
    ],
    "favicons": [
        {
            "rel": "icon",
            "sizes": "16x16",
            "href": "https://secure.example.com/favicon/favicon-16x16.png",
        },
        {
            "rel": "icon",
            "sizes": "32x32",
            "href": "favicon-32x32.png",
        },
        {
            "rel": "apple-touch-icon",
            "sizes": "180x180",
            "href": "apple-touch-icon-180x180.png"
        },
    ],
    "navbar_end": ["navbar-icon-links.html"]
}
html_context = {
    # "github_url": "https://github.com", # or your GitHub Enterprise interprise
    "github_user": "djmaity",
    "github_repo": "md-davis",
    "github_version": "master",
    "doc_path": "docs/source",
}

autosectionlabel_prefix_document = True

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False

apidoc_module_dir = '../../md_davis'
apidoc_output_dir = 'reference'
apidoc_excluded_paths = ['tests']
apidoc_separate_modules = True

todo_include_todos = True

coverage_show_missing_items = True
