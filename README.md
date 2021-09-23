MD DaVis
========

[![PyPI Package latest release][pypi-badge]][pypi-link]
[![Documentation Status][docs-badge]][docs-link]
[![PyPI Wheel][wheel-badge]][wheel-link]
[![Supported versions][python-badge]][python-link]
[![Commits since latest release][commits-badge]][commits-link]

A tool for comparative analysis of molecular dynamics simulations of proteins.

* Free software: [MIT license](https://github.com/djmaity/md-davis/blob/master/LICENSE)
* Documentation: [https://md-davis.readthedocs.io](https://md-davis.readthedocs.io)

Introduction
------------
MD DaVis is a tool to facilitate the comparative analysis of molecular dynamics simulations of proteins.

Features
--------
1. Free energy landscape
2. Residue properties plot
3. Surface electrostatics
4. Electric field dynamics
5. H-bond/Contact matrix

Installation
------------
The easiest installation method is with
[Anaconda](https://www.anaconda.com/products/individual) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html).
It is highly recommended to install MD DaVis in a virtual environment.
The [environment.yml](https://github.com/djmaity/md-davis/blob/master/environment.yml)
is provided to ease the process.
```
conda env create -f environment.yml -n md_davis_env
```
This automatically creates a conda environment called `md_davis_env` with all required dependencies.
Activate the environment and install MD DaVis in it using:
```
conda activate md_davis_env
pip install md-davis
```

For detailed installation instructions, see:
[https://md-davis.readthedocs.io/en/latest/install.html](https://md-davis.readthedocs.io/en/latest/install.html)

### Dependencies
* [PyMOL](https://pymol.org/2/)
* [Delphi](http://compbio.clemson.edu/delphi)
* [MSMS](http://mgltools.scripps.edu/downloads#msms)

Usage
-----

If MD DaVis is installed in a virtual environment, it should be activated before using MD DaVis.
For example, to activate the virtual environment created with conda, use:
```
conda activate md_davis_env
```

The MD DaVis CLI can be called with:
```
md_davis
```

Author
------
Dibyajyoti Maity - [www.djmaity.com](http://www.djmaity.com)

<!--  --------------------------------------------------------------------- -->
[pypi-badge]: https://img.shields.io/pypi/v/md-davis.svg
[pypi-link]: https://pypi.org/project/md-davis
[docs-badge]: https://readthedocs.org/projects/md-davis/badge/?style=flat
[docs-link]: https://md-davis.readthedocs.io/
[wheel-badge]: https://img.shields.io/pypi/wheel/md-davis.svg
[wheel-link]: https://pypi.org/project/md-davis
[python-badge]: https://img.shields.io/pypi/pyversions/md-davis.svg
[python-link]: https://pypi.org/project/md-davis
[commits-badge]: https://img.shields.io/github/last-commit/djmaity/md-davis
[commits-link]: https://github.com/djmaity/md-davis/
