MD DaVis
========

[![PyPI Package latest release][pypi-badge]][pypi-link]
[![Travis CI][travis-badge]][travis-link]
[![Documentation Status][docs-badge]][docs-link]
[![PyPI Wheel][wheel-badge]][wheel-link]
[![Supported versions][python-badge]][python-link]
[![Commits since latest release][commits-badge]][commits-link]

A tool for comparative analysis of molecular dynamics simulations of proteins.

* Free software: MIT license
* Documentation: [https://md_davis.readthedocs.io](https://md_davis.readthedocs.io)

Introduction
------------

MD DaVis is a tool to facilitate the comparative analysis of molecular 
dynamics trajectories


Features
--------

1. Residue data plot
2. Free energy landscape
3. Surface electrostatics
4. Additional tools


Installation
------------

For detailed installation instructions see: 
[https://md_davis.readthedocs.io](https://md_davis.readthedocs.io)

I highly recommend installing this package in a virtual environment created with 
[Anaconda](https://www.anaconda.com/products/individual) or 
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) distribution. 
Make sure to activate the environment before running the install commands
```console
conda create --name md_davis_env python>=3.7
conda activate md_davis_env
```
`md_davis_env` is the name of the virtual environment and you can choose any name you like.
The environment must be activated with `conda activate md_davis_env` before using MD DaVis.

On Windows machines the install command fails because `pip` tries to compile 
the dependencies. Please install the following dependencies with conda before 
running the `pip` command:
```console
conda install -c conda-forge mdtraj
conda install -c conda-forge pymol-open-source
```

To install MD DaVis, run this command in your terminal:
```shell
pip install md-davis
```

Dependencies
------------

* [Open-Source PyMOL](https://github.com/schrodinger/pymol-open-source) available from `conda-forge` is for 64-bit linux and 
Windows systems only and requires Python > 3.7

Commercial version of pymol can be installed with:
```console
conda install -c schrodinger pymol-bundle
```
This can also be used with an Educational PyMOL [license](https://pymol.org/edu/?q=educational)

mdtraj is available for linux-64, osx-64, win-32, and win-64

Python dependencies are automatically installed. However, electrostatic 
calculation requires on following two programs which must be downloaded 
and installed separately.
* [Delphi](http://compbio.clemson.edu/delphi)
* [MSMS](http://mgltools.scripps.edu/downloads#msms)

To Do
-----

This page is still under construction. Visit `https://djmaity.github.io/md_davis/ <https://djmaity.github.io/md_davis/>`_ to know more.

Author
------

Dibyajyoti Maity - [www.djmaity.com](www.djmaity.com)

<!--  --------------------------------------------------------------------- -->
[pypi-badge]: https://img.shields.io/pypi/v/md-davis.svg
[pypi-link]: https://pypi.org/project/md-davis
[travis-badge]: https://travis-ci.org/uiri/toml.svg?branch=master
[travis-link]: https://travis-ci.org/
[docs-badge]: https://readthedocs.org/projects/md_davis/badge/?style=flat
[docs-link]: https://md_davis.readthedocs.io/
[wheel-badge]: https://img.shields.io/pypi/wheel/md-davis.svg
[wheel-link]: https://pypi.org/project/md-davis
[python-badge]: https://img.shields.io/pypi/pyversions/md-davis.svg
[python-link]: https://pypi.org/project/md-davis
[commits-badge]: https://img.shields.io/github/commits-since/djmaity/md_davis/v0.3.0.svg
[commits-link]: https://github.com/djmaity/md_davis/compare/v0.3.0...master
