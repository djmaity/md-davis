MD DaVis
========

[![Supported versions][python-badge]][python-link]
[![PyPI Package latest release][pypi-badge]][pypi-link]
[![Documentation Status][docs-badge]][docs-link]
[![PyPI Wheel][wheel-badge]][wheel-link]
[![MIT License][license-badge]][license-link]
[![Commits since latest release][commits-badge]][commits-link]
[![DOI][zenodo-badge]][zenodo-link]

Introduction
------------
MD DaVis is a tool for comparative analysis of molecular dynamics simulations of proteins.

Documentation: [https://md-davis.readthedocs.io](https://md-davis.readthedocs.io)

Features
--------
1. Free energy landscape
2. Residue properties plot
3. Surface electrostatics
4. Electric field dynamics
5. H-bond/Contact matrix

System Requirements
-------------------

* A 64-bit operating system
* A Python 3 installation with version ≥ 3.9

Installation
------------
The easiest installation method is with
[Anaconda](https://www.anaconda.com/products/individual) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html).
Create a conda environment called `md-davis` with all required dependencies using the following command:
```
conda env create djmaity/md-davis
```
Activate the environment with:
```
conda activate md-davis
```
Install MD DaVis in this environment using:
```
pip install md-davis
```

For detailed installation instructions, see:
[https://md-davis.readthedocs.io/en/latest/install.html](https://md-davis.readthedocs.io/en/latest/install.html)

### Dependencies

The following tools required for full functionality will have to be obtained separately.
* [Delphi v8.1 or later](http://compbio.clemson.edu/delphi)
* [DSSP v2.3.0](https://github.com/cmbi/dssp)
* [MSMS v2.6.1](https://ccsb.scripps.edu/msms/downloads/)
* [GROMACS v5.1.4 or later](https://www.gromacs.org)

Usage
-----

Remember to activate the `md-davis` environment before using MD DaVis.
```
conda activate md-davis
```

The MD DaVis graphical user interface can be invoked with:
```
md-davis-gui
```
The MD DaVis command line tool can be called with:
```
md-davis
```
The MD DaVis package can be used in a Python script with a standard import
statement like:
```
import md_davis
```

Author
------
Dibyajyoti Maity - [www.djmaity.com](http://www.djmaity.com)

<!--  --------------------------------------------------------------------- -->
[pypi-badge]: https://img.shields.io/pypi/v/md-davis.svg
[pypi-link]: https://pypi.org/project/md-davis
[docs-badge]: https://readthedocs.org/projects/md-davis/badge/?version=latest
[docs-link]: https://md-davis.readthedocs.io/en/latest/
[wheel-badge]: https://img.shields.io/pypi/wheel/md-davis.svg
[wheel-link]: https://pypi.org/project/md-davis
[python-badge]: https://img.shields.io/pypi/pyversions/md-davis.svg
[python-link]: https://pypi.org/project/md-davis
[commits-badge]: https://img.shields.io/github/last-commit/djmaity/md-davis
[commits-link]: https://github.com/djmaity/md-davis/
[license-badge]: https://img.shields.io/pypi/l/md-davis?color=success
[license-link]: https://github.com/djmaity/md-davis/blob/master/LICENSE
[zenodo-badge]: https://zenodo.org/badge/186578728.svg
[zenodo-link]: https://zenodo.org/badge/latestdoi/186578728
