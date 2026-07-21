History
=======

Unreleased
----------

* Replace commercial PyMOL distributions with ``pymol-open-source``.
* Modernize packaging, supported Python versions, and Conda environments.
* Consolidate duplicated parsing, plotting, GUI, and PyMOL integration code.
* Reorganize the package into core, analysis, I/O, integration, visualization,
  CLI, and GUI domains and remove the legacy import hierarchy.
* Add unit, command-line, packaging, and environment validation tests.

0.1.0 (2019-09-06)
------------------

* First release on PyPI.

0.2.0 (2020-04-10)
------------------

* Second alpha release on GitHub.

0.3.0 (2021-10-01)
------------------

* First stable feature-rich release

0.4.0 (2022-02-23)
------------------

* GUI for MD DaVis added; invoke with md-davis-gui
* Fixed plotting xvg files without title
* Changed data type of resi from int to int32 in collate.py

0.4.1 (2022-03-02)
------------------

* Fixed installation on Linux due to open-source PyMOL 2.4.0
* Fixing readthedocs documentation compilation
* Bug Fix: landscape_animation was not working due change of landscapes from array to dict
* Removed docopt dependency
* Renamed md_davis to md-davis in documentations

