Contributing
============

Bug reports and focused pull requests are welcome. Include the operating
system, Python version, MD DaVis version, input format, and a minimal example
when reporting a problem.

Development environment
-----------------------

Clone the repository and create the development environment from the project
root:

.. code-block:: bash

   conda env create --file environment_dev.yml
   conda activate md-davis-dev
   pip install --no-deps --editable .

The environment uses ``pymol-open-source`` from conda-forge. There is no
commercial PyMOL environment variant.

Quality checks
--------------

Run the complete local checks with:

.. code-block:: bash

   pytest
   tox -e lint
   sphinx-build -W --keep-going -b html docs/source docs/_build/html
   python -m build
   twine check dist/*

``tox`` provides equivalent isolated test, lint, documentation, and package
build environments. New behavior should have tests, and bug fixes should add a
regression test that fails without the fix.

Pull requests
-------------

Keep changes focused, use the canonical domain packages, update the
documentation for user-visible behavior, and add a note to ``HISTORY.rst``.
The GitHub Actions matrix validates supported Python versions and performs
Conda solver checks on Windows, Linux, and macOS.
