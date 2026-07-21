Install
=======

MD DaVis supports 64-bit Windows, Linux, and macOS with Python 3.12 through
3.14. A Miniforge, Miniconda, or Anaconda installation is required for the
PyMOL and wxPython binary packages used by the complete application.

Runtime installation
--------------------

Clone or download the repository, change to its root directory, and create the
checked-in environment:

.. code-block:: bash

   conda env create --file environment.yml
   conda activate md-davis
   pip install --no-deps --editable .

``--no-deps`` keeps pip from replacing the compatible scientific packages
that Conda has already solved. The environment includes
``pymol-open-source`` from conda-forge; MD DaVis does not use a commercial
PyMOL bundle.

For a custom environment name, pass ``--name``:

.. code-block:: bash

   conda env create --file environment.yml --name my-md-davis
   conda activate my-md-davis

Development installation
------------------------

The development manifest adds tests, linting, documentation, and package
build tools:

.. code-block:: bash

   conda env create --file environment_dev.yml
   conda activate md-davis-dev
   pip install --no-deps --editable .
   pytest

The environment definitions are tested with Conda's dry-run solver on Windows,
Linux, and macOS in continuous integration.

PyMOL
-----

The runtime and development manifests install
`Open-Source PyMOL <https://github.com/schrodinger/pymol-open-source/>`_ from
`conda-forge <https://anaconda.org/conda-forge/pymol-open-source>`_. It can be
added to another Conda environment with:

.. code-block:: bash

   conda install --channel conda-forge "pymol-open-source>=3.1,<4"

PyMOL is imported only when a PyMOL-backed feature is requested. Core parsing,
analysis, and plotting code can therefore be imported without PyMOL.

External programs
-----------------

Several analyses call external scientific programs that are not distributed
with MD DaVis:

GROMACS
^^^^^^^

GROMACS provides trajectory processing, RMSD, radius-of-gyration, secondary
structure, contact, and H-bond calculations. Install a current supported
GROMACS release using its official packages or documentation.

DSSP
^^^^

Secondary-structure calculations require DSSP when invoked through GROMACS.
Configure the ``DSSP`` environment variable if the GROMACS installation does
not discover the executable automatically.

DelPhi and MSMS
^^^^^^^^^^^^^^^

Surface electrostatics require:

* `DelPhi C++ <http://compbio.clemson.edu/delphi>`_ 8.1 or later.
* `MSMS <https://ccsb.scripps.edu/msms/downloads/>`_ 2.6.1 or later.

Pass their executable paths to the ``md-davis electrostatics`` command.

Verification
------------

After installation, verify the package and command entry points:

.. code-block:: bash

   python -c "import md_davis; print(md_davis.__version__)"
   md-davis --help
   md-davis-gui

Uninstall
---------

Remove the complete environment without affecting other Python installations:

.. code-block:: bash

   conda deactivate
   conda env remove --name md-davis
