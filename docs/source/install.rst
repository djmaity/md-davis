Install
=======

System Requirements
-------------------

* A 64-bit operating system (We have tested the installation on Ubuntu 20.04,
  Windows 10, and Windows 11, but any modern operating system should work.)
* An `Anaconda <https://www.anaconda.com/products/individual>`_
  or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ Python
  distribution (3.9 or greater)
* An internet connection

Installation Instructions
-------------------------

**Step 1:** Download and install Anaconda
(https://www.anaconda.com/products/individual)
or Miniconda (https://docs.conda.io/en/latest/miniconda.html)
Python distribution.

.. image:: /_static/install/Anaconda.png
   :target: https://www.anaconda.com/products/individual

**Step 2:** Create the `conda environment <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
with MD DaVis dependencies.
Open a Terminal or Anaconda prompt (on Windows) and issue the following command:

.. code-block:: bash

   conda env create djmaity/md-davis

This creates a conda environment called ``md-davis`` with all
required dependencies.

**Step 3:** Activate the ``md-davis`` environment with:

.. code-block:: bash

   conda activate md-davis

**Step 4:** Install MD DaVis in this environment with pip:

.. code-block:: bash

   pip install md-davis

**Step 5:** Obtain the
:ref:`external dependencies <install:External Dependencies>`
as stated below.

Installation for Development
----------------------------

**Step 1:** Download and install Anaconda
(https://www.anaconda.com/products/individual)
or Miniconda (https://docs.conda.io/en/latest/miniconda.html) Python
distribution.

**Step 2:** Create the conda environment with MD DaVis dependencies. In this
case, you should use:

.. code-block:: bash

   conda env create djmaity/md-davis-dev

This would create a conda environment called ``md-davis-env`` with
the core dependencies and install additional libraries for packaging,
documentation, and linting. You can also change the name of the conda
environment by modifying the command above:

.. code-block:: bash

   conda env create djmaity/md-davis-dev --name my_env

This will name the conda environment to ``my_env``. Make sure that
there are no other conda environments with the same name.
The ``djmaity/md-davis-dev`` in the above commands uses the environment files
uploaded to https://anaconda.org/djmaity/environments. These are the same files
`environment.yml <https://github.com/djmaity/md-davis/blob/master/conda_environments/environment.yml>`_
and `dev_environment.yml <https://github.com/djmaity/md-davis/blob/master/conda_environments/environment.yml>`_
in the source code, which may also be used to create the conda environment and
install the dependencies:

.. code-block:: bash

   conda env create --file dev_environment.yml

See https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
to learn more about conda environments.

.. note::
   We recommend using the environment files over installing the dependencies
   individually because installing all the dependencies together allows conda
   to resolve the appropriate version of packages to install and avoid
   package conflicts.

**Step 3:** Activate the conda environment with:

.. code-block:: bash

   conda activate md-davis-dev

Remember to change ``md-davis-dev`` to the appropriate environment name if
you modify it.

**Step 4:** Obtain a copy of the source code by cloning the public repository:

.. code-block:: bash

   git clone git@github.com:djmaity/md-davis.git

Or, download the `tarball <https://github.com/djmaity/md-davis/tarball/master>`_:

.. code-block:: bash

   curl -OL https://github.com/djmaity/md-davis/tarball/master

**Step 5:** Once you have a copy of the source, it can be installed with:

.. code-block:: bash

   pip install path/to/source/code

The path to the directory containing the ``setup.py`` file has to be provided in
the command above. You may want to install it as an editable package with:

.. code-block:: bash

   pip install -e path/to/source/code

This will allow changes to the source code to be immediately reflected in the
installed package.

Or, go to the directory containing the ``setup.py`` file and use:

.. code-block:: bash

   python setup.py install

You can also install the development version directly with:

.. code-block:: bash

   pip install https://github.com/djmaity/md-davis/archive/master.zip

**Step 6:** Obtain the
:ref:`external dependencies <install:External Dependencies>`
as stated below.

External Dependencies
---------------------

GROMACS
^^^^^^^

Currently, most analyses have to be performed with GROMACS, and the output is
provided to MD DaVis for visualization. We have successfully used MD DaVis on
simulations and analysis performed with various GROMACS versions from 5.1.4
to 2021.3. Other analysis tools may be used as long as the input to MD DaVis
can be appropriately formatted. See
:ref:`introduction:Interfacing MD DaVis to Other Analysis Tools`.

DSSP
""""

The secondary structure of a trajectory is calculated by the GROMACS tool
``do_dssp``, which requires `DSSP <https://github.com/cmbi/dssp>`_.
The latest and only available version of
`DSSP is 2.3.0 <https://github.com/cmbi/dssp>`_.
The executable is called ``mkdssp`` now. Help GROMACS find it by any of the
following means:

* Rename ``mkdssp`` to ``dssp``
* Make a symlink called ``dssp`` to ``mkdssp``
* Set the ``DSSP`` environment

PyMOL
^^^^^

You do not need to obtain PyMOL separately if you use the installation methods
outlined above using ``conda``. Unfortunately, PyMOL is not available in the
`python package index <https://pypi.org/>`_. Therefore, it cannot be
automatically installed with ``pip``. However, Open-Source PyMOL and
Commercial/Educational PyMOL are available in conda channels
`conda-forge <https://anaconda.org/conda-forge/pymol-open-source>`_
and `schrodinger <https://anaconda.org/schrodinger>`_, respectively.

Open-Source PyMOL
"""""""""""""""""

The command in **Step 2** automatically installs
`Open-Source PyMOL <https://github.com/schrodinger/pymol-open-source/>`_
available from
`conda-forge <https://anaconda.org/conda-forge/pymol-open-source>`_.
It can also be installed with:

.. code-block:: bash

   conda install -c conda-forge pymol-open-source=2.5.0

.. warning::
   Open-Source PyMOL 2.4.0 has a bug where it cannot open Gaussian cube files.
   DelPhi to output phimap volumetric data is in this format, which
   is used in the :ref:`electrostatics <guides/electrostatics:Surface Electrostatic Potential Per Residue>`
   calculations.

Alternatively, on Linux, PyMOL can be installed with the system package
manager, e.g., ``apt`` in Ubuntu or ``dnf`` in Fedora. However, it is not
possible to install PyMOL into a virtual environment using this method.
Therefore, MD DaVis must be installed in the system python, which may
create conflicts with existing python packages on which many system programs
may depend.

On Windows, if you are not using ``conda``, then pre-built Open-Source PyMOL
can be downloaded from
`Christoph Gohlke's page <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pymol-open-source>`_
distributing unofficial windows binaries for python extension packages.
However, we have faced issues with using this package with other dependencies
of MD DaVis.

Commercial/Educational PyMOL
""""""""""""""""""""""""""""

The `Commercial/Educational PyMOL <https://pymol.org/2/buy.html?q=buy>`_
can be used instead of Open-Source PyMOL by installing the
`md-davis-pymol <https://anaconda.org/djmaity/md-davis-pymol>`_ or
`md-davis-pymol-dev <https://anaconda.org/djmaity/md-davis-pymol-dev>`_
environment in **Step 2** of the above-mentioned installation methods:

.. code-block:: bash

   conda env create djmaity/md-davis-pymol

Alternatively, it can be installed with:

.. code-block:: bash

   conda install -c schrodinger pymol-bundle

DelPhi and MSMS
^^^^^^^^^^^^^^^

Python dependencies are automatically installed. However, the electrostatics
calculation requires the following two programs, which must be obtained
separately.

* `Delphi C++ <http://compbio.clemson.edu/delphi>`_ version greater than or
  equal to 8.1
* `MSMS <https://ccsb.scripps.edu/msms/downloads/>`_

Uninstall
=========

MD DaVis can be easily uninstalled like any other python package, with:

.. code-block:: bash

   pip uninstall md-davis

As with any python package, this does not remove the dependencies installed
by MD DaVis. That is why installing MD DaVis in a conda environment is
recommended. Then, the whole conda environment may be entirely removed without
affecting other python packages on the system.

.. code-block:: bash

    conda env remove --name md-davis

where ``md-davis`` is the name of the conda environment. Modify the command
to provide the appropriate name for the conda environment if you change it.

.. note::
   On Linux, if MD DaVis was installed as root or with ``sudo`` (highly
   discouraged), the uninstall command should also be run with ``sudo``.
