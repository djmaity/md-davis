Install
=======

System Requirements
-------------------

* A 64-bit operating system
* A Python 3 installation with version ≥ 3.7

Conda Installation
------------------

The easiest method to install is with
`Anaconda <https://www.anaconda.com/products/individual>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which works on all operating systems.
It is highly recommended to install MD DaVis in a virtual environment.
The `environment.yml <https://github.com/djmaity/md-davis/blob/master/environment.yml>`_
is provided to ease the process.

.. code-block:: bash

    conda env create -f environment.yml -n md_davis_env

This automatically creates a conda environment called ``md_davis_env`` with all required dependencies.
Activate the environment and install MD DaVis in it using:

.. code-block:: bash

    conda activate md_davis_env
    pip install md-davis

If you don't have `pip <https://pip.pypa.io>`_ installed, this
`Python installation guide <http://docs.python-guide.org/en/latest/starting/installation/>`_
can guide you through the process.

Finally, obtain the `Delphi <http://compbio.clemson.edu/delphi>`_ and
`MSMS <http://mgltools.scripps.edu/downloads#msms>`_ if you intend to perform `electrostatics analysis <user_guides/electrostatics>`_.

Installation in the base conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you do not want to use the `environment.yml <https://github.com/djmaity/md-davis/blob/master/environment.yml>`_
file or if you want to install it in the base environment. Install mdtraj and pymol before installing md-davis.

.. code-block:: bash

    conda install -c conda-forge mdtraj pymol-open-source
    pip install md-davis

Linux Installation
------------------

MD DaVis can be installed on Linux with:

.. code-block:: bash

    pip install md-davis

.. note:: PyMOL is required to run MD DaVis, which may be challenging to install in a virtual environment without ``conda``. (see :ref:`External Dependencies`)

Windows Installation
--------------------

On Windows, ``pip`` may fail to install MD DaVis due to errors with
compiling dependencies. Please follow the instructions under :ref:`Conda
Installation`.

Development Version
-------------------

You can obtain the source code by cloning the public repository:

.. code-block:: bash

    git clone git@github.com:djmaity/md-davis.git

or downloading the `tarball <https://github.com/djmaity/md-davis/tarball/master>`_:

.. code-block:: bash

    curl -OL https://github.com/djmaity/md-davis/tarball/master

Once you have a copy of the source, it can be installed with:

.. code-block:: bash

    python setup.py install

OR

.. code-block:: bash

    pip install path/to/extracted/source/code

You can also install the development version directly with:

.. code-block:: bash

    pip install https://github.com/djmaity/md-davis/archive/master.zip

External Dependencies
---------------------

GROMACS
^^^^^^^

Currently, most analyses have to be performed with GROMACS,
and the output is provided to MD DaVis for visualization.
Other analysis tools may be used as long as the input to MD DaVis can be appropriately formatted.

PyMOL
^^^^^

PyMOL is not available in the `python package index <https://pypi.org/>`_.
Therefore, it cannot be automatically installed with ``pip``.

Open-Source PyMOL
"""""""""""""""""

`Open-Source PyMOL <https://github.com/schrodinger/pymol-open-source/>`_
available from `conda-forge <https://anaconda.org/conda-forge/pymol-open-source>`_
can be installed with:

.. code-block:: bash

    conda install -c conda-forge pymol-open-source

Alternatively, on Linux, PyMOL can be installed with the system package manager, e.g., ``apt`` in Ubuntu or ``dnf`` in Fedora.
However, it is not possible to install PyMOL into a virtual environment using this method.
Therefore, MD DaVis must be installed in the system python as well.

On Windows, if you are not using ``conda``, then pre-built Open-Source PyMOL can be downloaded from
`Christoph Gohlke's page <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pymol-open-source>`_
distributing unofficial windows binaries for python extension packages.

Commercial/Educational PyMOL
""""""""""""""""""""""""""""

`Commercial/Educational PyMOL <https://pymol.org/2/buy.html?q=buy>`_ from
`Schrödinger <https://pymol.org/2/#download>`_ can be installed with:

.. code-block:: bash

    conda install -c schrodinger pymol-bundle


DelPhi and MSMS
^^^^^^^^^^^^^^^

Python dependencies are automatically installed.
However, the electrostatics calculation requires the following two programs,
which must be obtained separately.

* `Delphi <http://compbio.clemson.edu/delphi>`_
* `MSMS <http://mgltools.scripps.edu/downloads#msms>`_

Uninstall
=========

MD DaVis can be easily uninstalled like any other python package, with:

.. code-block:: bash

    pip uninstall md-davis

As with any python package, this does not remove the dependencies installed by MD DaVis.
That is why it is recommended to install MD DaVis in a virtual environment.
Then, the virtual environment may be entirely removed without affecting other python packages on the system.

.. note:: On Linux, if MD DaVis was installed as root or with ``sudo``, the uninstall command should be run with ``sudo``.
