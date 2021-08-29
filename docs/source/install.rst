Installation
============

Stable release
--------------

To install MD DaVis, run this command in your terminal:

.. code-block:: bash

    pip install md-davis

This is the preferred method to install MD DaVis, as it will always
install the most recent stable release.

If you don't have ```pip`` <https://pip.pypa.io>`_ installed, this
`Python installation guide <http://docs.python-guide.org/en/latest/starting/installation/>`_
can guide you through the process.

Windows installation
^^^^^^^^^^^^^^^^^^^^

On Windows, ``pip`` may fail to install MD DaVis

when it tries to compile numpy while installing mdtraj

due to errors with compiling some dependencies. Please use the
`Anaconda <https://www.anaconda.com/products/individual>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ python
distribution to install the following dependencies:

.. code-block:: bash

    conda install -c conda-forge mdtraj pymol-open-source

Before running the next command to install MD DaVis

.. code-block:: bash

    pip install md-davis

Installation from source
------------------------

The sources for MD DaVis can be downloaded from the `Github
repo <https://github.com/djmaity/md_davis>`_.

You can either clone the public repository:

.. code-block:: bash

    git clone git://github.com/djmaity/md_davis

Or download the
`tarball <https://github.com/djmaity/md_davis/tarball/master>`_:

.. code-block:: bash

    curl  -OL https://github.com/djmaity/md_davis/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: bash

    python setup.py install

You can also install the development version with

.. code-block:: bash

    pip install https://github.com/djmaity/md_davis/archive/master.zip

OR follow the steps below:

Step 1: Clone this repository using:

.. code-block:: bash

    git clone git@github.com:djmaity/md_davis.git

Step 2: I highly recommend installing this package in a virtual
environment to avoid dependency conflicts with installed packages.
Create a virtual environment in
`Anaconda <https://www.anaconda.com/products/individual>`_ or
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
distribution using:

.. code-block:: bash

    conda create --name md_davis python-3

Activate the environment before running the install command

.. code-block:: bash

    conda activate md_davis

Install the dependencies using the following commands

.. code-block:: bash

    conda install scipy
    conda install psutil
    conda install -c schrodinger pymol
    conda install -c conda-forge mdtraj
    conda install -c plotly plotly-orca
    conda update -c conda-forge h5py

Install this package using pip:

.. code-block:: bash
    
    pip install path/to/md_davis