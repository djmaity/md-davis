Installation
============



Stable release
--------------

To install MD DaVis, run this command in your terminal:

.. code-block:: bash

    pip install md_davis

This is the preferred method to install MD DaVis, as it will always
install the most recent stable release.

If you don't have ```pip`` <https://pip.pypa.io>`_ installed, this
`Python installation guide <http://docs.python-guide.org/en/latest/starting/installation/>`_
can guide you through the process.



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

 Make sure to activate the environment before running the install commands

To install MD DaVis, run this command in your terminal:
```shell
pip install md_davis
```


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

    pip install md_davis

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

Dependencies
------------

.. code-block:: bash

    pip install path/to/md_davis

* [Open-Source PyMOL](https://github.com/schrodinger/pymol-open-source) available from `conda-forge` is for 64-bit Linux and Windows systems only and requires Python > 3.7
Commercial version of pymol can be installed with:
        ```console
        conda install -c schrodinger pymol-bundle
        ```
This can also be used with an Educational PyMOL [license](https://pymol.org/edu/?q=educational)

* mdtraj is available for linux-64, osx-64, win-32, and win-64

Python dependencies are automatically installed. However, electrostatic calculation requires on following two programs which must be downloaded and installed separately.
