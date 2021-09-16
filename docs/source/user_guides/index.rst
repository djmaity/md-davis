User Guide
==========

.. toctree::
   :maxdepth: 2

   quick_start
   landscape
   residue_property_plot
   secondary_structure
   electrostatics
   contacts
   xvg
   sequence

MD Davis is a Python package to perform exploratory and comparative data
analysis of molecular dynamics trajectories. The source code can be found
at: https://github.com/djmaity/md_davis

a python package to easily create helpful interactive visualizations to compare and analyze multiple molecular dynamics trajectories of proteins (https://github.com/djmaity/md_davis). The tool can increase the productivity of researchers using molecular dynamics simulations and make the analysis of such simulations accessible to other researchers.

The package provides the commandline tool ``md_davis`` for creating common visualizations and

Features of MD DaVis
--------------------

1. Create dynamic visualization of electric field around the molecule using PyMOL.

.. image:: /_static/2VH7_electrodynamics.webp
     :alt: Dynamics of electric field

2. Create interactive plot for time series data: root mean squared deviation (RMSD) and radius of gyration (R\ :subscript:\ G)

<iframe src="/_static/acylphosphatase_rmsg_rg.html" frameborder="0" width="100%" height="500px"></iframe>

3. Create Free energy landscapes using .

<iframe src="/_static/landscapes.html" frameborder="0" width="100%" height=500px></iframe>

It can also use an alignment to align the residue level data from different proteins along the x-axis. This ensures that the peaks line up properly for better interpretation.

<iframe src="/_static/acylphosphatase_residue_wise_data_aligned.html" frameborder="0" width="1200px" height=700px></iframe>

Some instructions for creating these types of plots are available here: [Instructions](http://djmaity.com/md_davis/instructions.html) and instructions for calculating surface electrostatic potentials can be found here: [Electrostatics](http://djmaity.com/md_davis/electrostatics.html)

help for each command in MD Davis can be accessed with ``-h``, for example

.. code-block:: bash

    md_davis -h

## Usage

To use MD DaVis in a project

.. code:: python

    import md_davis

<!-- ## Plot .xvg file
To use this script type the following command in a terminal or command prompt and press 'Enter':
.. code-block:: bash

    python plot_xvg.py <path/to/file.xvg>

Replace `<path/to/file.xvg>` with the location of your ``.xvg`` file.

## Plot DSSP file obtained from GROMACS
To obtain the input file for the script, run **do_dssp** command from GROMACS with the **-ssdump** option:

.. code-block:: bash

    gmx do_dssp -f <trajectory> -s <structure> -o <ss.xpm> -ssdump <ssdump.dat>

Use the following script to count and plot the percentage secondary structure per residue throughout the trajectory:

.. code-block:: bash

    python plot_do_dssp_per_residue.py <path/to/ssdump.dat>

Replace `<path/to/ssdump.dat>` depending on the location of your ``.dat`` file. This will show the plot on screen. To output the plot as an image instead of displaying it on screen use:

.. code-block:: bash

    python plot_do_dssp_per_residue.py <path/to/ssdump.dat> -o image.png

### Getting Help
Providing **-h** option to each python script will print out its help message.

.. code-block:: bash

    python <script_name.py> -h

NOTE: _Please replace the text wtih angular brackets < > by the respective filename or path._

# README #

These scripts extract dihedral angles from GROMACS (molecular dynamics) trajectory and analyse enregy landscapes.

### What is this repository for? ###

* Analysing GROMACS trajectories
* 0.1


* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines
