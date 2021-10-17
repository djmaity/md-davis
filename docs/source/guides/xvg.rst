.. click:: md_davis.xvg:main
   :prog: md_davis xvg
   :nested: full

The GROMACS molecular dynamics software bundles numerous analysis tools. The
output from these are in (.xvg) format to be viewed with xmgrace. Although
xmgrace has some powerful tools for quick calculation of mean standard
deviation, etc., the plots are cumbersome to deal with and difficult to
compare multiple files.

MD DaVis provides the command ``xvg`` to plot such files using plotly or
matplotlib.

The .xvg files are
space or tab delimited text files with time and data in the first and second
columns, respectively.

.. code-block::

    # gmx rms is part of G R O M A C S:
    #
    # Gromacs Runs On Most of All Computer Systems
    #
    @    title "RMSD"
    @    xaxis  label "Time (ps)"
    @    yaxis  label "RMSD (nm)"
    @TYPE xy
    @ subtitle "System after lsq fit to System"
       0.0000000    0.0000031
      10.0000000    0.0677038
      20.0000000    0.0837483
      30.0000000    0.0894995



It is generally intended that the when plotting multiple xvg files, they
contain the same kind of data. Therefore, the titles and axes labels in the
last supplied file are used.

MD DaVis can also plot multiple Grace (.xvg) files, which is the format for the output files from many GROMACS analysis tools. For example, an interactive plot with RMSD and RG from multiple trajectories can be created for quick comparison.

.. code-block:: bash

    md_davis xvg <path/to/file.xvg>

Replace `<path/to/file.xvg>` with the location of your ``.xvg`` file.

Create interactive plot for time series data: root mean squared deviation (RMSD) and radius of gyration (R\ :subscript:\ G)

<iframe src="/_static/acylphosphatase_rmsg_rg.html" frameborder="0" width="100%" height="500px"></iframe>
