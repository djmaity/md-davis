Plotting xmgrace (.xvg) files
=============================

The GROMACS molecular dynamics software bundles numerous analysis tools. The
output from these are in (.xvg) format to be viewed with xmgrace. Although
xmgrace has some powerful tools for quick calculation of mean standard
deviation, etc., the plots are cumbersome to deal with and difficult to
compare multiple files.

MD DaVis provides the command ``xvg`` to plot such files using plotly or
matplotlib.

It is generally intended that the when plotting multiple xvg files, they
contain the same kind of data. Therefore, the titles and axes labels in the
last supplied file are used.
