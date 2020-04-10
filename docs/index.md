# MD Davis
A Python package to perform preliminary analysis of molecular dynamics trajectories. The source code can be found at: https://github.com/djmaity/md_davis

## What can MD Davis do?
1. Create dynamic visualization of electric field around the molecule using PyMOL.
![Dynamics of electric field](2GV1_100frames.gif "Dynamics of electric field")

2. Create interactive plot for time series data: root mean squared deviation (RMSD) and radius of gyration

<iframe src="acylphosphatase_rmsg_rg.html" frameborder="0" width="100%" height="500px"></iframe>

3. Create Free energy landscapes using common ranges and number of bins.

<iframe src="landscapes.html" frameborder="0" width="100%" height=500px></iframe>

4. Create interactive plot of containing the following residue level properties:
    * root mean squared fluctuation (RMSF)
    * circular standard deviation of backbone dihedral angles
    * Secondary structure
    * Solvent accessible surface area
    * Mean and standard deviation for the total surface     electrostatic potential per residue
    * Mean and standard deviation for the mean surface electrostatic potential per residue

It can also use an alignment to align the residue level data from different proteins along the x-axis. This ensures that the peaks line up properly for better interpretation.

<iframe src="acylphosphatase_residue_wise_data_aligned.html" frameborder="0" width="1200px" height=700px></iframe>

Some instructions for creating these types of plots are available here: https://github.com/djmaity/md_davis/instructions.html
help for each command in MD&nbsp;Davis can be accessed with `-h`, for example
```shell
md_davis -h
```