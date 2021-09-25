Residue Property Plot
=====================

Understanding the function–dynamics relationship of protein often re-quires comparing two or more properties, e.g., the effect of dihedral fluctuation on the RMSF of the protein. The two properties, dihedral fluctuation and RMSF, can be easily plotted with any plotting program. Later, one might be interested in knowing the trend between the RMSF and solvent exposure of the residues. Traditionally, each time new prop-erties are compared, new plots need to be prepared. This repetitive pro-cess is alleviated by making overlaid plots of all the residue properties at once (Figure 1A). MD DaVis can plot the following quantities: RMSF, torsional flexibility, secondary structure, solvent accessible surface area, and surface electrostatic potential on each residue. Showing all the prop-erties as an overlaid plot can be overwhelming and cluttering. However, the option to interactively turn on or off the visualization for each data series clears the clutter and highlights the similarities and differences. Moreover, the labels and annotations that appear on hovering the cursor over certain regions improve the interpretability of these plots, thereby granting substantial time savings.
The data obtained from multiple trajectories can be overlaid in a single interactive plot for ease of comparison, which removes the need for making multiple plots and immediately brings out the distinguishing features between the trajectories. The novel feature of MD DaVis is that it can align the data for similar proteins by inserting required gaps along the x-axis using an alignment file. Aligning the residues on the x-axis aligns the peaks, highlighting the similarities between the datasets, which is incredibly powerful when comparing the dynamical information from similar proteins.

4. Create interactive plot of containing the following residue level properties:
    * root mean squared fluctuation (RMSF)
    * torsional flexibility (circular standard deviation of backbone dihedral angles)
    * Secondary structure
    * Solvent accessible surface area
    * Mean and standard deviation for the total surface     electrostatic potential per residue
    * Mean and standard deviation for the mean surface electrostatic potential per residue


.. note:: the paths in the input toml file is relative to the location where the md_davis command will be called from. To avoid any confusion try using absolute paths.

How to interact with the plot
-----------------------------
