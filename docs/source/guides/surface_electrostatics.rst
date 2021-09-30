Surface Electrostatic Potential Per Residue
===========================================

Usually, the distribution of electrostatic potential around a protein mole-cule is visualized by solving the Poisson-Boltzmann equation and color-ing the molecular surface with obtained electrostatic potentials, which can be compared qualitatively. MD DaVis provides a wrapper to calcu-late the electrostatic potential on the vertices of a triangulated molecular surface and extract the total and mean surface electrostatic potential per residue, enabling a quantitative comparison. DelPhi versions later than 8.0 (Li et al., 2019) can calculate the electrostatic potential at user-specified points, which was leveraged to calculate the surface electrostat-ic potential at each residue. The steps involved are as follows:
1.	Calculate the triangulated molecular surface using the MSMS program (Sanner et al., 1996).
2.	The electrostatic potentials at these vertices of the triangulated surface are calculated using DelPhi.
3.	The vertices corresponding to each residue are identified along with the electrostatic potentials at these points.
4.	The total and mean of the electrostatic potentials on the vertices for each residue are calculated.
The total gives information about the charge on the surface as well as the exposed surface area, while the mean nullifies the contribution from the exposed surface area. This calculation is repeated for every confor-mation in a sample from the trajectory with tens or hundreds of frames. The mean and the standard deviation of both total surface potential and the mean surface potential are included in the overlaid residue plot creat-ed by MD DaVis. Using a sample from the MD for the electrostatics calculations incorporates the dynamical information into the plot. There-fore, it has greater confidence compared to inferences based on a single structure or frame. It is advisable to align the frames and use the same grid size for each Delphi calculation. However, the variation in the mean and total surface electrostatic potentials due to the different orientations of the protein molecule in the grid is negligible as long as the triangulated molecular surface is high resolution.
The default parameters are used for running Delphi with 2000 linear interations and a maximum change in potential of 10-10 kT/e.  The salt concentration is set to 0.15 M, and a solvent dielectric value of 80 (die-lectric of water) was used. The atomic charges and radii from the CHARMM force field are provided by default. The output potential map is saved in the cube file format used by Gaussian software; the other formats do not load in molecular visualization software. A discrepancy between the molecular surface calculated by MSMS and the molecular surface used by DelPhi to detect the protein’s interior was observed for most buried residues in multimeric proteins. DelPhi uses a different die-lectric constant for the molecule’s interior, so the potential on these resi-dues was unreasonably high. This caveat should be borne in mind while analyzing the results of these calculations.


Delphi and MSMS must be installed for this to work. You may download these from: [Delphi](http://compbio.clemson.edu/delphi)  and [MSMS](http://mgl.scripps.edu/people/sanner/html/msms_home.html)

You will need the parameter files for delphi run, which can be downloaded from [paramter files](http://compbio.clemson.edu/downloadDir/delphi/parameters.tar.gz)

**Step 1:** Calculate the electrostatics using the following command

.. code-block:: bash

    md_davis electrostatics --surface_potential \
        -m ~/.opt/msms/ \
        -d ~/.opt/delphi8.4.2_sequential
        -c path/to/charge/file.crg \
        -r path/to/radius/file.siz \
        PDB_FILE.pdb \
        path/to/electrostatics/output/directory

The following charge and radius files are provided with the source code for MD&nbsp;DaVis.

.. code-block:: bash

    md_davis/md_davis/electrostatics/charmm.crg
    md_davis/md_davis/electrostatics/charmm.siz

If you receive a warning during Delphi run regarding missing charge or radius. Then the missing properties must be added to these files or whichever files you provide to `md_davis electrostatics`.

**Step 2:** Add the surface electrostatics data into the residue dataframe by modifying the command for creating residue dataframe to include the option `-d` and specifying the path to the directory containing the surface potential files. This will search for all  `.pot` files and include their mean and total in the residue dataframe.

.. code-block:: bash

    md_davis residue dataframe \
            --potential  path/to/electrostatics/output/directory \
            --annotation  annotations.json \
            --prefix  prefix \
            md_davis_collect_output_data.h5  \
            output_residue_wise_data.p

**Step 3:** Plot the residue dataframe with the usual command:

.. code-block:: bash

    md_davis plot residue output_residue_wise_data.p


The electrostatic potentials calculated for a sample of frames from the MD trajectories, as discussed above, can be visualized as a 3D animation of electric field lines (Figure 1C). This highlights the effect of dynamics of the electric field in the vicinity of the molecule. The 3D electric field dynamics is visualized as described below:
1.	The coordinates of the reference structure are translated to place the center of mass of the molecule at the origin and rotated so that the first, second, and third principal axes are along the x, y, and z-axes, respectively.
2.	The frames sampled from the trajectory are aligned to the refer-ence.
3.	The electrostatic potentials are obtained for each sampled struc-ture using Delphi. The box for each calculation is centered at the origin, and the number of grid points is manually set to the same value for each structure to ensure the same box size during each calculation.
4.	The surface electrostatic potentials calculated per residue or atom are written into the output PDB file’s B-factor or occupancy col-umn.
5.	The output PDB file and the corresponding electric field from the sample are visualized as frames in PyMOL (Schrödinger, LLC, 2015), which can animate the dynamics of the electric field lines.
