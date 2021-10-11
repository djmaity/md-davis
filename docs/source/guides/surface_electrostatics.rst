Surface Electrostatic Potential Per Residue
===========================================

Usually, the distribution of electrostatic potential around a protein
molecule is visualized by solving the Poisson-Boltzmann equation using a
grid base approach (`APBS <https://www.poissonboltzmann.org/>`_/`Delphi
<http://compbio.clemson.edu/delphi>`_) and coloring the molecular
surface with the obtained electrostatic potentials, as shown below.

.. image:: /_static/surface_electrostatics.png
   :align: center

This only allows for qualitative inspection of one structure or conformation
at a time. In contrast, MD DaVis can calculate the surface electrostatic
potential per residue from electrostatics calculations on a sample of
conformation from the MD trajectory. Thereby, incorporating dynamical
information into the analysis. The steps involved are as follows:

#.  Obtain a sample of conformations from the trajectory as PDB files.
    Ideally, all the conformations should be centered at the origin and
    aligned. Optionally, the reference conformation used for alignment may be
    oriented with its first, second and third `principal axis <https://en
    .wikipedia.org/wiki/Moment_of_inertia#Principal_axes>`_
    along the x, y, and z axis, respectively. When comparing samples from
    two or more trajectories ensure that the reference structures for each
    are aligned to reduce discrepancy caused by grid selection during
    electrostatics calculation.

#.  Use :ref:`md_davis electrostatics` to run MSMS and DelPhi on each PDB file.

    .. code-block::

         md_davis electrostatics --surface -m <MSMS> -d <DELPHI> -o <OUTPUT_DIR> [PDB_FILES]

    This calculates the triangulated molecular surface using the
    `MSMS <http://mgl.scripps.edu/people/sanner/html/msms_home.html>`_ program
    (`Sanner et al., 1996 <https://doi.org/10.1002/(SICI)1097-0282(199603)
    38:3%3C305::AID-BIP4%3E3.0.CO;2-Y>`_) and evaluates the electrostatic
    potentials at the vertices of the triangulated surface
    using `Delphi <http://compbio.clemson.edu/delphi>`_
    (`Li et al., 2019 <https://doi.org/10.1002/jcc.26006>`_). In addition,
    the volumetric data for electrostatic potential at grid points around the
    molecule are also saved by DelPhi in
    `Gaussian CUBE format <https://gaussian.com/cubegen/>`_. This can be
    used for 3D vislization of the surface electrostatic potentials (Figure
    above) or :ref:`electric field dynamics <Electric Field Dynamics>`.
    It is advisable to align the PDB files and specify the grid size for
    Delphi calculation, so that the same grid points are used for each
    calculation.

#.  Specify the path to the directory containing the surface potential
    files in the input `TOML <https://toml.io/en/>`_ file used by
    :ref:`md_davis collate`.

    .. code-block:: toml

        [residue_property]
            surface_potential = '<OUTPUT_DIR>'

#.  The total and mean electrostatic potential per residue is calculated
    when the data is :ref:`collated <md_davis collate>` into an HDF file.

    .. code-block::

        md_davis collate input.toml

    This will search for all ``.pot`` files in the specified
    ``<OUTPUT_DIRECTORY>``. The vertices corresponding to each residue are
    identified along with the corresponding electrostatic potential. The
    total and mean of the electrostatic potentials on the vertices for each
    residue are calculated. The total gives information about the charge on
    the surface as well as the exposed surface area, while the mean
    nullifies the contribution from the exposed surface area. The mean and
    standard deviation of both are calculated and saved in the
    output HDF file.

#.  Follow the steps for creating a :ref:`residue property plot <Residue
    Property Plot>` to visualize the data.

    .. code-block::

        md_davis plot residue output_residue_wise_data.p

.. note:: A discrepancy between the molecular surface calculated by MSMS and
    the molecular surface used by DelPhi to detect the protein's interior
    was observed for buried residues in multimeric proteins. DelPhi uses a
    different dielectric constant for the molecule's interior, so the potential
    on these residues was very high.
