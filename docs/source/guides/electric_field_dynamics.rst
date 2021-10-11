Electric Field Dynamics
=======================

The electrostatic potentials calculated in :ref:`Surface Electrostatic
Potential Per Residue` can be visualized as a 3D animation of
electric field lines using:

.. code-block:: bash

    md_davis electrodynamics --ss_color --surface --name Human_AcP 2VH7/2VH7_electrostatics

This creates a `PyMOL <https://pymolwiki.org/>`_ session with the
conformations as frames in the animation as shown below:




1. The coordinates of the reference structure are translated to place the
   center of mass of the molecule at the origin and rotated so that the first,
   second, and third principal axes are along the x, y, and z-axes,
   respectively.

2. The frames sampled from the trajectory are aligned to the reference.

3. The electrostatic potentials are obtained for each sampled structure
   using Delphi. The box for each calculation is centered at the origin, and
   the number of grid points is manually set to the same value for each
   structure to ensure the same box size during each calculation.

4. The surface electrostatic potentials calculated per residue or atom are
   written into the output PDB file's B-factor or occupancy column.

5. The output PDB file and the corresponding electric field from the sample
   are visualized as frames in PyMOL (Schrödinger, LLC, 2015), which can
   animate the dynamics of the electric field lines.
