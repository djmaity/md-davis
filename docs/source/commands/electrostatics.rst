md_davis electrostatics
=======================

.. code-block:: bash

    md_davis electrostatics --surface -m <MSMS_EXECUTABLE> -d <DELPHI_EXECUTABLE> -o <OUTPUT_DIRECTORY> [PDB_FILES]

``md_davis electrostatics`` is a wrapper for running
`Delphi <http://compbio.clemson.edu/delphi>`_ and reporting
the electrostatic potential at the vertices of a triangulated surface obtained using
`MSMS <http://mgltools.scripps.edu/downloads#msms>`_. Therefore, these must
be installed for this command to work, which may be downloaded from
http://compbio.clemson.edu/delphi and
http://mgltools.scripps.edu/downloads#msms, respectively.

The default parameters used for running Delphi are:

* 2000 linear interations
* maximum change in potential of 10\ :sup:`-10` kT/e.
* The salt concentration is set to 0.15 M
* solvent dielectric value of 80 (dielectric of water)

If you want to use different set of parameters, you may run Delphi without
the md

The atomic charges and radii from the CHARMM force field are provided by
default. The parameter files containing the atomic charges and radii for
DelPhi are available at:
http://compbio.clemson.edu/downloadDir/delphi/parameters.tar.gz

The output potential map is saved in the cube file format (used by Gaussian
software); the other formats do not load in molecular visualization software.



The following charge and radius files are provided with the source code for
MD DaVis.

    ``md_davis/md_davis/electrostatics/charmm.crg``
    ``md_davis/md_davis/electrostatics/charmm.siz``

If you receive a warning during Delphi run regarding missing charge or
radius. Then the missing properties must be added to these files or
whichever files you provide to ``md_davis electrostatics``.
