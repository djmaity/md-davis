Quickstart
==========

Python Module
-------------

To use MD DaVis as Python module:

.. code:: python

    import md_davis

Command Line Interface
----------------------

Getting Help
^^^^^^^^^^^^

The list of available commands can be seen with ``--help`` or ``-h``:

.. code-block:: bash

    md_davis --help

Similarly, the help for each command can be accessed with ``--help`` or ``-h``:

.. code-block:: bash

    md_davis <command> -h

The following steps walks you through the process of using the MD DaVis
command line interface (CLI) to assimilate the analysis data and creating
interactive plots.

Step 1: Calculate the required quantities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MD DaVis has been designed to natively work with the output from GROMACS
analysis tools. However, the same analysis performed with other tools may be
visualized by MD DaVis as long as the input files can be be appropriately
formatted. Root mean squared deviation (RMSD) and radius of gyration (R\
:sub:`G`\ ) will be used to create free energy landscapes and the rest will
be used for residue property plot. Ensure that the trajectories and structures
are processed and corrected for periodic boundary conditions.

Root mean squared deviation (RMSD)
""""""""""""""""""""""""""""""""""

The structure provided will be used as the reference for calculating RMSD of
the trajectory.

.. code-block:: bash

    gmx rms -f trajectory.trr -s structure.pdb -o rmsd.xvg

Radius of gyration (R\ :sub:`G`\)
"""""""""""""""""""""""""""""""""

.. code-block:: bash

    gmx gyrate -f trajectory.trr -s structure.pdb -o rg.xvg

Root mean squared fluctuation (RMSF)
""""""""""""""""""""""""""""""""""""

.. code-block:: bash

    gmx rmsf -f trajectory.trr -s structure.pdb -res -o rmsf.xvg

This would use the whole molecule for RMSF calculation. You may calculate
the RMSF of each chain separately in a protein with multiple chains. This
would reflect the fluctuations only due to backbone motion and eliminate the
fluctuation from relative motion of chains.

Solvent accessible surface area (SASA)
""""""""""""""""""""""""""""""""""""""

To calculate the average SASA per residue pass ``-or`` option to ``gmx sasa``.

.. code-block:: bash

    gmx sasa -f trajectory.trr -s structure.pdb -o sasa.xvg -or resarea.xvg


Secondary structure using DSSP
""""""""""""""""""""""""""""""

.. code-block:: bash

    gmx do_dssp -f trajectory.trr -s structure.pdb -o dssp -ssdump sec_str -sc sec_str_count

``-ssdump`` option will output a file ``sec_str.dat`` containing the secondary
structure for the full trajectory as single letter codes.

+------+----------------------+
| Code | Secondary structure  |
+======+======================+
| H    | α-helix              |
+------+----------------------+
| G    | 3\ :sub:`10`\ -helix |
+------+----------------------+
| I    | π-helix              |
+------+----------------------+
| B    | β-bridge             |
+------+----------------------+
| E    | β strand             |
+------+----------------------+
| T    | Turn                 |
+------+----------------------+
| S    | Bend                 |
+------+----------------------+
| ~    | Loop                 |
+------+----------------------+

This would be processed by MD DaVis to calculate the percentage of occurrence
of each secondary structure at each residue location.

.. note:: For the GROMACS command ``do_dssp`` to work the DSSP binary must
    be available on your system. Download DSSP from
    ftp://ftp.cmbi.ru.nl/pub/software/dssp/

Step 2: Collate the data into an HDF file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Step 2a:** Create a TOML file (input.toml) as shown below specifying the
output files from Step 1.

.. code-block:: toml

    name = 'protein_name_no_spaces'
    output = 'protein_data.h5'
    label = 'Protein from <i>some organism</i>'
    text_label = 'Protein Name'

    trajectory = 'trajectory.xtc'
    structure = 'structure.pdb'

    [timeseries]
        rmsd = 'rmsd.xvg'
        rg = 'rg.xvg'

    [dihedral]
        chunk = 101         # Number of frames to read at a time

    [residue_property]
        secondary_structure = 'sec_str.dat'
        sasa = 'resarea.xvg'

        [residue_property.rmsf]
            rmsf_files = 'rmsf.xvg'
            start = 0               # Starting time for RMSF calculation reference
            end = 100               # End time for RMSF calculation reference

If the ``chunk`` under ``dihedral`` is provided, MD DaVis will calculate
backbone dihedral angles for all frames and the torsional flexibility
(circular standard deviation) of each dihedral angle.

**Step 2b:** Collect all the data from the output files into a single HDF
file using the following command:

.. code-block:: bash

    md_davis collate input.toml

Step 3: Plot the free energy landscape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    md_davis landscape -T 300 --common --select backbone output1.h5 output2.h5 -s landscapes.h5

This command will create an html file with the interactive landscapes. It
will not open the file like other plotting commands, so check the working
directory for the output html file.

Step 4: Plot the residue property plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Step 4a:** Create a pickle file with the residue dataframe using:

.. code-block:: bash

    md_davis residue protein_data.h5 -o residue_dataframe.p

**Step 4b:** Plot the residue data pickle file from the previous command using:

.. code-block:: bash

    md_davis plot_residue residue_dataframe.p -o plot.html



