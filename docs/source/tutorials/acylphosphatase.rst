Acylphosphatase Homologs
========================

In this tutorial we shall analyze the molecular dynamics (MD) trajectories of four acylphosphatase (AcP) homologs:

+------------+-----------+
| PDB ID     | Organism  |
+============+===========+
| 2VH7       | Human     |
+------------+-----------+
| 2GV1       | E\. coli  |
+------------+-----------+
| 2ACY       | Bovine    |
+------------+-----------+
| 1URR       | Fruit Fly |
+------------+-----------+

If you are new to running MD simulations using GROMACS, please refer to the legendary GROMACS tutorial `Lysozyme in Water <http://www.mdtutorials.com/gmx/lysozyme/index.html>`_ by Dr. Justin A. Lemkul at http://www.mdtutorials.com/

RMSD and R\ :sub:`G`\
---------------------
Often the first step after a successful MD simulation is to calculate the root-mean-square deviation (RMSD) and radius of gyration (R\ :sub:`G`\ ). In GROMACS, these can be calculated using `gmx rms`_ and `gmx gyrate`_, respectively.

.. code-block:: bash

    gmx rms -f 2VH7_trajectory.xtc -s 2VH7_structure.pdb -o 2VH7_rmsd
    gmx gyrate -f 2VH7_trajectory.xtc -s 2VH7_structure.pdb -o 2VH7_rg

2VH7_rmsd.xvg and 2VH7_rg.xvg are text files containing the RMSD and R\ :sub:`G`\ , respectively. Repeat the process for all the other trajectories.

.. note:: If you have installed MD DaVis in a virtual or conda environment as suggested in the installation instructions, make sure to activate it before running the ``md_davis`` commands.

Plot 2VH7_processed_rmsd.xvg using:

.. code-block:: bash

    md_davis xvg 2VH7_rmsd.xvg -o 2VH7_rmsd.html

You should obtain a plot like this:

To plot the RMSD from the four trajectories together use:

.. code-block:: bash

    md_davis xvg 2VH7_rmsd.xvg 2GV1_rmsd.xvg 2ACY_rmsd.xvg 1URR_rmsd.xvg -o AcP_rmsd.html

.. image:: /_static/AcP_rmsd.png
    :target: AcP_rmsd.html

Similarly, the R\ :sub:`G`\  can also be plotted using the ``md_davis xvg`` command.

Create Free Energy Landscapes
-----------------------------
Next we will create the free energy landscape using RMSD and R\ :sub:`G`\  as the `x` and `y` variables. The sample trajectories provided with this tutorial only contain 100 frames to keep their sizes small. Thus, the RMSD and R\ :sub:`G`\  files created in the beginning of this tutorial only contain 100 time steps each. Using 100 points to create a free energy landscape would not be accurate, therefore, the RMSD and R\ :sub:`G`\  files calculated from the full 1 μs trajectory containing 100,000 frames are provided with the tutorial. Use these files to create the  free energy landscape. For human acylphosphstase the free energy landscape can be created with:

.. code-block:: bash

    md_davis landscape_xvg -T 300 -x 2VH7_rmsd_full.xvg -y 2VH7_rg_full.xvg -n "2VH7" -l "Human AcP" -o 2VH7_landscape.html

Here, ``-T 300`` specifies 300 K as the temperature of the system. Now to plot all four landscapes together:

.. code-block:: bash

    md_davis landscape_xvg -T 300 --common -x 2VH7_rmsd_full.xvg -y 2VH7_rg_full.xvg --name "2VH7" --label "Human AcP" -x 2GV1_rmsd_full.xvg -y 2GV1_rg_full.xvg --name "2GV1" --label "E. coli AcP" -x 2ACY_rmsd_full.xvg -y 2ACY_rg_full.xvg --name "2ACY" --label "Bovine AcP" -x 1URR_rmsd_full.xvg -y 1URR_rg_full.xvg --name "1URR" --label "Fruit Fly AcP" -o AcP_FEL.html

In the command above, remember to provide ``-x``, ``-y``, ``--name``, and  ``--label`` together before those for the next trajectory. The option ``--common`` instructs MD DaVis to create the four landscapes using identical ranges and binning. This allows us to reliably compare the landscapes. The output from the above command is shown below, click the image to view the interactive html file.

.. image:: /_static/AcP_FEL.png
   :target: AcP_FEL.html

Surface Elecrtostatic Using Delphi
----------------------------------

1. Create a sample of frames for PBSA calculation using Delphi


.. code-block:: bash

    gmx trjconv -f 2VH7/2VH7_trajectory.xtc -s 2VH7/2VH7_structure.pdb -o 2VH7/2VH7_electrostatics/2VH7_frame.pdb -dt 10000 -sep


.. code-block:: bash

    md_davis electrostatics -m msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1 2VH7/2VH7_structure.pdb ./ -d ./delphicpp_v8.4.5_serial


Residue Properties Plot
-----------------------

1. Calculate the root-mean-square fluctuation, solvent accessible surface area, and secondary structure using GROMACS:

.. code-block:: bash

    gmx rmsf -res -f 2VH7_trajectory.xtc -s 2VH7_structure.pdb -o 2VH7_rmsf
    gmx sasa -f 2VH7_trajectory.xtc -s 2VH7_structure.pdb -o 2VH7_sasa.xvg -or 2VH7_resarea.xvg
    gmx do_dssp -f 2VH7_trajectory.xtc -s 2VH7_structure.pdb -o 2VH7_dssp -ssdump 2VH7_dssp -sc 2VH7_dssp_count

Repeat for the remaining trajectories. We are also going to plot the torsional flexibility, but that will be calculated by MD DaVis later.

.. note:: For the ``gmx do_dssp`` command to work the ``dssp`` or ``mkdssp`` binary must be available on your system. Download it from ftp://ftp.cmbi.ru.nl/pub/software/dssp/ and ensure GROMACS can find it by setting the ``DSSP`` environment variable to point to its location on your system.

2. Collect and store all the calculated properties into a HDF file. To do that first create a TOML file as shown below telling MD DaVis the location of each file

.. code-block:: toml

    name = '2VH7'
    output = '2VH7_data.h5'
    label = 'Human AcP'
    text_label = 'Human AcP'

    trajectory = '2VH7_trajectory.xtc'
    structure = '2VH7_structure.pdb'

    [timeseries]
        rmsd = '2VH7_rmsd_full.xvg'
        rg = '2VH7_rg_full.xvg'

    [dihedral]
        chunk = 101

    [residue_property]
        secondary_structure = '2VH7_dssp.dat'
        sasa = '2VH7_resarea.xvg'
        surface_potential = '2VH7/2VH7_electrostatics'

        [residue_property.rmsf]
            rmsf_files = '2VH7_rmsf.xvg'
            start = 0
            end = 100

Input toml files for each trajectory is provided with the tutorial files.

Next, collate all the data using MD DaVis. MD DaVis can process multiple such TOML file at once for each system and create the respective HDF file.

.. code-block:: bash

    md_davis collate 2VH7_input.toml 2GV1_input.toml 2ACY_input.toml 1URR_input.toml

3. Combine the data from the HDF file into a pandas dataframe with:

.. code-block:: bash

    md_davis residue 2VH7_data.h5 2GV1_data.h5 2ACY_data.h5 1URR_data.h5 -o AcP_residue_data.p

4. Plot the residue properties:

.. code-block:: bash

    md_davis plot_residue AcP_residue_data.p -o AcP_residue_data.html

Now, we can also align the residues of the different trajectories to align the peaks in the data.

1. obtain the sequence of residues in FASTA format from each PDB file using the ``sequence`` command in MD DaVis:

.. code-block:: bash

    md_davis sequence 2VH7\2VH7_structure.pdb -r fasta

2. Use a sequence alignment program or webserver like `Clustal Omega <https://www.ebi.ac.uk/Tools/msa/clustalo/>`_ or `T-coffee <https://www.ebi.ac.uk/Tools/msa/tcoffee/>`_ to obtain the alignment of these sequences in ClustalW format.

.. code-block:: text

    CLUSTAL O(1.2.4) multiple sequence alignment


    2GV1_structure      ---MSKVCIIAWVYGRVQGVGFRYTTQYEAKRLGLTGYAKNLDDGSVEVVACGEEGQVEK	57
    1URR_structure      -VAKQIFALDFEIFGRVQGVFFRKHTSHEAKRLGVRGWCMNTRDGTVKGQLEAPMMNLME	59
    2VH7_structure      ----TLISVDYEIFGKVQGVFFRKHTQAEGKKLGLVGWVQNTDRGTVQGQLQGPISKVRH	56
    2ACY_structure      AEGDTLISVDYEIFGKVQGVFFRKYTQAEGKKLGLVGWVQNTDQGTVQGQLQGPASKVRH	60
                              ..:   ::*:**** **  *. *.*:**: *:  *   *:*:    .   :: .

    2GV1_structure      LMQWLKSGGPRSARVERVLSEPH--HPSGELTDFRIR-	92
    1URR_structure      MKHWLENNRIPNAKVSKAEFSQIQEIEDYTFTSFDIKH	97
    2VH7_structure      MQEWLETRGSPKSHIDKANFNNEKVILKLDYSDFQIVK	94
    2ACY_structure      MQEWLETKGSPKSHIDRASFHNEKVIVKLDYTDFQIVK	98
                        : .**:.    .:::.:.         .   :.* *

3. create a TOML file to specify which alignment file corresponds to which chain and which sequence label corresponds to which data as shown below:

.. code-block:: toml

    [names]
    2GV1 = '2GV1_structure'
    1URR = '1URR_structure'
    2VH7 = '2VH7_structure'
    2ACY = '2ACY_structure'

    [alignment]
    'chain 0' = 'AcP_alingment.clustal_num'

4. Run the ``md_davis residue`` command passing the TOML file with the ``--alignment`` option to generate the pandas dataframes.

.. code-block:: bash

    md_davis residue 2VH7_data.h5 2GV1_data.h5 2ACY_data.h5 1URR_data.h5 --alignment Acp_alignment_input.toml -o AcP_residue_data_aligned.p

5. Plot the aligned data frames.

.. code-block:: bash

    md_davis plot_residue AcP_residue_data_aligned.p -o AcP_residue_data_aligned.html


.. _gmx rms: https://manual.gromacs.org/documentation/current/onlinehelp/gmx-rms.html
.. _gmx gyrate: https://manual.gromacs.org/documentation/current/onlinehelp/gmx-gyrate.html