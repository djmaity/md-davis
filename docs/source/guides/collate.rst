Collate Analysis Data
=====================

The collate function of MD Davis uses `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ file format to store the heterogeneous data obtained from consolidating and organizing the data from multiple calculations using the h5py Python module. `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ being an open binary format, allows users to open these files directly in C, C++, or FORTRAN programs. It is much smaller to store and faster to load than text files.  Moreover, the data can be inspected with the HDF View GUI.

**Step 1:** Create a TOML file (input.toml) as shown below specifying the
output files from GROMACS.

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
        rmsf = 'rmsf.xvg'
        sasa = 'resarea.xvg'

If the ``chunk`` under ``dihedral`` is provided, MD DaVis will calculate
backbone dihedral angles for all frames and the torsional flexibility
(circular standard deviation) of each dihedral angle.

In case, the RMSF for each chain was calculated separately, the files may be provided with ``rmsf`` key as follows:

.. code-block:: toml

    [residue_property]
        rmsf = ['rmsf_chain_A.xvg', 'rmsf_chain_B.xvg']

**Step 2:** Collect all the data from the output files into a single HDF
file using the following command:

.. code-block:: bash

    md-davis collate input.toml

Multiple input files can also be provided.

.. code-block:: bash

    md-davis collate input1.toml input2.toml

.. note:: The command can be run multiple times to incrementally add data into a HDF5 file by modifying the TOML file. For example, first adding the RMSD and Rg then adding dihedral and so no. However, sometimes if a variable or group already exists in the HDF5 file then h5py will fail to create the file. Please the delete the exisintg file and rerun the command by providing all data at once.






