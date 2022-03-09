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

It can also use an alignment to align the residue level data from different proteins along the x-axis. This ensures that the peaks line up properly for better interpretation.

Traditionally, each time the analysis of MD trajectories are compared new statis plots have to be created.

Conventional analysis requires plotting the data each time new properties are compared. This repetitive process is alleviated by making overlaid plots of all the informative residue propertie

.. note:: the paths in the input toml file is relative to the location where the md-davis command will be called from. To avoid any confusion try using absolute paths.

Note that each data is optional in the residue property
plot. Therefore, a plot can be created even if the secondary structure is not
available, albeit devoid of rich structural information.

Secondary Structure
-------------------

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

How to interact with the plot
-----------------------------

## Step 3: Plotting overlaid residue data
**Step 3a:** Create a pickle file with the residue dataframe using:

.. code-block:: bash

    md-davis residue dataframe --prefix name1 output1.h5 data1.p

The optional argument ``-a annotations.json`` can be provided to place a mark at certain residue locations. The contents of ``annotations.json`` should be of the following form:

.. code-block:: toml

    {
        "chain 0": {"Active Site": [23, 41], "Substrate Binding Site": [56]},
        "chain 1": {"Nucleotide Binding Regions": [15, 18]}
    }

Each type of annotation is rendered with a different mark. Following annotations are available at present:
* Active Site
* Nucleotide Binding Regions
* NADP Binding Site
* Substrate Binding Site
* Metal Binding Site
* Cofactor Binding Site
* Mutation

**Step 3b:** If your proteins are of different lengths and you need the peaks to be aligned, create a JSON file as shown below.

.. code-block:: toml

    {
        "alignment": "path/to/alignment_file.clustal_num",
        "locations": {
            "name1": "name1_residue_wise_data.p",
            "name2": "name2_residue_wise_data.p",
            "name3": "name3_residue_wise_data.p"
        },
        "output": "acylphosphatase_residue_wise_data_aligned.p"
    }


The contents of the alignment file, ``alignment_file.clustal_num`` must be
in CLUSTAL format.


**Step 3b:** Plot the residue data pickle file from the previous command using:

.. code-block:: bash

    md-davis plot residue data1.p data2.p


..
    Annotations
    -----------

    .. code-block:: toml

        {
            "chain 0": {"Active Site": [23, 41], "Substrate Binding Site": [56]},
            "chain 1": {"Nucleotide Binding Regions": [15, 18]}
        }

    Each type of annotation is rendered with a different mark. Following annotations are available at present:
    * Active Site
    * Nucleotide Binding Regions
    * NADP Binding Site
    * Substrate Binding Site
    * Metal Binding Site
    * Cofactor Binding Site
    * Mutation


Note the numbers at the end of the ``--rmsf`` options are the start and end
time for the RMSF calculation in nanosecond. These will be inserted as
attributes in the HDF file and must be provided. In case, the RMSF for each
chain was calculated separately, the files may be provided to ``--rmsf``
option in the correct order followed by the start and end times.

The optional argument ``-a annotations.json`` can be provided to place a mark at certain residue locations. The contents of ``annotations.json`` should be of the following form:

.. code-block:: toml

    {
        "chain 0": {"Active Site": [23, 41], "Substrate Binding Site": [56]},
        "chain 1": {"Nucleotide Binding Regions": [15, 18]}
    }

Each type of annotation is rendered with a different mark. Following annotations are available at present:
* Active Site
* Nucleotide Binding Regions
* NADP Binding Site
* Substrate Binding Site
* Metal Binding Site
* Cofactor Binding Site
* Mutation

**Step 3b:** If your proteins are of different lengths and you need the peaks to be aligned, create a JSON file as shown below.

.. code-block:: toml

    {
        "alignment": "path/to/alignment_file.clustal_num",
        "locations": {
            "name1": "name1_residue_wise_data.p",
            "name2": "name2_residue_wise_data.p",
            "name3": "name3_residue_wise_data.p"
        },
        "output": "acylphosphatase_residue_wise_data_aligned.p"
    }

The contents of the alignment file, ``alignment_file.clustal_num`` must be in CLUSTAL format; for example::

    CLUSTAL O(1.2.4) multiple sequence alignment

    name1      --STARPLKSVDYEVFGRVQGVCFRMYAEDEARKIGVVGWVKNTSKGTVTGQVQGPEEKV	58
    name2      --------PRLVALVKGRVQGVGYRAFAQKKALELGLSGYAENLPDGRVEVVAEGPKEAL	52
    name3      ---VAKQIFALDFEIFGRVQGVFFRKHTSHEAKRLGVRGWCMNTRDGTVKGQLEAPMMNL	57
                            : *:**** :*  .  :. .  : *:  *   * *     .    :

    name1      NSMKSWLSKVGSPSSRIDRTNFSNEKTISKLEYSNFSVRY	98
    name2      ELFLHHLKQ--GPRLARVEAVEVQWGEE--AGLKGFHVY-	87
    name3      MEMKHWLENNRIPNAKVSKAEFSQIQEIEDYTFTSFDIKH	97
                :   :     *          :           * :

