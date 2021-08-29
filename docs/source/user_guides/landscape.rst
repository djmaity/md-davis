Free Energy Landscapes
======================

MD DaVis can create alchemical free energy landscapes using two variables and
multiple such landscapes can be plotted together. The landscapes can be
plotted using common ranges and same number of bins. The wonderful feature of
the plot created using MD DaVis is that it is an interactive html, and
rotating one subplot rotates all others to the same view. These features
facilitate the quick and efficient comparison of the free energy landscapes.

.. note::
    Plotting free energy landscapes requires a large number of frames,
    therefore, always save enough number of frames during the simulation.
    The number of frames required to obtain a reasonable representation of the
    free energy landscape depends on flexibility of the protein. If the protein
    accesses a large number of conformational states a greater number of
    frames would be required. As a ball park, at least 10000 frames are
    required.

.. warning::
    These commands do not open the output html file in a web browser like
    other plotting commands, so please check the working directory for the
    output html file. If a file with the same name exists, it will be
    overwritten.

Free Energy Landscapes from xvg files
-------------------------------------

The quickest way to plot the free energy landscape is from the .xvg files output by GROMACS.

.. code-block:: bash

    md_davis landscape_xvg -c -T 300 -x rmsd_file_1.xvg -y rg_file_1.xvg -n "name_1" -l "Free Energy Landscape for 1"

Here, the RMSD and R\ :sub:`G`\ precalculated using GROMACS is provided as the ``-x`` and ``-y`` input files. The .xvg files are space or tab delimited text files with time and data in the first and second columns, respectively. The output from any other analysis tool may be used as input once the file is appropriately formatted. Therefore, other properties may be plotted on the x and y axes, for example, the first and the second principal components from a principal component analysis.

If the output filename is not provided the default name is ``landscape.html``.

``-c`` specifies that the ranges of the all the landscapes must be the same

``-T 300`` specifies to use 300 K as the temperature for the boltzmann inversion. If this option is not provided, only the 2D histrogram of the data is plotted.

To plot multiple free energy landscapes as subplots just provide the inputs one after the other. The ``-x``, ``-y``, ``-n``, and ``-l`` for one trajectory must be provided together before the next set of inputs from the other trajectory.

.. code-block:: bash

    md_davis landscape_xvg -c -T 300 -x rmsd_file_1.xvg -y rg_file_1.xvg -n "name_1" -l "Free Energy Landscape for 1" -x rmsd_file_1.xvg -y rg_file_1.xvg -n "name_2" -l "Free Energy Landscape for 2" -x rmsd_file_1.xvg -y rg_file_1.xvg -n "name_3" -l "Free Energy Landscape for 3"

Plot Free Energy Landscapes from HDF Files
------------------------------------------

If the RMSD and R\ :sub:`G`\ have been collated into HDF files for each trajectory. These may be plotted using:

.. code-block:: bash

    md_davis landscape rmsd_rg -c -T 300 file1.h5 file2.h5 file3.h5

Plot Free Energy Landscapes Overlaid with Trajectory Points
-----------------------------------------------------------

One must save the landscape created by ``landscape_xvg`` or ``landscape`` with the ``-s`` option before this one can be used. Since the output generated for single landscape is big, visualization of multiple landscapes becomes impractical. So, it only plots one landscape at a time. Select the desired landscape in ``landscapes.h5`` by providing its index with ``-i``. By default only the first landscape is plotted

.. code-block:: bash

    md_davis landscape animation landscapes.h5 -i 0 --static -o trajectory.html


