Contact/H-bond matrix
=====================

Motivation
----------

GROMACS has an excellent tool to calculate hydrogen bonds (H-bonds) or contacts from trajectories called `gmx hbond`_. It is written in C++, and the loops are parallelized. Therefore, its performance cannot be matched by any Python code. The default output from `gmx hbond`_ is a ``.xvg`` file containing the total number of H-bonds/contacts between the selected groups over all frames, which is not very informative. The option ``-hbm`` outputs the existence matrix for all H-bonds/contacts over all frames in an `X PixMap`_ (``.xpm``) image:

.. image:: /_static/gmx_hbond_xpm.png

Only a few image viewers support ``.xpm`` file format, and the image cannot be easily interpreted because it does not contain any labels for the participating atoms and residues in the H-bonds/contacts. To overcome these limitations, use MD DaVis for analyzing the output files created by ``gmx hbond``.

Steps
-----
1. Calculate the H-bonds/contacts using ``gmx hbond`` with the options  ``-hbm`` and ``-hbn``:

.. code-block:: bash

    gmx hbond -f trajectory.xtc -s structure.tpr -num hbnum.xvg -hbm hb_matrix -hbn hb_index

To calculate the contacts provide ``--contact`` option. ``gmx hbond`` requires the ``.tpr`` file to be supplied to the ``-s`` option, and other structure formats like ``.pdb`` or ``.gro`` are not accepted. The order of rows in the hb_matrix.xpm image is the same as the order of atoms in the hb_index.ndx file. MD DaVis utilizes this to determine the atoms in each H-bond/contact (rows in the ``.xpm`` file).

2. Create and save the H-bonds/contacts data:

.. code-block:: bash

    md_davis contacts -f hb_matrix.xpm --index hb_index.ndx -s structure.tpr --pickle hb_data.p -g hbonds_Protein

3. Plot the H-bonds/contacts matrix:


Bonus Tip
---------

Scroll down the hb_index.ndx file to find the name of the group containing the atoms participating in H-bonds/contacts. It contains three or two columns of data for H-bonds or contacts, respectively.

.. code::

     1235  1243  1244  1248  1249  1260  1261  1263  1264  1280  1285  1286  1301  1302  1320
     1321  1339  1340  1356  1361  1362  1380  1381  1389  1390  1392  1393  1410  1413  1414
     1421  1424  1425  1433  1434  1436  1437  1456  1457  1468  1469  1473  1474  1492  1493
     1508  1509  1525  1530  1531
    [ hbonds_Protein ]
          9     10     35
         36     37    773
         55     56   1285
         62     63    725
         62     63    726

.. _gmx hbond: https://manual.gromacs.org/documentation/current/onlinehelp/gmx-hbond.html
.. _X PixMap: https://en.wikipedia.org/wiki/X_PixMap


