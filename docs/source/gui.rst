Graphical User Interface
========================

The GUI (graphical user interface) can be invoked by running the following
command from a terminal:

.. code-block:: bash

   md-davis-gui

Each MD DaVis functionality is in a separate tab of the GUI.
Select the appropriate tab and click the image to see detailed usage.

.. tabs::

   .. tab:: Landscape

      .. tabs::

         .. group-tab:: Linux

            .. image:: /_static/gui/gui-landscape-linux.png
               :target: guides/free_energy_landscapes.html#plot-free-energy-landscapes-using-the-md-davis-gui

         .. group-tab:: Windows

            .. image:: /_static/gui/gui-landscape-windows.png
               :target: guides/free_energy_landscapes.html#plot-free-energy-landscapes-using-the-md-davis-gui

   .. tab:: Electrodynamics

      .. tabs::

         .. group-tab:: Linux

            .. image:: /_static/gui/gui-electrodynamics-linux.png
               :target: guides/electrostatics.html#electric-field-dynamics

         .. group-tab:: Windows

            .. image:: /_static/gui/gui-electrodynamics-windows.png
               :target: guides/electrostatics.html#electric-field-dynamics

   .. tab:: Collate

      .. tabs::

         .. group-tab:: Linux

            .. image:: /_static/gui/gui-collate-linux.png
               :target: guides/collate.html

         .. group-tab:: Windows

            .. image:: /_static/gui/gui-collate-windows.png
               :target: guides/collate.html

   .. tab:: Sequence

      .. tabs::

         .. group-tab:: Linux

            .. image:: /_static/gui/gui-sequence-linux.png
               :target: guides/sequence.html

         .. group-tab:: Windows

            .. image:: /_static/gui/gui-sequence-windows.png
               :target: guides/sequence.html

.. note::

   The GTK+ windowing system in Linux may issue a warning as shown below.
   This does not affect the functionality of MD DaVis and can be safely ignored.

   .. image:: /_static/gui/gui-warning-linux.png

   This may be annoying and bury some the important messages and output from MD DaVis.
   Please pay attention to the terminal used to launch the MD DaVis GUI.
   We will redirect the output to a better interface in a future release.

HDFView
-------

The HDFView (https://www.hdfgroup.org/downloads/hdfview/) GUI program from
`HDF Group <https://www.hdfgroup.org/>`_ can be used to inspect and modify
the HDF files created by MD DaVis. The modified file can be then provided to
MD DaVis plotting commands. See :ref:`guides/hdf:HDF Files` for details.

.. image:: /_static/hdfview/hdfview-windows.png
