Introduction
============

Molecular dynamics (MD) simulations are indispensable for gaining atomistic insights into the biomolecular function.
MD of increasingly larger proteins has become accessible to researchers with recent advancements in computing and MD algorithms.
However, the analysis of MD trajectories still remains tedious.
MD DaVis (Molecular Dynamics Data Visualizer) is a tool and Python 3 package to perform comparative data analysis of MD trajectories of similar proteins or the same protein under different conditions.

There are many :ref:`MD analysis tools <Molecular Dynamics Analysis Tools>`.
However, the output from most has to be visualized using another plotting library requiring a significant amount of coding.
The MD DaVis package provides a command-line tool to create helpful interactive visualizations easily.
The tool can increase the productivity of researchers simulating protein using MD and make the analysis of such simulations accessible to everyone.

Features of MD DaVis
--------------------

.. panels::
   :body: text-center

   .. link-button:: guides/free_energy_landscapes
      :type: ref
      :text: Free Energy Landscapes
      :classes: stretched-link font-weight-bold

   .. image:: /_static/free_energy_landscape.png
      :alt: Free energy landscape

   ---

   .. link-button:: guides/residue_property_plot
      :type: ref
      :text: Residue Properties Plot
      :classes: stretched-link font-weight-bold

   .. image:: /_static/residue_property_plot.png
      :alt: Residue property plot

   ---

   .. link-button:: Surface Electrostatic Potential Per Residue
      :type: ref
      :text: Surface Electrostatics
      :classes: stretched-link font-weight-bold

   .. image:: /_static/surface_electrostatics.png
      :alt: Surface electrostatics

   ---

   .. link-button:: Electric Field Dynamics
      :type: ref
      :text: Electric Field Dynamics
      :classes: stretched-link font-weight-bold

   .. image:: /_static/2VH7_electrodynamics.webp
      :alt: Electric field dynamics

   ---

   .. link-button:: guides/hbond_or_contact_matrix
      :type: ref
      :text: H-bond / Contact Matrix
      :classes: stretched-link font-weight-bold

   .. image:: /_static/2VH7_hbond_matrix.png
      :alt: H-bond matrix

Molecular Dynamics Analysis Tools
---------------------------------

MD DaVis does not implement analysis tools and functions available elsewhere.
Presently, MD DaVis natively supports analysis performed using GROMACS tools.
However, other analysis tools may be used as well, depending on the format of the MD trajectory.
In such cases, the output may be formatted to that accepted by MD DaVis.
The following is a non-exhaustive list of tools for analysis of MD simulations:

* MDTraj (McGibbon et al., 2015)
* MDAnalysis (Michaud-Agrawal et al., 2011)
* MD-TASK (Brown et al., 2017)
* TRAVIS (Brehm and Kirch-ner, 2011)
* MDTRA (Popov et al., 2013)
* Bio3d (Grant et al., 2006)
* MDplot (Margreitter and Oostenbrink, 2017)

* GROMACS (Berendsen et al., 1995)
* AMBER (Case et al., 2020)

