Architecture
============

MD DaVis uses domain-oriented packages. These are the only supported import
locations.

.. code-block:: text

   md_davis/
   |-- core/              # sequence, residue, and collation domain logic
   |-- analysis/          # numerical and scientific analyses
   |   |-- electrostatics/
   |   |-- landscape/
   |   `-- structure/
   |-- io/                # configuration and scientific file formats
   |-- integrations/      # optional external programs such as PyMOL/GROMACS
   |-- visualization/     # reusable plotting and contact visualization
   |-- gui/               # wxPython user interface
   `-- cli/               # argparse command registration and entry point

Dependencies point inward: the CLI and GUI may use visualization and analysis;
analysis may use core and I/O; core must not depend on either user interface.
External tools are isolated behind ``integrations`` so importing the main
package does not require those optional programs to be initialized.

Public imports must use the domain packages directly, for example
``md_davis.core.collation``, ``md_davis.io.xvg``, or
``md_davis.visualization.plotting``.
