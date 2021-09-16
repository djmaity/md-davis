"""md_davis, a tool for comparative analysis of molecular dynamics trajectories

md_davis provides the following commands:

Commands
--------
sequence    display information about the current install
xvg         plot xmgrace (.xvg) files
landscape

collect     list packages linked into a specified environment
create      print information about a specified package
plot        display a list of available conda commands and their help strings
Electrostatics
--------------
calculate   create a new conda environment from a list of specified packages
Landscapes
----------
create      create a free energy landscape

Additional help for each command can be accessed by using::

    md_davis <command> -h

"""

import click
import sys
import md_davis

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=md_davis.__version__, prog_name='md_davis')
def main():
    """MD DaVis: A python package for comparative analysis of molecular
    dynamics trajectories
    """
    pass


main.add_command(md_davis.sequence.main)
main.add_command(md_davis.xvg.main)
main.add_command(md_davis.collate.main)
main.add_command(md_davis.residue.main)
main.add_command(md_davis.landscape.landscape_xvg.landscape_xvg)
main.add_command(md_davis.landscape.landscape_animate.main)
main.add_command(md_davis.hbond.main)
main.add_command(md_davis.plotting.plot_hbond.main)
main.add_command(md_davis.plotting.plot_residue_dataframe.main)
main.add_command(md_davis.electrostatics.electrostatics.main)
main.add_command(md_davis.electrostatics.electrodynamics.main)



#     elif args['<command>'] == 'surface':
#         from .structure import surface
#         surface.main(argv=argv)

#             plot_dipoles.main(argv=argv)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
