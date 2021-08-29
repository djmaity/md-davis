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
main.add_command(md_davis.contacts.main)
main.add_command(md_davis.plotting.plot_residue_dataframe.main)
main.add_command(md_davis.electrostatics.electrostatics.main)


# def main():
#     """Console script entrypoint for md_davis."""
#     args = docopt(__doc__,
#                   version=__version__,
#                   options_first=True)

#     argv = [args['<command>']] + args['<args>']

#     if args['<command>'] == 'collect':
#         from . import collect
#         collect.main(argv=argv)

#     elif args['<command>'] == 'contacts':
#         from .utils import contacts
#         contacts.main(argv=argv)

#     elif args['<command>'] == 'hbonds':
#         from .utils import hbonds
#         hbonds.main(argv=argv)

#     elif args['<command>'] == 'surface':
#         from .structure import surface
#         surface.main(argv=argv)

#     elif args['<command>'] == 'electrostatics':
#         from .electrostatics import surface_electrostatics
#         surface_electrostatics.main(argv=argv)

#     elif args['<command>'] == 'residue':
#         if len(args['<args>']) > 1:
#             if args['<args>'][0] == 'dataframe':
#                 import md_davis.collect_data.create_residue_dataframe as res_df
#                 arguments = docopt(res_df.__doc__, argv=argv)
#                 res_df.main(arguments)
#                 return
#             if args['<args>'][0] == 'aligned':
#                 import md_davis.collect_data.create_aligned_residue_dataframe as res_aligned_df
#                 res_aligned_df.main(argv=argv)
#                 return
#         print('Invalid command. The available commands are:')
#         print('  md_davis residue dataframe')
#         print('  md_davis residue aligned')


#     elif args['<command>'] == 'plot':
#         if len(args['<args>']) < 1:
#             print('Choose a plotting command')
#             return
#         if args['<args>'][0] == 'dipoles':
#             from .plot import plot_dipoles
#             plot_dipoles.main(argv=argv)
#         elif args['<args>'][0] == 'rmsd_rg':
#             from .plot import plot_rmsd_rg
#             plot_rmsd_rg.main()
#         elif args['<args>'][0] == 'residue':
#             import md_davis.plot.plot_residue_dataframe
#             arguments = docopt(
#                 doc = md_davis.plot.plot_residue_dataframe.__doc__,
#                 argv = argv
#             )
#             md_davis.plot.plot_residue_dataframe.main(arguments)
#         else:
#             pass
#     elif args['<command>'] == 'landscape':
#         if len(args['<args>']) > 1:
#             if args['<args>'][0] == 'rmsd_rg':
#                 from .landscape import rmsd_rg_landscape
#                 rmsd_rg_landscape.main(argv=argv)
#                 return
#             if args['<args>'][0] == 'animation':
#                 from .landscape import landscape_animation
#                 landscape_animation.main(argv=argv)
#                 return
#         print('Invalid command. The available commands are:')
#         print('  md_davis landscape rmsd_rg')
#         print('  md_davis landscape animation')
#     else:
#         exit("%r is not a md_davis command. See 'md_davis --help'." % args['<command>'])
#     return


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
