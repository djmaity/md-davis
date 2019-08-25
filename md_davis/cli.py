"""
MD DaVis: A python package to analyze molecular dynamics trajecotires
          of proteins

Usage:
  md_davis <command> [<args>...]
  md_davis -h | --help
  md_davis --version

Options:
  -h --help     Show this screen.
  --version     Show version.

Available Commands:
  sequence      get the sequence from a PDB file
  plot_xvg      Plot an .xvg file created by GROMACS
  collect       Collect data from various calculations into one HDF binary file
  plot_dipoles          Plot commands
"""

from docopt import docopt
import subprocess

def main():
    args = docopt(__doc__,
                  version='0.0.2',
                  options_first=True)
    # print('global arguments:')
    # print(args)
    # print('command arguments:')

    argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'sequence':
        import md_davis.structure.sequence as seq
        arguments = docopt(seq.__doc__, argv=argv)
        seq.main(arguments)
    elif args['<command>'] == 'collect':
        import md_davis.collect
        arguments = docopt(md_davis.collect.__doc__, argv=argv)
        md_davis.collect.main(arguments)
    elif args['<command>'] == 'plot_xvg':
        import md_davis.plot_xvg
        arguments = docopt(md_davis.plot_xvg.__doc__, argv=argv)
        md_davis.plot_xvg.main(arguments)
    elif args['<command>'] == 'plot_dipoles':
        import md_davis.plotly_plots.plot_dipoles
        arguments = docopt(md_davis.plotly_plots.plot_dipoles.__doc__, argv=argv)
        md_davis.plotly_plots.plot_dipoles.main(arguments)
    else:
        pass


if __name__ == "__main__":
    main()
    