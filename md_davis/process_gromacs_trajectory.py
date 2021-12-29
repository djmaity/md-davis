import subprocess
import re
import pathlib
import click


def first_frame(trajectory, gmx_path='gmx'):
    """Get the time of the first frame of a GROMACS trajectory file"""
    process = subprocess.run([gmx_path, 'dump', '-f', trajectory],
                             capture_output=True, encoding='utf-8')
    m = re.search('time=([\d\.e\+\-]+)', process.stdout)
    time = float(m.group(1))
    return time


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(name='process', context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--index', type=click.Path(exists=True),
              help='Index file for processing the system.')
@click.option(
    '-o', '--output', type=click.Path(), default='./', help='Output directory')
@click.option(
    '-s', '--selection', default='non-Water_&_!Ion',
    help='Name of the group defining the system atoms in the index file.')
@click.option('-g', '--gmx_path', default='non-Water_&_!Ion',
              help='If gmx command is not available in the terminal. You may provide the full path to gmx executable')
@click.argument('trajectory', required=True, type=click.Path(exists=True))
@click.argument('tpr', required=True, type=click.Path(exists=True))
def process_gromacs_trajectory(trajectory, tpr, output='./',
                               index=None, prefix='',
                               selection='non-Water_&_!Ion',
                               gmx_path='gmx'):
    """Process a GROMACS trajecotry to correct for periodic boundary condition.
    """

    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    if index is None:
        index = f'{output}/{prefix}index.ndx'
        if selection == 'non-Water_&_!Ion':
            subprocess.run([gmx_path, 'make_ndx', '-f', tpr, '-o', index],
                           input='"non-Water"&!"Ion"\nq\n',
                           encoding='ascii')
        if selection == 'Protein':
            subprocess.run([gmx_path, 'make_ndx', '-f', tpr, '-o', index],
                           input='"Protein"\nq\n',
                           encoding='ascii')

    begin = str(first_frame(trajectory, gmx_path=gmx_path))

    subprocess.run([gmx_path, 'trjconv',
                    '-f', trajectory,
                    '-s', tpr,
                    '-o', f'{output}/{prefix}whole.pdb',
                    '-n', index,
                    '-pbc', 'whole',
                    '-dump', begin],
                   input=selection, encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', trajectory,
                    '-s', tpr,
                    '-o', f'{output}/{prefix}whole.xtc',
                    '-n', index,
                    '-pbc', 'whole'],
                   input=selection, encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}whole.xtc',
                    '-s', f'{output}/{prefix}whole.pdb',
                    '-o', f'{output}/{prefix}nojump.pdb',
                    '-pbc', 'nojump',
                    '-dump', begin],
                   input='System', encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}whole.xtc',
                    '-s', f'{output}/{prefix}whole.pdb',
                    '-o', f'{output}/{prefix}nojump.xtc',
                    '-pbc', 'nojump'],
                   input='System', encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}nojump.xtc',
                    '-s', tpr,
                    '-o', f'{output}/{prefix}center_mol.pdb',
                    '-n', index, '-center',
                    '-pbc', 'mol', '-ur', 'compact',
                    '-dump', begin],
                   input=f'{selection}\n{selection}', encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}nojump.xtc',
                    '-s', tpr,
                    '-o', f'{output}/{prefix}center_mol.xtc',
                    '-n', index, '-center',
                    '-pbc', 'mol', '-ur', 'compact'],
                   input=f'{selection}\n{selection}', encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}center_mol.xtc',
                    '-s', f'{output}/{prefix}center_mol.pdb',
                    '-o', f'{output}/{prefix}processed.pdb',
                    '-fit', 'progressive',
                    '-dump', begin],
                   input='System  System', encoding='ascii')

    subprocess.run([gmx_path, 'trjconv',
                    '-f', f'{output}/{prefix}center_mol.xtc',
                    '-s', f'{output}/{prefix}center_mol.pdb',
                    '-o', f'{output}/{prefix}processed.xtc',
                    '-fit', 'progressive'],
                   input='System  System', encoding='ascii')


if __name__ == '__main__':
    process_gromacs_trajectory()
