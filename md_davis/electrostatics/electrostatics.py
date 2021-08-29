""" This module finds the electrostatics from a 3D potential map in Gaussian
    (Software) cube format at the points on the surface of protein. The points
    on the surface of the protein are obtained by MSMS program from MGL Tools

    Author: Dibyajyoti Maity

Usage:  md_davis electrostatics [options] PDB_FILE OUTPUT_DIRECTORY

Options:
  -m, --msms PATH                   full path to MSMS executable [default: ~/]
  -d, --delphi PATH                 full path to Delphi executable [default: ~/]
  -r, --radius PATH                 path to radius file
  -c, --charge PATH                 path to charge file
  -g, --grid_size <odd_int>         Grid size to use for Delphi calculation
  -v, --vertices <filename.vert>    vertices file generated by MSMS
  -s, --surface <surface.pdb>       file containing surface vertices created by
                                    running 'md_davis vert2pdb' on vertex file
                                    created by MSMS
  --surface_potential               Whether to calculate the electrostatic
                                    potential on the surface or not
  --center                          Center the grid for Delphi at the origin
"""

import os
import click
import subprocess
import pandas


def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise NotADirectoryError(path)


# TODO: Remove hard coded application paths and names
def run_msms(pdb_file, output_directory=None, msms_path=None, executable='msms'):
    """ Calculate triangulated surface using MSMS program """
    basename, _ = os.path.splitext(os.path.basename(pdb_file))

    if not output_directory:
        output_directory = './'
    elif output_directory[-1] != '/':
        output_directory = output_directory + '/'
    else:
        pass

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    with open(output_directory + basename + '.xyz', 'w') as xyz_file:
        subprocess.run([msms_path + "pdb_to_xyzrn", pdb_file], stdout=xyz_file)

    subprocess.run([msms_path + executable,
                     "-if", basename + ".xyz",
                     "-of", basename,
                     "-probe_radius", "1.4"], cwd=output_directory)
    return output_directory + basename + '.vert'


def run_delphi(pdb_file, output_directory, output_filename,
    delphi_path, radius_file, charge_file, grid_size=101, surface=None, center=False):
    """ Run Delphi on protein surface created by MSMS program """
    # TODO: Rewrite using template string
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    parameters = [
        f'in(pdb,file="{pdb_file}")',
        f'in(siz,file="{radius_file}")',
        f'in(crg,file="{charge_file}")',
        f'gsize={grid_size}',
        f'salt=0.10',
        f'exdi=80',
        f'linit=2000',
        f'maxc=0.0000000001',
        f'out(phi, file="{output_directory}/{output_filename}.cub", format="cube")',
    ]
    if center:
        parameters += ['acenter(0,0,0)']
    if surface:
        parameters += [
            f'in(frc,file="{surface}")',
            f'out(frc, file="{output_directory}/{output_filename}.pot")',
            f'site(Atom, Potential, Reaction, Coulomb, Field)',
        ]
    print('\n'.join(parameters) + '\n', file=open(f'{output_filename}_tmp.prm', 'w'))
    subprocess.run([delphi_path, f'{output_filename}_tmp.prm'])
    subprocess.run(['rm', f'{output_filename}_tmp.prm'])


def parse_electrostatic_potential(potential_file):
    df = pandas.read_fwf(potential_file, skiprows=12, skipfooter=2,
                         dtype={'resSeq': int}, engine='python',
                         names=['name', 'resName', 'chainID', 'resSeq',
                                'potential', 'reaction', 'coulomb',
                                'Ex', 'Ey', 'Ez'
                         ],
                         widths=[5, 3, 3, 9, 10, 10, 10, 10, 10, 10]
    )
    df['chainID'].fillna('A', inplace=True)
    output = {}
    chain = 0
    for _, data in df.groupby(['chainID'], as_index=False):
        grouped_df = data.groupby(['resSeq'], as_index=False)['potential']
        potential = grouped_df.sum()
        potential.rename(columns={'potential':'total'}, inplace=True)
        potential['mean'] = grouped_df.mean()['potential']
        output[f'chain {chain}'] = potential
        chain += 1
    return output


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='electrostatics', context_settings=CONTEXT_SETTINGS)
@click.option('-m', '--msms', metavar='PATH', type=click.Path(exists=True),
              default='./', help='full path to MSMS executable')
@click.option('-d', '--delphi', metavar='PATH', type=click.Path(exists=True),
              default='./', help='full path to Delphi executable')
@click.option('-r', '--radius', metavar='PATH', type=click.Path(exists=True),
              help='path to radius file')
@click.option('-c', '--charge', metavar='PATH', type=click.Path(exists=True),
              help='path to charge file')
@click.option('-g', '--grid_size', metavar='<odd_int>',
              help='Grid size to use for Delphi calculation')
@click.option('-v', '--vertices', metavar='<filename.vert>',
              help='vertices file generated by MSMS')
@click.option('-s', '--surface', metavar='<surface.pdb>',
              help='file containing surface vertices created by running '
              "'md_davis vert2pdb' on vertex file created by MSMS")
@click.option('--surface_potential', help='Whether to calculate the '
              'electrostatic potential on the surface or not')
@click.option('--center', help='Center the grid for Delphi at the origin')
@click.argument('pdb_file')
@click.argument('outdir')
def main(pdb_file, outdir, msms, delphi, radius, charge, grid_size, vertices, surface, surface_potential, center):
    """ Get the electrosatic potential on the surface points generated by MSMS """

    this_script_path = os.path.dirname(os.path.realpath(__file__))

    pdb_file = pdb_file
    output_filename = os.path.splitext(os.path.basename(pdb_file))[0]

    surface_file = None
    if surface_potential:
        if surface:
            surface_file = surface
        else:
            if vertices:
                vert_file = vertices
            else:
                vert_file = run_msms(pdb_file=pdb_file,
                                     output_directory=outdir,
                                     msms_path=msms)
            surface_file = f"{outdir}/{output_filename}_surf.pdb"
            # TODO: remove this subprocess run with a function
            subprocess.run(['python', f'{this_script_path}/vert2pdb.py', vert_file, pdb_file, '-o', surface_file])

    # run_delphi(pdb_file=pdb_file,
    #            output_directory=outdir,
    #            output_filename=output_filename,
    #            surface=surface_file,
    #            delphi_path=delphi,
    #            radius_file=radius,
    #            charge_file=charge,
    #            grid_size=grid_size,
    #            center=center,
    # )


if __name__ == '__main__':
    main()
