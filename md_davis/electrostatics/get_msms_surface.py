#! /usr/bin/env python

""" This module finds the electrostatics from a 3D potential map in Gaussian
    (Software) cube format at the points on the surface of protein. The points
    on the surface of the protein are obtained by MSMS program from MGL Tools

    Author: Dibyajyoti Maity
"""

import os
import argparse
import subprocess


def get_arguments():
    """ Get filename and arguments from the commandline """
    parser = argparse.ArgumentParser(description='Calculate triangluated surface using MSMS program')
    parser.add_argument('pdb',  metavar='filename.pdb',
                        help='PDB structure file')
    parser.add_argument('-o', '--output', dest='output', help='Output directory')
    return parser.parse_args()


def run_msms(pdb_file, output_directory=None, msms_path=None):
    basename, _ = os.path.splitext(os.path.basename(pdb_file))

    if msms_path:
        if msms_path[-1] != '/':
            msms_path = msms_path + '/'

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

    subprocess.run([msms_path + "msms.x86_64Linux2.2.6.1",
                     "-if", output_directory + basename + ".xyz",
                     "-of", output_directory + basename,
                     "-probe_radius", "1.4"])
    return output_directory + basename + '.vert'


def main():
    args = get_arguments()
    run_msms(args.pdb, output_directory=args.output, msms_path='/home/djmaity/.opt/msms')


if __name__ == '__main__':
    main()
