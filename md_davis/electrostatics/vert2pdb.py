#! /usr/bin/env python
""" Convet .vert file obtained from MSMS (Michael F. Sanner) to .pdb 
    file to be supplied to Delphi v8.0 as frc file 
    
    This script calculates the closest atom to each vertex and writes 
    that in the PDB output.
"""

import argparse
import numpy
from sklearn.neighbors import KDTree
from biopandas.pdb import PandasPdb


def main():
    parser = argparse.ArgumentParser(
        description='Convert triangluated surface from MSMS program to PDB for input to Delphi')
    parser.add_argument('vert',  metavar='filename.vert',
                        help='MSMS output vertices file')
    parser.add_argument('pdb',  metavar='filename.pdb',
                        help='PDB file used for vertex calculation')
    parser.add_argument('-o', '--output',  metavar='surface.pdb', default=None, type=argparse.FileType('w'),
                        help='Output PDB file with the surface')
    args = parser.parse_args()

    # Some Input files for hemoglobin simulation were raising UnicodeDecodeError
    vertex_file = open(args.vert, 'r', errors='replace')
    vertices = numpy.loadtxt(vertex_file, usecols=(0, 1, 2), skiprows=3)

    pdb_file = PandasPdb().read_pdb(args.pdb)
    df = pdb_file.df['ATOM']
    coordinates = df[['x_coord', 'y_coord', 'z_coord']]

    kdt = KDTree(coordinates, metric='euclidean')
    atom_indices = kdt.query(vertices, return_distance=False, dualtree=True)

    for vertex, index in zip(vertices, atom_indices):
        residue = df.at[int(index), 'residue_name']
        chain = df.at[int(index), 'chain_id']
        name = df.at[int(index), 'atom_name']
        resSeq = df.at[int(index), 'residue_number']
        x, y, z = vertex
        print(f'ATOM        {name:^4} {residue:3} {chain:1}{resSeq:4}    {x:8.3f}{y:8.3f}{z:8.3f}',
            file=args.output)

if __name__ == "__main__":
    main()