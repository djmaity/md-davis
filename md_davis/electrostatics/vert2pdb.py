#! /usr/bin/env python
""" Convet .vert file obtained from MSMS (Michael F. Sanner) to .pdb 
    file to be supplied to Delphi v8.0 as frc file 
    
    This script calculates the closest atom to each vertex and writes 
    that in the PDB output.
"""

import argparse
import numpy
import Bio.PDB
from scipy.spatial import distance_matrix


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

    pdb_parser = Bio.PDB.PDBParser()
    structure = pdb_parser.get_structure("structure", args.pdb)
    model = structure[0]    # Use the first model

    atoms, coordinates = [], []
    for atom in model.get_atoms():
        atoms.append(atom)
        coordinates.append(atom.get_coord())

    coordinates = numpy.array(coordinates)

    distances = distance_matrix(vertices, coordinates)
    atom_indices = numpy.argmin(distances, axis=1)

    for vertex, index in zip(vertices, atom_indices):
        residue = atoms[index].get_parent()
        chain = residue.get_parent()._id
        name = atoms[index].get_fullname()
        resName = residue.get_resname()
        resSeq = residue.id[1]
        x, y, z = vertex
        print(f'ATOM        {name:^4} {resName:3} {chain:1}{resSeq:4}    {x:8.3f}{y:8.3f}{z:8.3f}',
            file=args.output)

if __name__ == "__main__":
    main()