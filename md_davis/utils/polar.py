""" This module converts cartesian coordinates to spherical polar coordinates
    of a structure or trajectory and plots them. """
# TODO: Convert this module to docopt and figure out its utility
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import DSSP
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, atan2, degrees
import argparse
import os


def spherical_numpy(xyz):
    """ Transform cartesian coordinates to spherical polar coordinates """
    ptsnew = np.zeros(xyz.shape)
    r_xy = xyz[:, 0]**2 + xyz[:, 1]**2
    ptsnew[:, 0] = np.sqrt(r_xy + xyz[:, 2]**2)
    ptsnew[:, 1] = np.arctan2(xyz[:, 1], xyz[:, 0])
    # for elevation angle defined from Z-axis down
    ptsnew[:, 2] = np.arctan2(np.sqrt(r_xy), xyz[:, 2])
    # # for elevation angle defined from XY-plane up
    # ptsnew[:, 3] = np.arctan2(xyz[:, 2], np.sqrt(r_xy))
    return ptsnew


def spherical(x, y, z):
    """ Transform cartesian coordinates to spherical polar coordinates """
    xy = x**2 + y**2
    r = sqrt(xy + z**2)
    theta = degrees(atan2(x, y))
    # for elevation angle defined from Z-axis down
    phi = degrees(atan2(sqrt(xy), z))
    # # for elevation angle defined from XY-plane up
    # phi = degrees(atan2(z, sqrt(xy)))
    return r, theta, phi


def main():
    """ main function """
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('pdb', help='Input PDB file')
    parser.add_argument('-o', '--output', help='Output directory')
    args = parser.parse_args()

    basename = os.path.splitext(os.path.basename(args.pdb))[0]

    parser = PDBParser()
    structure = parser.get_structure('structure', args.pdb)

    atom_output = []
    residue_output = []
    c_alpha_output = []
    model = structure[0]

    dssp = DSSP(model, args.pdb, dssp='mkdssp')

    for chain in model:
        for residue in chain:
            residue_sph_coord = []
            c_alpha_output.append(spherical(*residue['CA'].get_coord()))
            for atom in residue:
                sph_coord = spherical(*atom.get_coord())
                atom_output.append(sph_coord)
                residue_sph_coord.append(sph_coord)
            residue_output.append(np.mean(residue_sph_coord, axis=0))

    atom_output = np.array(atom_output).T
    residue_output = np.array(residue_output).T
    c_alpha_output = np.array(c_alpha_output).T

    if args.output:
        output_file = os.path.join(args.output, basename)
    else:
        output_file = basename
    np.savez(output_file, atom_wise=atom_output, residue_wise=residue_output,
             c_alpha=c_alpha_output, dssp=dssp)


if __name__ == '__main__':
    main()
