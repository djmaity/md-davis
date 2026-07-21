from __future__ import annotations

import argparse
from os import PathLike

import mdtraj
import numpy
from scipy.constants import Boltzmann, R, hbar, nano, u

CONSTANT = (Boltzmann * numpy.exp(1) ** 2 * 298.15 * nano**2 * u) / (hbar**2)


def get_arguments() -> argparse.Namespace:
    """Get filename and arguments from the commandline"""
    parser = argparse.ArgumentParser(description="Calculate configurational entropy using Schlitter's Method")
    parser.add_argument("-f", "--trajectory", dest="trajectory", help="Trajectory file", required=True)
    parser.add_argument("-s", "--structure", dest="structure", help="Structure file", required=True)
    parser.add_argument(
        "-c", "--chunk", dest="chunk", help="Chunk size: Number of frames to load at a time adn calculate its entropy"
    )
    return parser.parse_args()


def calculate_entropy(
    trajectory: str | PathLike[str], structure: str | PathLike[str], chunk_size: int = 1000
) -> None:
    traj_chunks = mdtraj.iterload(trajectory, top=structure, chunk=chunk_size)
    reference = mdtraj.load(structure)

    top = reference.topology
    n_coords = top.n_atoms * 3

    mass = numpy.zeros(n_coords)
    for i, atom in enumerate(top.atoms):
        mass_sqrt = numpy.sqrt(atom.element.mass)
        mass[3 * i] = mass_sqrt
        mass[3 * i + 1] = mass_sqrt
        mass[3 * i + 2] = mass_sqrt

    print("# Entropy in J / (mol K)")
    for chunk in traj_chunks:
        superposed = chunk.superpose(reference)
        coords = superposed.xyz.reshape((chunk_size, n_coords)) * mass
        covariance = numpy.cov(coords.T)
        eigen_values = numpy.linalg.eigvalsh(covariance)
        temp = numpy.ones((1, n_coords)) + CONSTANT * eigen_values
        entropy = 0.5 * R * numpy.sum(numpy.log(temp))
        print(entropy)


if __name__ == "__main__":
    args = get_arguments()
    calculate_entropy(args.trajectory, args.structure)
