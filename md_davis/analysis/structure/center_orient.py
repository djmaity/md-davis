from __future__ import annotations

import argparse
from typing import Any

import numpy
from numpy.typing import NDArray

from md_davis.integrations.pymol import initialize_pymol


def inertia_tensor(
    selection: str, name: str = "tensor", state: int = 1, command: Any | None = None
) -> tuple[NDArray, list[float]]:
    """https://pymolwiki.org/index.php/Inertia_tensor"""
    if command is None:
        command = initialize_pymol(window=False)
    totmass = 0.0
    x_com, y_com, z_com = 0, 0, 0

    model = command.get_model(selection, state)

    for a in model.atom:
        if a.symbol != "e":
            x_com += a.coord[0] * a.get_mass()
            y_com += a.coord[1] * a.get_mass()
            z_com += a.coord[2] * a.get_mass()
            totmass += a.get_mass()

    x_com /= totmass
    y_com /= totmass
    z_com /= totmass

    inertia = numpy.zeros(9)

    for a in model.atom:
        if a.name == "FE":
            a.symbol = "Fe"
        temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
        temp_x -= x_com
        temp_y -= y_com
        temp_z -= z_com
        inertia[0] += a.get_mass() * (temp_y**2 + temp_z**2)
        inertia[4] += a.get_mass() * (temp_x**2 + temp_z**2)
        inertia[8] += a.get_mass() * (temp_x**2 + temp_y**2)
        inertia[1] -= a.get_mass() * temp_x * temp_y
        inertia[3] -= a.get_mass() * temp_x * temp_y
        inertia[2] -= a.get_mass() * temp_x * temp_z
        inertia[6] -= a.get_mass() * temp_x * temp_z
        inertia[5] -= a.get_mass() * temp_y * temp_z
        inertia[7] -= a.get_mass() * temp_y * temp_z

    tensor = inertia.reshape((3, 3))
    _, eigenvectors = numpy.linalg.eigh(tensor)
    return eigenvectors.T, [x_com, y_com, z_com]


def main() -> None:
    """Center PDB and oritent the principal axes to x, y and x directions"""
    parser = argparse.ArgumentParser(description="Orient and center a molecule")
    parser.add_argument("pdb_file", help="Input file", metavar="molecule.pdb")
    parser.add_argument("-o", "--output", dest="output", help="Output filename", metavar="output.pdb")
    args = parser.parse_args()

    cmd = initialize_pymol(window=False)
    cmd.set("retain_order", 1)

    cmd.load(args.pdb_file, object="molecule")
    rot_mat, com = inertia_tensor("molecule", command=cmd)

    if numpy.linalg.det(rot_mat) < 0:
        reflect = numpy.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        inverse = numpy.linalg.inv(reflect @ rot_mat)
    else:
        inverse = numpy.linalg.inv(rot_mat)

    transform = numpy.zeros((4, 4))
    transform[:3, :3] = inverse
    transform[3, :3] = [-i for i in com]
    transform[3, 3] = 1
    cmd.transform_object("molecule", matrix=transform.flatten())
    if args.output:
        cmd.save(filename=args.output, selection="molecule")


if __name__ == "__main__":
    main()
