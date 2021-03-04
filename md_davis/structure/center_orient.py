import argparse
import numpy
import pymol
from pymol.cgo import *
from pymol import cmd


def initialize_pymol(window=False):
    """ The commands necessary to get pymol running """
    if window:
        pymol.pymol_argv = ['pymol']
    else:
        pymol.pymol_argv = ['pymol','-cQ']
    pymol.finish_launching()
    return pymol.cmd


def inertia_tensor(selection, name="tensor", state=1):
    """ https://pymolwiki.org/index.php/Inertia_tensor """
    totmass = 0.0
    x_com, y_com, z_com = 0, 0, 0

    model = cmd.get_model(selection, state)

    for a in model.atom:
        if a.symbol != 'e':
            x_com += a.coord[0] * a.get_mass()
            y_com += a.coord[1] * a.get_mass()
            z_com += a.coord[2] * a.get_mass()
            totmass += a.get_mass()

    x_com /= totmass
    y_com /= totmass
    z_com /= totmass

    I = []

    for index in range(9):
        I.append(0)

    for a in model.atom:
        if a.name == 'FE':
            a.symbol = 'Fe'
        temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
        temp_x -= x_com
        temp_y -= y_com
        temp_z -= z_com
        I[0] += a.get_mass() * (temp_y ** 2 + temp_z ** 2)
        I[4] += a.get_mass() * (temp_x ** 2 + temp_z ** 2)
        I[8] += a.get_mass() * (temp_x ** 2 + temp_y ** 2)
        I[1] -= a.get_mass() * temp_x * temp_y
        I[3] -= a.get_mass() * temp_x * temp_y
        I[2] -= a.get_mass() * temp_x * temp_z
        I[6] -= a.get_mass() * temp_x * temp_z
        I[5] -= a.get_mass() * temp_y * temp_z
        I[7] -= a.get_mass() * temp_y * temp_z


    tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])
    vals, vects = numpy.linalg.eig(tensor)  # they come out unsorted, so the command below is needed

    eig_ord = numpy.argsort(vals)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    ord_vals = vals[eig_ord]
    ord_vects = vects[:, eig_ord].T

    return ord_vects, [x_com, y_com, z_com]


def main():
    """ Center PDB and oritent the principal axes to x, y and x directions """
    parser = argparse.ArgumentParser(description='Orient and center a molecule')
    parser.add_argument('pdb_file', help='Input file',
                        metavar='molecule.pdb')
    parser.add_argument('-o', '--output', dest='output',
                        help='Output filename', metavar='output.pdb')
    args = parser.parse_args()

    cmd = initialize_pymol(window=False)
    # cmd = initialize_pymol(window=True)

    cmd.set('retain_order', 1)

    cmd.load(args.pdb_file, object='molecule')
    # cmd.run('axes.py')

    rot_mat, com = inertia_tensor('molecule')

    # import inertia_tensor as it
    # rot_mat = it.tensor('molecule')
    # com = cmd.centerofmass('molecule')

    if numpy.linalg.det(rot_mat) < 0:
        reflect = numpy.matrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
        inverse = numpy.linalg.inv(numpy.array(reflect * rot_mat))
    else:
        inverse = numpy.linalg.inv(numpy.array(rot_mat))

    print(numpy.linalg.det(inverse))

    transform = numpy.zeros((4, 4))
    transform[:3, :3] = inverse
    transform[3, :3] = [-i for i in com]
    transform[3, 3] = 1
    cmd.transform_object('molecule', matrix=transform.flatten())
    if args.output:
        cmd.save(filename=args.output, selection='molecule')


if __name__ == "__main__":
    main()
