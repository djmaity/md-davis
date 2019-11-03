import argparse
import mdtraj as md
import csv


class Atom(object):
    """ Wrapper for MDTraj topology.atom object """

    def __init__(self, index, topology):
        atom = topology.atom(index)
        self.segment = atom.residue.segment_id
        self.chain = atom.residue.chain.index
        self.residue = atom.residue.name
        self.res_id = atom.residue.resSeq
        self.name = atom.name

    def __repr__(self):
        return f'{self.segment:2} {self.chain_char(self.chain)} ' + \
            f'{self.residue} {self.res_id:3} {self.name:4}'

    @staticmethod
    def chain_char(index):
        return chr(index + 65)

    def __iter__(self):
        """ Segment Chain Residue ResSeq Atom """
        return iter([self.segment, self.chain_char(self.chain), self.residue,
                     self.res_id, self.name])


class Contact(object):
    """ Defines a contact between two atoms """

    def __init__(self, atom_indices, topology):
        self.atom1 = Atom(atom_indices[0], topology)
        self.atom2 = Atom(atom_indices[-1], topology)
        self.count = 0
        self.time_series = []

    def __repr__(self):
        return f'\n{self.atom1} -- {self.atom2}: {self.count}'

    def __iter__(self):
        return iter(list(self.atom1) + list(self.atom2) + [self.count])


class Contacts(object):
    """ All hydrogen bonds """

    def __init__(self, topology):
        self.topology = topology
        self.bonds = []

    def __repr__(self):
        return self.bonds

    def __str__(self):
        out_str = 'Segment Chain Residue ResSeq Atom -- Segment Chain Residue ResSeq Atom : Count\n'
        for bond in self.bonds:
            out_str += bond.__repr__()  # + '\n'
        return out_str

    def __iter__(self):
        return iter([list(_) for _ in self.bonds])

    def __len__(self):
        return len(self.bonds)

    def parse_indices(self, index_file, group):
        save = False
        with open(index_file) as ndx_file:
            for line in ndx_file:
                if save:
                    indices = [int(_) - 1 for _ in line.strip().split()]
                    self.bonds.append(Contact(indices, self.topology))
                if line == f'[ {group} ]\n':
                    save = True

    def add_counts(self, xpm_file):
        bond_idx = 0
        with open(xpm_file, 'r') as xpm_data:
            for line in xpm_data:
                if line.startswith('"') and len(line) > 1000:  # Read data
                    bond_seires_array = [1 if char ==
                                         'o' else 0 for char in line]
                    self.bonds[bond_idx].time_series = bond_seires_array
                    self.bonds[bond_idx].count = sum(bond_seires_array)
                    bond_idx += 1

    def write_csv(self, csv_file):
        with open(csv_file, 'w') as csvfile:
            bondwriter = csv.writer(
                csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
            bondwriter.writerow(['Segment1', 'Chain1', 'Residue1', 'ResSeq1', 'Atom1',
                                 'Segment2', 'Chain2', 'Residue2', 'ResSeq2', 'Atom2', 'Count'])
            for bond in list(self):
                bondwriter.writerow(bond)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('xpm', help='xpm file')
    parser.add_argument('ndx', help='ndx file')
    parser.add_argument('gro', help='gro file')
    parser.add_argument(
        'group', help='group to match from index file to get the list of h-bonds')
    parser.add_argument('-o', '--out', dest='output', help='Save Hydrogen bonds '
                        'data in csv')
    parser.add_argument('-p', '--prefix', dest='prefix', help='Prefix of '
                        'output filenames')
    parser.add_argument('-b', '--begin', dest='begin', type=int,
                        help='Frame to start calculation from')
    return parser.parse_args()


def main(args):
    topology = md.load(args.gro).topology
    hydrogen_bonds = Contacts(topology)
    hydrogen_bonds.parse_indices(index_file=args.ndx, group=args.group)
    hydrogen_bonds.add_counts(xpm_file=args.xpm)
    if args.output:
        hydrogen_bonds.write_csv(args.output)
    print(hydrogen_bonds)


if __name__ == "__main__":
    arguments = parse_args()
    main(arguments)
