import argparse
import mdtraj as md
import csv
# My modules
import contacts

## Subclass the Contacts to get Hbonds Class


class Hbond(object):
    """ Defines a hydrogen bond """

    def __init__(self, atom_indices, topology):
        self.donor = Atom(atom_indices[0], topology)
        self.hydrogen = Atom(atom_indices[1], topology)
        self.recipient = Atom(atom_indices[2], topology)
        self.count = 0
        self.time_series = []

    def __repr__(self):
        return f'\n{self.donor} -- {self.hydrogen} -> {self.recipient}: {self.count}'

    def __iter__(self):
        return iter(list(self.donor) + list(self.recipient) + [self.count])


class Hbonds(object):
    """ All hydrogen bonds """

    def __init__(self, topology):
        self.topology = topology
        self.bonds = []

    def __repr__(self):
        return self.bonds

    def __str__(self):
        out_str = 'Seg Ch Res ResSeq Atom -- Seg Ch Res ResSeq Atom --> ' \
            'Seg Ch Res ResSeq Atom : Count\n'
        for bond in self.bonds:
            out_str += bond.__repr__()  # + '\n'
        return out_str

    def __iter__(self):
        return iter([list(_) for _ in self.bonds])

    def parse_indices(self, index_file, group):
        save = False
        with open(index_file) as ndx_file:
            for line in ndx_file:
                if save:
                    indices = [int(_) - 1 for _ in line.strip().split()]
                    self.bonds.append(Hbond(indices, self.topology))
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
            bondwriter.writerow(['Segment', 'Chain', 'Residue', 'ResID', 'Atom',
                                 'Segment', 'Chain', 'Residue', 'ResID', 'Atom', 'Count'])
            for bond in list(self):
                bondwriter.writerow(bond)


def main(args):
    topology = md.load(args.gro).topology
    hydrogen_bonds = Hbonds(topology)
    hydrogen_bonds.parse_indices(index_file=args.ndx, group=args.group)
    hydrogen_bonds.add_counts(xpm_file=args.xpm)
    if args.output:
        hydrogen_bonds.write_csv(args.output)
    print(hydrogen_bonds)


if __name__ == "__main__":
    arguments = contacts.parse_args()
    main(arguments)
