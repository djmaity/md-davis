# -*- coding: utf-8 -*-
"""
Parse contacts evaluated by gmx hbonds

Usage:
  md_davis contacts [options] (--file <.xpm>)
                              (--index <.ndx>)
                              (--structure <.pdb/.gro>)
                              (--group <string>)

  md_davis contacts -h | --help

Options:

  -f, --file <.xpm>             Contact file obtained from GROMACS
  -i, --index <.ndx>            Index file
  -s, --structure <.pdb/.gro>   Structure file
  -g, --group <string>          Group to match from index file to get the list of contacts

  -b, --begin <int>             Frame to start calculation from

  --pickle FILENAME             Save the output to a pickle file
  --hdf FILENAME                Save the output to a HDF file
  --csv FILENAME                Save the output to a CSV file

  -h, --help                    Show this screen.
"""

import argparse
import mdtraj as md
import csv
import re
import pickle
import docopt
import pandas


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

    def __str__(self):
        return f'{self.segment:2} {self.chain_char(self.chain)} ' + \
            f'{self.residue} {self.res_id:3} {self.name:4}'

    def __gt__(self, other):
        for prop1, prop2 in zip(list(self), list(other)):
            if prop1 == prop2:
                continue
            elif prop1 > prop2:
                return True
            else:
                return False
        return False

    def __ge__(self, other):
        if self == other or self > other:
            return True
        return False

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
        return f'\n{repr(self.atom1)} -- {repr(self.atom2)}: {self.count}'

    def __str__(self):
        return f'\n{str(self.atom1)} -- {str(self.atom2)}: {self.count}'

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
        toggle = False
        with open(xpm_file, 'r') as xpm_data:
            for line in xpm_data:
                if line.startswith('/* x-axis:'):
                    toggle = True
                if toggle:
                    if line.startswith('"'):  # Read data
                        contacts = re.search('"(.*)"', line).group(1)
                        bond_seires_array = [
                            1 if char == 'o' else 0 for char in contacts
                        ]
                        self.bonds[bond_idx].time_series = bond_seires_array
                        self.bonds[bond_idx].count = sum(bond_seires_array)
                        bond_idx += 1

    @property
    def to_df(self):
        columns = ['Segment1', 'Chain1', 'Residue1', 'ResSeq1', 'Atom1',
                   'Segment2', 'Chain2', 'Residue2', 'ResSeq2', 'Atom2', 'Count']
        df = pandas.DataFrame(data=list(self), columns=columns)
        return df


def main(argv):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    topology = md.load(args['--structure']).topology
    molecular_contacts = Contacts(topology)
    molecular_contacts.parse_indices(index_file=args['--index'], group=args['--group'])
    molecular_contacts.add_counts(xpm_file=args['--file'])

    print(molecular_contacts)

    df = molecular_contacts.to_df
    if args['--hdf']:
        df.to_hdf(args['--hdf'], key='contacts')
    if args['--csv']:
        df.to_csv(args['--csv'])
    if args['--pickle']:
        df.to_pickle(args['--pickle'])


if __name__ == '__main__':
    main()
