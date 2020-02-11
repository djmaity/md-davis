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
import csv
import re
import pickle
import docopt
import pandas
from biopandas.pdb import PandasPdb
# Local imports
import .contacts


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

# Subclass the Contacts to get Hbonds Class
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

def main(argv):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    structure = PandasPdb().read_pdb(args['--structure'])
    molecular_contacts = Hbonds(structure)
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

if __name__ == "__main__":
    main()
