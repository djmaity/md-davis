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

  -f, --file <.xpm>             Contact file obtained from GROMACS -i,
  --index <.ndx>                Index file
  -s, --structure <.pdb/.gro>   Structure file
  -g, --group <string>          Group to match from index file to get the list of contacts

  -b, --begin <int>             Frame to start calculation from

  --pickle FILENAME             Save the output to a pickle file
  --hdf FILENAME                Save the output to a HDF file
  --csv FILENAME                Save the output to a CSV file

  --matrix FILENAME             Save a file containing the contact matrix

  --pdb FILENAME                Save a pdb file with percentage of frames
                                in the B-factor column

  -h, --help                    Show this screen.
"""

import click
import re
import pickle
import pandas
import collections
import numpy
from biopandas.pdb import PandasPdb


# TODO: Directly use biopandas data frame and remove the atom class

class Atom():
    """ Wrapper for MDTraj structure.atom object """

    def __init__(self, index, structure):
        atom = structure.df['ATOM']
        self.index = index
        self.segment = atom.iloc[index]['segment_id']
        self.chain = atom.iloc[index]['chain_id']
        self.residue = atom.iloc[index]['residue_name']
        self.res_id = atom.iloc[index]['residue_number']
        self.name = atom.iloc[index]['atom_name']

    def __repr__(self):
        return f'{self.segment:2} {self.chain} {self.residue} {self.res_id:3} {self.name:4}'

    def __str__(self):
        return f'{self.segment:2} {self.chain} {self.residue} {self.res_id:3} {self.name:4}'

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

    def __iter__(self):
        """ Segment Chain Residue ResSeq Atom """
        return iter([self.segment, self.chain, self.residue, self.res_id, self.name])


class Contact():
    """ Defines a contact between two atoms """

    def __init__(self, atom_indices, structure):
        self.atom1 = Atom(atom_indices[0], structure)
        self.atom2 = Atom(atom_indices[-1], structure)
        self.count = 0
        self.time_series = []

    def __repr__(self):
        return f'\n{repr(self.atom1)} -- {repr(self.atom2)}: {self.count}'

    def __str__(self):
        return f'\n{str(self.atom1)} -- {str(self.atom2)}: {self.count}'

    def __iter__(self):
        return iter(list(self.atom1) + list(self.atom2) + [self.count])


class Contacts():
    """ All hydrogen bonds """

    def __init__(self, structure) -> object:
        self.structure = PandasPdb().read_pdb(structure)
        self.bonds = []
        self.nframes = 1

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
                    self.bonds.append(Contact(indices, self.structure))
                if line == f'[ {group} ]\n':
                    save = True
        if not save:
            print('Specified group not found in index file')

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

    def to_pdb(self, filename):
        contact_dict = collections.defaultdict(list)
        for bond in self.bonds:
            contact_dict[bond.atom1.index].append(bond.time_series)
            contact_dict[bond.atom2.index].append(bond.time_series)
        for index, contacts in contact_dict.items():
            contacts = numpy.array(contacts).any(axis=0)
            b_factor = round(numpy.sum(contacts) / len(contacts) * 99.99, 2)
            self.structure.df['ATOM'].loc[index,'b_factor'] = b_factor
        self.structure.to_pdb(path=filename,
            records=None,
            gz=False,
            append_newline=True)

    def plot_matrix(self, df=None, n_frames=100, cutoff=0):

        def text(atoms):
            output = []
            for atom in atoms:
                if atom:
                    output.append('-'.join([str(_) for _ in atom]))
                else:
                    output.append('')
            return output

        if df:
            df = pandas.read_pickle(pickle_file)
        else:
            df = self.to_df

        # Convert to percentage of frames; default is to return counts
        df['Count'] = df['Count'] * 100 / n_frames

        # counts or percentage above which to show contacts
        df = df[df['Count'] > cutoff]

        atoms1 = list(df.groupby(['Chain1', 'Residue1', 'ResSeq1', 'Atom1']).groups.keys())
        atoms1 = sorted(atoms1, key=operator.itemgetter(0, 2, 1))
        atom_group2 = df.groupby(['Chain2', 'Residue2', 'ResSeq2', 'Atom2'])
        atoms2 = list(atom_group2.groups.keys())
        atoms2 = sorted(atoms2, key=operator.itemgetter(0, 2, 1))

        matrix = numpy.empty(shape=(len(atoms1), len(atoms2)))
        matrix[:] = numpy.NaN
        for group, atoms in atom_group2:
            j = atoms2.index(group)
            for index, row in atoms.iterrows():
                i = atoms1.index((row['Chain1'], row['Residue1'], row['ResSeq1'], row['Atom1']))
                matrix[i, j] = row['Count']

        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=text(atoms1),
            y=text(atoms2),
            colorscale='Greys', zmin=0, zmax=100,
            transpose=True,
        ))
        fig.update_layout(
            title=f'Contact Matrix for the interface between {sel1} and {sel2}',
        )
        py.plot(fig, filename=f'{prefix}_{sel1}-{sel2}_contact_matrix.html', auto_open=False)


class Hbond():
    """ Defines a hydrogen bond """

    def __init__(self, atom_indices, structure):
        self.donor = Atom(atom_indices[0], structure)
        self.hydrogen = Atom(atom_indices[1], structure)
        self.recipient = Atom(atom_indices[2], structure)
        self.count = 0
        self.time_series = []

    def __repr__(self):
        return f'\n{self.donor} -- {self.hydrogen} -> {self.recipient}: {self.count}'

    def __iter__(self):
        return iter(list(self.donor) + list(self.hydrogen) + list(self.recipient) + [self.count])

# Subclass the Contacts to get Hbonds Class
class Hbonds():
    """ All hydrogen bonds """

    def __init__(self, structure):
        self.structure = PandasPdb().read_pdb(structure)
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
                    self.bonds.append(Hbond(indices, self.structure))
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

    @property
    def to_df(self):
        columns = ['Segment1', 'Chain1', 'Residue1', 'ResSeq1', 'Atom1',
                   'Segment1-H', 'Chain1-H', 'Residue1-H', 'ResSeq1-H', 'Atom1-H',
                   'Segment2', 'Chain2', 'Residue2', 'ResSeq2', 'Atom2', 'Count']
        df = pandas.DataFrame(data=list(self), columns=columns)
        return df


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='hbond', context_settings=CONTEXT_SETTINGS)
@click.option('-x', '--xpm_file', metavar='<.xpm>', required=True,
              help='Contact/H-bond file obtained from GROMACS')
@click.option('-i', '--index_file', metavar='<.ndx>', required=True, help='Index file')
@click.option('-s', '--structure', metavar='<.pdb/.gro>', required=True,
              help='Structure file')
@click.option('-g', '--group', metavar='<string>', required=True,
              help='Group to match from index file to get the list of Contacts/H-bonds')
@click.option('-o', '--output', metavar='<.p>', help='Output Pickle file')
@click.option('--hbonds/--contacts',
              help='Select the type of input: contacts or H-bonds')
@click.option('--silent/--verbose', default=False, help='Choose silent to '
              'suppress the displaying to result to standard output')
@click.option('--save_pickle', help='Save the output to a pickle file')
@click.option('--save_hdf', help='Save the output to a HDF file')
@click.option('--save_csv', help='Save the output to a CSV file')
@click.option('--save_matrix',
              help='Save a file containing the contact matrix')
@click.option('--save_pdb', help='Save a pdb file with percentage of '
              'frames in the B-factor column')
def main(xpm_file, index_file, structure, group, hbonds, silent=False,
         output=None, save_pickle=None, save_hdf=None, save_csv=None,
         save_pdb=None, save_matrix=None):
    """Parse contacts evaluated by gmx hbonds"""
    if hbonds:
        molecular_contacts = Hbonds(structure)
    else:
        molecular_contacts = Contacts(structure)
    molecular_contacts.parse_indices(index_file=index_file, group=group)
    molecular_contacts.add_counts(xpm_file=xpm_file)

    if not silent:
        print(molecular_contacts)

    if output:
        pickle.dump(molecular_contacts, open(output, 'wb'))
    df = molecular_contacts.to_df
    if save_hdf:
        df.to_hdf(save_hdf, key='contacts')
    if save_csv:
        df.to_csv(save_csv)
    if save_pickle:
        df.to_pickle(save_pickle)
    if save_pdb:
        molecular_contacts.to_pdb(save_pdb)
    if save_matrix:
        molecular_contacts.to_pdb(save_matrix)


if __name__ == '__main__':
    main()
