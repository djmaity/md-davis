"""
Create a pandas dataframe from the information in HDF5 file for plotting

Usage:  md_davis create [options] HDF_FILE OUTPUT

Options:
  -h, --help                        Show this screen.
  -v, --version                     Show version.
  --prefix <string>                 Prefix used in the alignment file
  -a, --annotations <JSON_FILE>     JSON file containing annotations to mark on the plot
  --potentials <directory>          Dicrectory containing potential files
  --pdb_potentials <pot_file>       .pot file obtained from electrostatic calculation


"""

import argparse
import numpy
import h5py
import pandas
import json
# import re
# import sys
# import scipy.stats as statistics
import os
import pickle
import json
import collections
import docopt

# Import from my packages
from md_davis.utils import secStr_counts


def parse_potential(potential_file):
    df = pandas.read_csv(potential_file, skiprows=12, delim_whitespace=True, skipfooter=2,
        dtype={'resSeq': int}, engine='python',
        names=['name', 'resName', 'chainID', 'resSeq', 'potential',
            'reaction', 'coulomb', 'Ex', 'Ey', 'Ez'],
    )
    df = pandas.read_fwf(potential_file, skiprows=12, skipfooter=2,
        dtype={'resSeq': int}, engine='python',
        names=['name', 'resName', 'chainID', 'resSeq', 'potential',
            'reaction', 'coulomb', 'Ex', 'Ey', 'Ez'],
        widths=[5, 3, 3, 9, 10, 10, 10, 10, 10, 10]
    )
    output = {}
    chain = 0
    for _, data in df.groupby(['chainID'], as_index=False):
        # grouped_df = data.groupby(['resSeq', 'resName'], as_index=False)['potential']
        grouped_df = data.groupby(['resSeq'], as_index=False)['potential']
        potential =  grouped_df.sum()
        potential.rename(columns={'potential':'total'}, inplace=True)
        potential['mean'] =  grouped_df.mean()['potential']
        output[f'chain {chain}'] = potential
        chain += 1
    return output


def residue_dataframes(hdf_file, potentials=None, pdb_potentials=None):
    """ Create a pandas dataframe for Residue data to plot """
    output_dict = {}
    with h5py.File(hdf_file, 'r') as hdf5_file:
        for ch, chain_sequence in enumerate(hdf5_file.attrs['sequence'].split('/')):
            array_to_concat = []
            keys_to_concat = []
            # Add sequence
            sequence = pandas.DataFrame(data=enumerate(chain_sequence, start=1),
                                        columns=['resi', 'resn'])
            array_to_concat.append(sequence)
            keys_to_concat.append('sequence')
            # Add scondary structure to dtaframe from HDF5 file
            if f'secondary_structure/counts_per_residue/chain {ch}' in hdf5_file:
                ss_counts = hdf5_file[f'secondary_structure/counts_per_residue/chain {ch}']
                secondary_structure_dict = {}
                for code in secStr_counts.SECSTR_CODES.keys():
                    secondary_structure_dict[code] = ss_counts[code]
                secondary_structure_df = pandas.DataFrame(secondary_structure_dict)
                array_to_concat.append(secondary_structure_df)
                keys_to_concat.append('secondary_structure')
            # Add RMSF to dtaframe from HDF5 file
            if f'rmsf/chain {ch}' in hdf5_file:
                rmsf_df = pandas.DataFrame(hdf5_file[f'rmsf/chain {ch}'][:, 1])
                array_to_concat.append(rmsf_df)
                keys_to_concat.append('rmsf')
            # Add dihedral standard deviation to dtaframe from HDF5 file
            if f'dihedral_standard_deviation/chain {ch}' in hdf5_file:
                dih_sd_dict = {}
                for angle in ['phi', 'psi', 'omega']:
                    dih_sd_dict[angle] = hdf5_file[f'dihedral_standard_deviation/chain {ch}'][angle]
                dih_sd_df = pandas.DataFrame(dih_sd_dict)
                array_to_concat.append(dih_sd_df)
                keys_to_concat.append('dihedral_standard_deviation')
            # Combine all of the above
            output_dict[f'chain {ch}'] = pandas.concat(array_to_concat,
                                                     keys=keys_to_concat,
                                                     axis=1)
    # Add potential
    if pdb_potentials:
        pdb_potential_df = parse_potential(pdb_potentials)

    if potentials:
        surf_pot = collections.defaultdict(list)
        for file in os.listdir(potentials):
            if file.endswith('.pot'):
                print('Parsing: ', file)
                pot_file = os.path.join(potentials, file)
                for chain, pot_df in parse_potential(pot_file).items():
                    surf_pot[chain].append(pot_df)

        for chain, data_frame in output_dict.items():
            total_df = pandas.DataFrame(data_frame['sequence', 'resi'].values,
                                            columns=['resi'])
            mean_df = pandas.DataFrame(data_frame['sequence', 'resi'].values,
                                            columns=['resi'])
            total_column = 0
            mean_column = 0
            for potential_df in surf_pot[chain]:
                total_df[str(total_column)] = 0
                mean_df[str(mean_column)] = 0
                for _, row in potential_df.iterrows():
                    total_df.loc[total_df.resi == row.resSeq, str(total_column)] = row['total']
                    mean_df.loc[mean_df.resi == row.resSeq, str(mean_column)] = row['mean']
                total_column += 1
                mean_column += 1

            del total_df['resi']
            del mean_df['resi']

            surface_potential_df = pandas.DataFrame()
            surface_potential_df['mean_total'] = total_df.mean(skipna=True, axis=1)
            surface_potential_df['std_total'] = total_df.std(skipna=True, axis=1)
            surface_potential_df['mean_mean'] = mean_df.mean(skipna=True, axis=1)
            surface_potential_df['std_mean'] = mean_df.std(skipna=True, axis=1)

            columns = [
                ('surface_potential', 'mean_total'),
                ('surface_potential', 'std_total'),
                ('surface_potential', 'mean_mean'),
                ('surface_potential', 'std_mean'),
            ]

            # Add PDB potential
            if pdb_potentials:
                surface_potential_df['pdb_total'] = 0
                surface_potential_df['pdb_mean'] = 0
                for _, pdb_row in pdb_potential_df[chain].iterrows():
                    surface_potential_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_total'] = pdb_row['total']
                    surface_potential_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_mean'] = pdb_row['mean']
                columns += [
                    ('surface_potential', 'pdb_total'),
                    ('surface_potential', 'pdb_mean'),
                ]
            surface_potential_df.columns = pandas.MultiIndex.from_tuples(columns)
            output_dict[chain] = data_frame.join(surface_potential_df)
    return output_dict


def main(args):
    output = {}
    if args['--prefix']:
        output['prefix'] = args['--prefix']
    output['data'] = residue_dataframes(hdf_file=args['HDF_FILE'],
                                        potentials=args['--potentials'],
                                        pdb_potentials=args['--pdb_potentials'])
    if args['--annotations']:
        with open(args['--annotations'], 'r') as json_file:
            output['annotations'] = json.load(json_file)
    pickle.dump(output, open(args['OUTPUT'], 'wb'))


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    main(arguments)
