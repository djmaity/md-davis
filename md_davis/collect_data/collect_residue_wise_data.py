import argparse
import numpy
import h5py
import pandas
import json
import re
import scipy.stats as statistics
import sys
import os
import pickle
import json
import collections

# Import from my packages
from md_davis import secStr_counts


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


def main():
    parser = argparse.ArgumentParser(description='Collect residue wise data into pandas dataframes')
    parser.add_argument('-f', '--file', help='HDF5 file containing data')
    parser.add_argument('-p', '--prefix', default='', help='Prefix for the sequence in alignment file')
    parser.add_argument('-d', '--directory', help='Directory containing surface potential files')
    parser.add_argument('-e', '--electrostatics', help='Surface potentials for the PDB structure')
    parser.add_argument('-a', '--annotations', type=argparse.FileType('r'), help='Residues to mark on the plot')
    parser.add_argument('-o', '--output', type=argparse.FileType('wb'), help='Output pickle file')
    args = parser.parse_args()

    output = {'prefix': args.prefix}
    output['data'] = {}
    if args.annotations:
        output['annotations'] = json.load(args.annotations)

    with h5py.File(args.file, 'r') as hdf5_file:
        for ch, chain_sequence in enumerate(hdf5_file.attrs['sequence'].split('/')):
            # Add sequence
            sequence = pandas.DataFrame(data=enumerate(chain_sequence, start=1),
                                        columns=['resi', 'resn'])
            # Add scondary structure to dtaframe from HDF5 file
            ss_counts = hdf5_file[f'secondary_structure/counts_per_residue/chain {ch}']
            secondary_structure_dict = {}
            for code in secStr_counts.SECSTR_CODES.keys():
                secondary_structure_dict[code] = ss_counts[code]
            secondary_structure_df = pandas.DataFrame(secondary_structure_dict)
            # Add RMSF to dtaframe from HDF5 file
            rmsf_df = pandas.DataFrame(hdf5_file[f'rmsf/chain {ch}'][:, 1])
            # Add dihedral standard deviation to dtaframe from HDF5 file
            dih_sd_dict = {}
            for angle in ['phi', 'psi', 'omega']:
                dih_sd_dict[angle] = hdf5_file[f'dihedral_standard_deviation/chain {ch}'][angle]
            dih_sd_df = pandas.DataFrame(dih_sd_dict)
            # Combine all of the above
            output['data'][f'chain {ch}'] = pandas.concat([sequence, secondary_structure_df, rmsf_df, dih_sd_df],
                keys=['sequence', 'secondary_structure', 'rmsf', 'dihedral_standard_deviation'], axis=1)

    # Add potential
    potentials = collections.defaultdict(list)
    for file in os.listdir(args.directory):
        if file.endswith('.pot'):
            print('Parsing: ', file)
            pot_file = os.path.join(args.directory, file)
            for chain, pot_df in parse_potential(pot_file).items():
                potentials[chain].append(pot_df)

    pdb_potential_df = parse_potential(args.electrostatics)

    for chain, data_frame in output['data'].items():
        total_df = pandas.DataFrame(data_frame['sequence', 'resi'].values,
                                           columns=['resi'])
        mean_df = pandas.DataFrame(data_frame['sequence', 'resi'].values,
                                           columns=['resi'])
        total_column = 0
        mean_column = 0
        for potential_df in potentials[chain]:
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
        surface_potential_df['pdb_total'] = 0
        surface_potential_df['mean_mean'] = mean_df.mean(skipna=True, axis=1)
        surface_potential_df['std_mean'] = mean_df.std(skipna=True, axis=1)
        surface_potential_df['pdb_mean'] = 0
        # Add PDB potential
        for _, pdb_row in pdb_potential_df[chain].iterrows():
            surface_potential_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_total'] = pdb_row['total']
            surface_potential_df.loc[data_frame['sequence']['resi'] == pdb_row.resSeq, 'pdb_mean'] = pdb_row['mean']

        columns=[('surface_potential', 'mean_total'),
                 ('surface_potential', 'std_total'),
                 ('surface_potential', 'pdb_total'),
                 ('surface_potential', 'mean_mean'),
                 ('surface_potential', 'std_mean'),
                 ('surface_potential', 'pdb_mean'),
                ]
        surface_potential_df.columns = pandas.MultiIndex.from_tuples(columns)
        output['data'][chain] = data_frame.join(surface_potential_df)
    pickle.dump(output, args.output)


if __name__ == "__main__":
    main()
