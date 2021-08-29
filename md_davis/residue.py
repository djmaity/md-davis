"""
Create a pandas dataframe from the information in HDF5 file for plotting

Usage:  md_davis residue dataframe [options] HDF_FILE OUTPUT

Options:
  -h, --help                        Show this screen.
  -v, --version                     Show version.
  --prefix <string>                 Prefix used in the alignment file
  -a, --annotations <JSON_FILE>     JSON file containing annotations to mark on the plot
  --potentials <directory>          Dicrectory containing potential files
  --pdb_potentials <pot_file>       .pot file obtained from electrostatic calculation


"""

import h5py
import pandas
import os
import pickle
import numpy
import re
import toml
import collections
import click

# Import from my packages
from md_davis.common import secStr_counts

pandas.option_context('display.max_rows', None)


def residue_dataframe(hdf_file, potentials=None, pdb_potentials=None):
    """ Create a pandas dataframe for Residue data to plot """
    output_dict = {}
    for ch, chain_sequence in enumerate(hdf_file.attrs['sequence'].split('/')):
        array_to_concat = []
        keys_to_concat = []
        # Add sequence
        sequence = pandas.DataFrame(data=enumerate(chain_sequence, start=1),
                                    columns=['resi', 'resn'])
        array_to_concat.append(sequence)
        keys_to_concat.append('sequence')
        # Add scondary structure to dtaframe from HDF5 file
        if f'secondary_structure/counts_per_residue/chain {ch}' in hdf_file:
            ss_counts = hdf_file[
                f'secondary_structure/counts_per_residue/chain {ch}']
            secondary_structure_dict = {}
            for code in secStr_counts.SECSTR_CODES.keys():
                secondary_structure_dict[code] = ss_counts[code]
            secondary_structure_df = pandas.DataFrame(secondary_structure_dict)
            array_to_concat.append(secondary_structure_df)
            keys_to_concat.append('secondary_structure')
        # Add RMSF to dtaframe from HDF5 file
        if f'residue_property/rmsf/chain {ch}' in hdf_file:
            rmsf_df = pandas.DataFrame(
                hdf_file[f'residue_property/rmsf/chain {ch}'][:, 1])
            array_to_concat.append(rmsf_df)
            keys_to_concat.append('rmsf')
        # Add SASA to dtaframe from HDF5 file
        if f'residue_property/sasa/chain {ch}' in hdf_file:
            sasa_dict = {}
            for measure in ['average', 'standard_deviation']:
                sasa_dict[measure] = \
                hdf_file[f'residue_property/sasa/chain {ch}'][measure]
            sasa_df = pandas.DataFrame(sasa_dict)
            array_to_concat.append(sasa_df)
            keys_to_concat.append('sasa')
        # Add dihedral standard deviation to dtaframe from HDF5 file
        if f'dihedral_standard_deviation/chain {ch}' in hdf_file:
            dih_sd_dict = {}
            for angle in ['phi', 'psi', 'omega']:
                dih_sd_dict[angle] = \
                hdf_file[f'dihedral_standard_deviation/chain {ch}'][angle]
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
                    total_df.loc[
                        total_df.resi == row.resSeq, str(total_column)] = row[
                        'total']
                    mean_df.loc[mean_df.resi == row.resSeq, str(mean_column)] = \
                    row['mean']
                total_column += 1
                mean_column += 1

            del total_df['resi']
            del mean_df['resi']

            surface_potential_df = pandas.DataFrame()
            surface_potential_df['mean_total'] = total_df.mean(skipna=True,
                                                               axis=1)
            surface_potential_df['std_total'] = total_df.std(skipna=True,
                                                             axis=1)
            surface_potential_df['mean_mean'] = mean_df.mean(skipna=True,
                                                             axis=1)
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
                    surface_potential_df.loc[data_frame['sequence'][
                                                 'resi'] == pdb_row.resSeq, 'pdb_total'] = \
                    pdb_row['total']
                    surface_potential_df.loc[data_frame['sequence'][
                                                 'resi'] == pdb_row.resSeq, 'pdb_mean'] = \
                    pdb_row['mean']
                columns += [
                    ('surface_potential', 'pdb_total'),
                    ('surface_potential', 'pdb_mean'),
                ]
            surface_potential_df.columns = pandas.MultiIndex.from_tuples(
                columns)
            output_dict[chain] = data_frame.join(surface_potential_df)
    return output_dict


def parse_alignment(alignment_file):
    """ Parse and alignment file in CLUSTAL format"""
    alignment = {}
    with open(alignment_file) as alnFile:
        for line in alnFile.readlines()[1:]:
            m = re.match('(\S+)\s+(\S+)', line)
            if m:
                if m.group(1) in alignment:
                    alignment[m.group(1)] = alignment[m.group(1)] + m.group(2)
                else:
                    alignment[m.group(1)] = m.group(2)
    # assert len(set([len(_) for _ in alignment.values()])) == 1
    return alignment


# def align_residues(residue_dataframes, alignment, mapping=None):
#     """ Align residue dataframes based on sequence alignment """
#
#     for chain, alignment_file in alignment['alignment'].items():
#         aln_dict = parse_alignment(alignment)
#         if aln_dict is None:
#             return
#
#         # Length of sequences in alignment file
#         aln_seq_length = len(next(iter(aln_dict.values())))
#
#         aligned_dataframes = collections.defaultdict(dict)
#         for name, data in residue_dataframes.items():
#             aligned_dataframes[name] = data
#
#
#
#             out_df = pandas.DataFrame(
#                 data=numpy.zeros((aln_seq_length, df.shape[1])),
#                 columns=df.columns)
#             out_df['sequence', 'resn'] = '-'
#             i, j = 0, 0
#             for residue in aln_dict[name]:
#                 if residue != '-':
#                     if residue == df.iloc[j].sequence.resn:
#                         out_df.iloc[i] = df.iloc[j]
#                         j += 1
#                     else:
#                         print('Mismatch between sequence in HDF file'
#                               ' and alignment file.')
#                         return
#                 i += 1
#             aligned_dataframes[name][chain] = out_df
#     return aligned_dataframes


def align_residues(residue_dataframes, alignment, mapping=None):
    """ Align residue dataframes based on sequence alignment """

    print(residue_dataframes.keys())

    aligned_dataframes = residue_dataframes
    for chain, alignment_file in alignment['alignment'].items():
        aln_dict = parse_alignment(alignment_file)
        if aln_dict is None:
            return
        # Length of sequences in alignment file
        aln_seq_length = len(next(iter(aln_dict.values())))

        for name, seq_label in alignment['names'].items():
            df = residue_dataframes[name][chain]
            out_df = pandas.DataFrame(
                data=numpy.zeros((aln_seq_length, df.shape[1])),
                columns=df.columns)
            out_df['sequence', 'resn'] = '-'
            i, j = 0, 0
            for residue in aln_dict[seq_label]:
                if residue != '-':
                    if residue == df.iloc[j].sequence.resn:
                        out_df.iloc[i] = df.iloc[j]
                        j += 1
                    else:
                        print('Mismatch between sequence in HDF file'
                              ' and alignment file.')
                        return
                i += 1
            aligned_dataframes[name][chain] = out_df
    return aligned_dataframes



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='residue', context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--output', default='residue_dataframe.p')
@click.option('-a', '--annotations', type=str, default=None,
              help="Annotations")
@click.option('--alignment',
              help="Alignment file to align the residue data")
@click.argument('input_files', nargs=-1, required=True,
                type=click.Path(exists=True))
def main(input_files, output, alignment=None, annotations=None):
    res_dfs = collections.defaultdict(dict)
    for input_file in input_files:
        with h5py.File(input_file, 'r') as hdf_file:
            res_df = residue_dataframe(hdf_file=hdf_file)
            res_dfs['data'][hdf_file.attrs['name']] = res_df
    if annotations:
        with open(annotations, 'r') as annotations_file:
            res_dfs['annotations'] = toml.load(annotations_file)
    if alignment:
        alignment = toml.load(alignment)
        res_dfs['data'] = align_residues(res_dfs['data'], alignment=alignment)
    pickle.dump(res_dfs, open(output, 'wb'))


if __name__ == "__main__":
    main()
