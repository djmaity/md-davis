import collections
import pandas
import numpy
import pickle
import string
import sys

def get_ss_chain_slices(secondary_structure):
    slices = {}
    prev = 0
    last = len(secondary_structure)
    chain = iter(string.ascii_uppercase)
    for current, counter in enumerate(secondary_structure):
        if '=' in counter:
            slices[next(chain)] = (prev, current)
            prev = current + 1
    slices[next(chain)] = (prev, last)
    return slices


def add_secondary_structure(df_dict, ss_counts, slice_dict):
    for chain, (i, j) in slice_dict.items():
        for code, name in secStr_counts.SECSTR_CODES.items():
            ss_count = [_[code] for _ in ss_counts[i:j] ]
            df_dict[chain][name] = pandas.Series(ss_count + [0], index=df_dict[chain].index)


def main():
    for chain, seq in hemoglobin.HBA_SEQUENCE.items():
        df = hemoglobins['deoxy-HbA'][chain]
        df['ResSeq'] = df['Residue']        
        df['Residue'] = pandas.Series(list(seq), index=df.index)
        hemoglobins['deoxy-HbA'][chain] = df[['ResSeq', 'Residue', 'RMSF']]
        df = hemoglobins['deoxy-HbA (Crystal)'][chain]
        df['ResSeq'] = df['Residue']
        df['Residue'] = pandas.Series(list(seq), index=df.index)
        hemoglobins['deoxy-HbA (Crystal)'][chain] = df[['ResSeq', 'Residue', 'RMSF']]

    for chain, seq in hemoglobin.HBS_SEQUENCE.items():
        df = hemoglobins['deoxy-HbS'][chain]
        df['ResSeq'] = df['Residue']        
        df['Residue'] = pandas.Series(list(seq), index=df.index)
        hemoglobins['deoxy-HbS'][chain] = df[['ResSeq', 'Residue', 'RMSF']]
        df = hemoglobins['deoxy-HbS (Crystal)'][chain]
        df['ResSeq'] = df['Residue']        
        df['Residue'] = pandas.Series(list(seq), index=df.index)
        hemoglobins['deoxy-HbS (Crystal)'][chain] = df[['ResSeq', 'Residue', 'RMSF']]
        
    dihedral_sd_file = folder + '15_hemoglobin_dihedral_sd/hemoglobin_crystal_dihedral_sd.p'
    dihedral_sd = pickle.load(open(dihedral_sd_file, 'rb'))
    dihedral_slices = {'A': (0, 280),
                    'B': (280, 570),
                    'C': (570, 850),
                    'D': (850, 1140),
    }
    for chain, (start, end) in dihedral_slices.items():
        for name in hemoglobins.keys():
            phi = numpy.concatenate([[0], dihedral_sd[name][start+1:end:2]])
            psi = numpy.concatenate([dihedral_sd[name][start:end:2], [0]])
            df = hemoglobins[name][chain]
            df['phi'] = pandas.Series(phi, index=df.index)
            df['psi'] = pandas.Series(psi, index=df.index)


    ss_data = pickle.load(open(folder + '16_hemoglobin_crystal_secondary_structure/2hhb_crystal_ss_counts.p', 'rb'))
    add_secondary_structure(hemoglobins['deoxy-HbA (Crystal)'], ss_data, get_ss_chain_slices(ss_data))

    ss_data = pickle.load(open(folder + '16_hemoglobin_crystal_secondary_structure/2hbs_crystal_ss_counts.p', 'rb'))
    add_secondary_structure(hemoglobins['deoxy-HbS (Crystal)'], ss_data, get_ss_chain_slices(ss_data))

    pickle.dump(hemoglobins, open('hemoglobin_crystal_plot_data.p', 'wb'))

if __name__ == "__main__":
    main()