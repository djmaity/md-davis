#!/usr/bin/env python

"""
    This module parses the strings in .dat output of gmx do_dssp and
    plots the secondary structure exhibited by each residue.

    File name: secStr_per_residue.py
    Author: Dibyajyoti Maity
    Date created: 21 February 2019
    Python Version: 3.7
"""
import argparse
from collections import Counter
import pickle
import statistics

SECSTR_CODES = {'H':'α-helix',
                'G':'3_10-helix',
                'I':'π-helix',
                'B':'β-bridge',
                'E':'β strand',
                'T':'Turn',
                'S':'Bend',
                '~':'Loop',
}


def parse_dat(filename):
    """ Parse the DSSP data file """
    secondary_structure = []
    with open(filename, 'r') as dat_file:
        next(dat_file)
        for line in dat_file:
            secondary_structure.append(list(line.strip()))
    return secondary_structure


def get_counts_per_residue(ss_strings):
    output = []
    for ss in map(list, zip(*ss_strings)):
        output.append(Counter(ss))
    return output


def get_ss_percentage(ss_strings, start=0, end=None):
    output = {'Helix': [],
              'Sheet': [],
              'Turn + Bend': [],
              'Loop': [],
    }
    for ss in map(list, ss_strings):
        counter = Counter(ss[start:end])
        total = sum(counter.values())
        output['Helix'].append( (counter['H'] + counter['G'] + counter['I']) * 100 / total )
        output['Sheet'].append( (counter['E'] + counter['B']) * 100 / total )
        output['Turn + Bend'].append( (counter['T'] + counter['S']) * 100 / total )
        output['Loop'].append( (counter['~']) * 100 / total )
    output_string = ''
    for structure in output.keys():
        mean = statistics.mean(output[structure])
        std = statistics.stdev(output[structure])
        output_string += f'{structure:<12}: {mean:.2f} % +/- {std:.2f} %\n'
    return output_string


def get_arguments():
    """ Get input and output filename from the commandline """
    parser = argparse.ArgumentParser(description='Calculate and plot the'
        ' secondary structure exhibited by each residue in the output of'
        ' gmx do_dssp')
    parser.add_argument('filename',
                        help='.dat file output of gmx do_dssp')
    parser.add_argument('-o', '--output',
                        help='name and extension for output pickle file')
    return parser.parse_args()


def main():
    """ Parse and plot the DSSP data """
    args = get_arguments()
    data = parse_dat(args.filename)

    ss_percentage = get_ss_percentage(data)
    print(ss_percentage)

    if args.output:
        ss_counts = get_counts_per_residue(data)
        pickle.dump(ss_counts, open(args.output, 'wb'))


if __name__ == '__main__':
    main()
