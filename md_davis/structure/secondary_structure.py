#!/usr/bin/env python

"""
    This module parses the strings in .dat output of `gmx dssp` and
    plots the secondary structure exhibited by each residue.

    File name: secondary_strucutre.py
    Author: Dibyajyoti Maity
    Date created: 10 June 2024
    Python Version: 3.11
"""

import argparse
import pickle
import statistics
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Single letters used by DSSP to represent each secondary structure
DSSP_CODES = {'H': 'α-helix',
              'G': '3_10-helix',
              'I': 'π-helix',
              'P': 'Polyproline helix',
              'B': 'β-bridge',
              'E': 'β strand',
              'T': 'Turn',
              'S': 'Bend',
              '~': 'Loop'}

NUM_DSSP_CODES = 9  # Number of DSSP secondary structures

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


########################################################################################################################


def get_arguments():
    """ Get input and output filename from the commandline """
    parser = argparse.ArgumentParser(description='Calculate and plot the'
        ' secondary structure exhibited by each residue in the output of'
        ' gmx do_dssp')
    parser.add_argument('filename',
                        help='.dat file output of gmx do_dssp')
    parser.add_argument('-o', '--output',
                        help='name and extension for output image')
    return parser.parse_args()


def secondary_structure_percent(secstr, indices):
    """ Count the secondary structure exhibited by each residue throughout the trajectory. """
    frames = len(secstr)
    output = np.zeros((len_dssp_code, len(indices)))
    for i, ind in enumerate(indices):
        counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0,
                  'I': 0, 'T': 0, 'S': 0, '~': 0}
        for f in range(frames):
            counts[secstr[f][ind]] = counts[secstr[f][ind]] + 1
        for j in range(len_dssp_code):
            output[j][i] = counts[code[j]]
    return np.divide(output, frames / 100)


def do_dssp_count(filename, start=0, end=None):
    """ Parse the DSSP data """
    with open(filename, 'r') as dat_file:
        lines = dat_file.readlines()
        lines.pop(0)
        num_lines = len(lines)
    num_chains = len(lines[0].split('='))
    secstr = [['' for s in range(num_lines)] for n in
                range(num_chains)]
    for l, line in enumerate(lines):
        split_lines = line.strip().split('=')
        for chain in range(num_chains):
            if not end:
                secstr[chain][l] = split_lines[chain][start: end]
            else:
                secstr[chain][l] = split_lines[chain][start:]
    output = []
    for chain_secstr in secstr:
        output.append(secondary_structure_percent(chain_secstr, range(len(chain_secstr[0]))))
    return output


def plot_structure_counts(axes, counts, scale=1):
    """ Plot the stacked bar graph for the secondary structure counts """
    labels = [r'\textbf{$\alpha$-helix}',
              r'\textbf{$3_{10}$ helix}',
              r'\textbf{$\pi$-helix}',
              r'\textbf{residue in isolated $\beta$-bridge}',
              r'\textbf{$\beta$ strand}',
              r'\textbf{hydrogen bonded turn}',
              r'\textbf{bend}',
              r'\textbf{loop or irregular}',
              ]
    alpha = 0.4
    colors = [(1, 0, 0, alpha), (1, 0, 1, alpha), (0, 0, 1, alpha), (0, 1, 0, alpha),
              (0, 1, 1, alpha), (1, 1, 0, alpha), (0, 0, 0, alpha), (0.8, 0.8, 0.8, alpha)]
    width = 1
    cumulative = counts[0] * scale
    ind = list(range(len(cumulative)))
    axes.bar(ind, cumulative, width, label=labels[0], fc=colors[0],
             ec=(0, 0, 0, 0))
    for i in range(1, 8):
        axes.bar(ind, counts[i] * scale, width, bottom=cumulative,
                 label=labels[i], fc=colors[i], ec=(0, 0, 0, 0))
        cumulative = cumulative + (counts[i] * scale)
    axes.set_xlim(-0.5, ind[-1] + 1)


def main():
    """ Parse and plot the DSSP data """
    args = get_arguments()
    counts = do_dssp_count(args.filename)

    font = {'family' : 'serif',
            'weight' : 'bold',
           }
    matplotlib.rc('font', **font)
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

    if len(counts) > 1:
        fig, axes = plt.subplots(len(counts), 1, sharey=True)
        axis = axes.flatten()
        fig.set_size_inches(8, 2.5*len(counts) + 1)
    else:
        fig, axes = plt.subplots(1, 1, sharey=True)
        axis = [axes]
        fig.set_size_inches(8, 4.5)

    for i, count in enumerate(counts):
        plot_structure_counts(axis[i], count)
        handles, lbs = axis[i].get_legend_handles_labels()

    axis[-1].set_xlabel(r'\textbf{Residue Numbers}')
    lgd = axis[0].legend(handles, lbs, ncol=4, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    fig.text(0.05, 0.5, r'\textbf{\% of Frames with Specified Secondary Structure}',
             rotation='vertical', va='center')

    if args.output:
        fig.savefig(args.output, dpi=600, bbox_extra_artists=([lgd]), bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    main()
