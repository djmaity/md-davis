"""
Align the dataframes created by 'md_davis residue dataframe' for plotting.

Usage:  md_davis residue aligned --alignment FILE --map <dictionary> OUTPUT

Options:
  -h, --help                        Show this screen.
  -a, --alignment FILE              Alignment file in CLUSTAL format
  -m, --map <dictionary>            A dictionary mapping the labels in the
                                    alignment file with the corresponding
                                    pickle file containing residue dataframe.
"""

import docopt
import json
import pandas
import pickle
import re
import numpy
import collections
import os
import json

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

    assert len(set([len(_) for _ in alignment.values()])) == 1

    return alignment


def main(argv=None):
    if argv:
        args = docopt.docopt(__doc__, argv=argv)
    else:
        args = docopt.docopt(__doc__)

    if os.path.exists(args['--map']):
        mapping = json.load(open(args['--map'], 'r'))
    else:
        json_acceptable_string = args['--map'].replace("'", "\"")
        mapping = json.loads(json_acceptable_string)

    aln_dict = parse_alignment(args['--alignment'])
    aln_seq_length = len(next(iter(aln_dict.values())))  # Length of sequences in alignment file

    aligned_dataframe = collections.defaultdict(dict)
    for alignment_label, pickle_file in mapping.items():
        data = pickle.load(open(pickle_file, 'rb'))
        prefix = data['prefix']
        aligned_dataframe[prefix] = data
        if aln_dict:
            if alignment_label not in aln_dict:
                print(f'{alignment_label} not in alignment file')
                return
            for chain, df in data['data'].items():
                out_df = pandas.DataFrame(
                    data=numpy.zeros((aln_seq_length, df.shape[1])),
                    columns=df.columns)
                out_df['sequence', 'resn'] = '-'
                i, j = 0, 0
                for ressidue in aln_dict[alignment_label]:
                    if ressidue != '-':
                        assert ressidue == df.iloc[j].sequence.resn
                        out_df.iloc[i] = df.iloc[j]
                        j += 1
                    i += 1
                aligned_dataframe[prefix]['data'][chain] = out_df
    pickle.dump( aligned_dataframe, open(args['OUTPUT'], 'wb') )


if __name__ == "__main__":
    main()
