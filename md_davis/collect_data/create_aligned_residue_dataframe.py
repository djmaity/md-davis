import argparse
import json
import pandas
import pickle
import re
import numpy
import collections


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


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('json', type=argparse.FileType('r'),
                        help='JSON file containing input filenames')
    args = parser.parse_args()

    json_dict = json.load(args.json)

    if 'alignment' in json_dict:
        aln_dict = parse_alignment(json_dict['alignment'])
        aln_seq_length = len(next(iter(aln_dict.values())))  # Length of sequences in alignment file
    else:
        aln_dict = None
        aln_seq_length = 0

    aligned_dataframe = collections.defaultdict(dict)
    for prefix, location in json_dict['locations'].items():
        data = pickle.load(open(location, 'rb'))
        assert prefix == data['prefix']
        aligned_dataframe[prefix] = data
        if aln_dict:
            if prefix not in aln_dict:
                print(f'{prefix} not in alignment file')
                return
            for chain, df in data['data'].items():
                out_df = pandas.DataFrame(
                    data=numpy.zeros((aln_seq_length, df.shape[1])),
                    columns=df.columns)
                out_df['sequence', 'resn'] = '-'
                i, j = 0, 0
                for ressidue in aln_dict[prefix]:
                    if ressidue != '-':
                        assert ressidue == df.iloc[j].sequence.resn
                        out_df.iloc[i] = df.iloc[j]
                        j += 1
                    i += 1
                aligned_dataframe[prefix]['data'][chain] = out_df
    pickle.dump( aligned_dataframe, open(json_dict['output'], 'wb') )


if __name__ == "__main__":
    main()
