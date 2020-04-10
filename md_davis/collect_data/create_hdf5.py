#! /usr/bin/env python
"""
Create an HDF file for collection analysis data

Usage:
  md_davis collect <HDF_file> <JSON> 
  md_davis collect -h | --help
  md_davis collect --version

Options:
  -h --help     Show this screen.
  --version     Show version.
"""

import h5py
import os
import docopt
import json

EPILOG = """
The attributes required are given in the example below:

{
    "label": "MD Simulation",
    "short_label": "MD",
    "html": "<i>MD Simulation</i>",
    "short_html": "<i>MD Simulation</i>",
    "protein": "protein name",
    "scientific_name": "some organism",
    "common_name": "common name",
    "sequence": "PUT/YOUR/SEQUENCE/HERE"
}

This information is primarily parsed to create labels for plots with 
this data file. 'sequence' is required to split the data into
chains.
"""

def main(arguments):
    if os.path.exists(arguments['<JSON>']):
        attributes = json.load(open(arguments['<JSON>'], 'r'))
    else:
        json_acceptable_string = arguments['<JSON>'].replace("'", "\"")
        attributes = json.loads(json_acceptable_string)

    with h5py.File(arguments['<HDF_file>'], 'a') as datafile:
        for key, value in attributes.items():
            datafile.attrs[key] = value

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    main(args)
