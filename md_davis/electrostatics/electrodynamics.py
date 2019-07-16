#! /usr/bin/env python
import sys
import argparse
import os
import fnmatch
import pymol


def initialize_pymol(window=False):
    """ The commands necessary to get pymol running """
    if window:
        pymol.pymol_argv = ['pymol']
    else:
        pymol.pymol_argv = ['pymol','-cQ']
    pymol.finish_launching()
    return pymol.cmd


def main():
    parser = argparse.ArgumentParser(description='Create pymol session showing dynamics of electric field lines.')
    parser.add_argument('-p', '--potential', required=True,
                        help='Directory for potential files')
    parser.add_argument('-s', '--structure', required=True,
                        help='Directory for PDB structures')
    parser.add_argument('-n', '--name', required=True, help='Molecule name')
    parser.add_argument('--surface', default=False, action='store_true',
                        help='Show molecular surface')
    parser.add_argument('--ss_color', default=False, action='store_true',
                        help='Color the structure based on secondary structure')
    parser.add_argument('--spacing', default=4, type=int,
                        help='Spacing of field lines')
    parser.add_argument('-t', '--time_step', default=1, type=int,
                        help='time step for frames to show')
    parser.add_argument('-l', '--length', default=10, type=int,
                        help='length of field lines')
    parser.add_argument('-o', '--output', help='output session file name')
    parser.add_argument('--hide', default=True, action='store_false', help='Hide PyMOL window')
    args = parser.parse_args()

    cmd = initialize_pymol(args.hide)
    cmd.space('cmyk')
    cmd.set('bg_rgb', 'white')
    cmd.set('transparency', 0.5)
    cmd.set('movie_fps', 20)
    cmd.set('gradient_spacing', args.spacing)
    cmd.set('gradient_min_length', args.length)
    cmd.set('gradient_max_length', args.length + 2)

    state = 1
    states = []
    for root, dirs, files in os.walk(args.potential):
        for fname in files:
            if fnmatch.fnmatch(fname, "*.cub"):
                basename = os.path.splitext(os.path.basename(fname))[0]
                frame = int(basename.split('_')[-1])
                if frame % args.time_step == 0:
                    potential = os.path.join(root, fname)
                    structure =f'{args.structure}/{basename}.pdb'
                    cmd.load(structure, object=f'{args.name}_structure', state=state)
                    cmd.load(potential, object=f'{args.name}_potential_{state}')
                    cmd.ramp_new(name=f'{args.name}_ramp_{state}',
                                map_name=f'{args.name}_potential_{state}',
                                range=[-10, -1, 0, 1, 10],  # Cutoff for potential
                                color=['red', 'orange', 'gray', 'cyan', 'blue'],
                                )
                    cmd.gradient(f'{args.name}_field_lines',
                                f'{args.name}_potential_{state}',
                                state=state)
                    cmd.color(f'{args.name}_ramp_{state}', f'{args.name}_field_lines')
                    states.append(state)
                    state += 1
    if args.ss_color:
        cmd.color('orange', 'ss h')
        cmd.color('green', 'ss s')
        cmd.color('cyan', 'ss l+''')
    if args.surface:
        # TODO: Fix color
        cmd.copy(f'{args.name}_surface', f'{args.name}_structure')
        cmd.show_as('surface', f'{args.name}_surface')
        for state in states:
            cmd.set('surface_color', f'{args.name}_ramp_{state}', f'state {state}')

    cmd.group(f'{args.name}_potentials', f'{args.name}_potential*')
    cmd.group(f'{args.name}_ramps', f'{args.name}_ramp*')
    cmd.group(args.name, f'{args.name}_structure {args.name}_surface {args.name}_field_lines {args.name}_potentials {args.name}_ramps')
    if args.output:
        cmd.save(args.output, format='pse')


if __name__ == "__main__":
    main()
