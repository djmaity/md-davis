#! /usr/bin/env python
import argparse
import fnmatch
import os
import re
import click

import pymol


def initialize_pymol(window=False):
    """ The commands necessary to get pymol running """
    if window:
        pymol.pymol_argv = ['pymol']
    else:
        pymol.pymol_argv = ['pymol','-cQ']
    pymol.finish_launching()
    return pymol.cmd


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(name='electrodynamics', context_settings=CONTEXT_SETTINGS)
@click.option('-n', '--name', required=True, help='Name for the PyMOL obeject to be created for the molecule. (No spaces!)')
@click.option('--surface', default=False, is_flag=True,
                        help='Show molecular surface')
@click.option('--ss_color', default=False, is_flag=True,
                        help='Color the structure based on secondary structure')
@click.option('--spacing', default=4, type=int,
                        help='Spacing of field lines')
@click.option('-t', '--time_step', default=1, type=int,
                        help='time step for frames to show')
@click.option('-l', '--length', default=10, type=int,
                        help='length of field lines')
@click.option('--light/--dark', help='Enable dark or light mode')
@click.option('-o', '--output', help='Name for saving the PyMOL session', type=click.Path())
@click.option('--hide', default=True, is_flag=True, help='Hide PyMOL window')
@click.argument('electrostatics_directory', metavar='DIRECTORY')
def main(electrostatics_directory, name, surface, ss_color, spacing, time_step, length, output, hide, light):
    """ Create pymol session showing dynamics of electric field lines.
    DIRECTORY   directory containing structures and electrostatics calculation from md_davis electrostatics
    """
    cmd = initialize_pymol(hide)

    # cmd.feedback("disable", "all", "actions")
    # cmd.feedback("disable", "all", "results")

    cmd.space('cmyk')
    cmd.set('transparency', 0.5)
    cmd.set('movie_fps', 5)
    cmd.set('gradient_spacing', spacing)
    cmd.set('gradient_min_length', length)
    cmd.set('gradient_max_length', length + 2)
    cmd.set('line_radius', 20)

    if light:
        cmd.set('bg_rgb', 'white')
        ramp_colors = ['red', 'gray', 'blue']
    else:
        cmd.set('bg_rgb', 'black')
        ramp_colors = ['red', 'white', 'blue']

    state = 1
    states = []
    regex = re.compile(r'\d+')
    for root, dirs, files in os.walk(electrostatics_directory):
        for fname in files:
            if fnmatch.fnmatch(fname, "*.cub"):
                basename = os.path.splitext(os.path.basename(fname))[0]
                frame = int(regex.findall(basename)[-1])

                if frame % time_step == 0:
                    potential = os.path.join(root, fname)
                    structure =f'{electrostatics_directory}/{basename}.pdb'
                    cmd.load(structure, object=f'{name}_structure', state=state)
                    cmd.load(potential, object=f'{name}_potential_{state}')
                    cmd.ramp_new(name=f'{name}_ramp_{state}',
                                map_name=f'{name}_potential_{state}',
                                range=[-5, 0, 5],  # Cutoff for potential
                                color=ramp_colors,
                                )
                    cmd.gradient(f'{name}_field_lines',
                                f'{name}_potential_{state}',
                                state=state)
                    cmd.color(f'{name}_ramp_{state}', f'{name}_field_lines')
                    states.append(state)
                    state += 1
    if ss_color:
        cmd.color('orange', 'ss h')
        cmd.color('green', 'ss s')
        cmd.color('cyan', 'ss l+''')
    if surface:
        # TODO: Fix color
        cmd.copy(f'{name}_surface', f'{name}_structure')
        cmd.show_as('surface', f'{name}_surface')
        for state in states:
            cmd.set('surface_color', f'{name}_ramp_{state}', f'state {state}')

    cmd.group(f'{name}_potentials', f'{name}_potential*')
    cmd.group(f'{name}_ramps', f'{name}_ramp*')
    cmd.group(name, f'{name}_structure {name}_surface {name}_field_lines {name}_potentials {name}_ramps')
    if output:
        cmd.save(output, format='pse')


if __name__ == "__main__":
    main()
