import argparse
import string
from Bio import Phylo


def get_y_positions(tree):
    """Create a mapping of each clade to its vertical position.

    Dict of {clade: y-coord}.
    Coordinates are negative, and integers for tips.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = dict((tip, maxheight - i)
                for i, tip in enumerate(tree.get_terminals()))

    # Internal nodes: place at midpoint of children
    def calc_row(clade):
        for subclade in clade:
            if subclade not in heights:
                calc_row(subclade)
        # Closure over heights
        heights[clade] = (heights[clade.clades[0]] +
                        heights[clade.clades[-1]]) / 2.0

    if tree.root.clades:
        calc_row(tree.root)
    return heights


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    if len(node_path) == 1:
        return tree.root
    elif len(node_path) > 1:
        return node_path[-2]
    else:
        return


def tikz_phylogram(tree, x_scale=20, y_scale=1, text_gap=0.3):
    """ Make tikz/pgf phylogram for the provided tree """
    letter = iter(string.ascii_lowercase)
    output_nodes, output_edges, output_text = '', '', ''

    leaves = tree.get_terminals()
    shapes = dict((_, 'in' + string.ascii_lowercase[i]) for i, _ in enumerate(tree.get_nonterminals()))

    x_dict = tree.depths()
    y_dict = get_y_positions(tree)
    for node, x in x_dict.items():
        x = x*x_scale
        y = y_dict[node]*y_scale
        label = '' if not node.name else node.name
        if node in shapes:
            shape = shapes[node]
        else:
            shape = next(letter)
        output_nodes += f'\\node[anchor=west] ({shape}) at ({x:.3f}, {y:.3f}) {{{label}}};\n'
        parent = get_parent(tree, node)
        # Draw branch to the parent node
        if parent in shapes:
            if node in leaves:
                output_edges += f'\draw  ({shapes[parent]}.center) |- ({shape});\n'
            else:
                output_edges += f'\draw  ({shapes[parent]}.center) |- ({shape}.center);\n'
            # Place distances on the branches
            x_text = (x + x_dict[parent]*x_scale) / 2.0
            y_text = y + text_gap
            output_text += f'\\node () at ({x_text:.3f}, {y_text:.3f}) {{{node.branch_length}}};\n'

    return '\\begin{tikzpicture}[sloped]\n' + output_nodes + output_edges +\
        output_text + '\end{tikzpicture}\n'


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('file', help='Input phylogentic tree in Newick format')
    parser.add_argument('-o', '--output', help='output filename')
    args = parser.parse_args()
    output = tikz_phylogram(Phylo.read(args.file, 'newick'))
    print(output)

if __name__ == "__main__":
    main()