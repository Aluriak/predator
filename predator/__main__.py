"""entry point for predator.

"""

import os
import argparse
from . import graph as graph_module
from . import predator


def existant_file(filepath:str) -> str:
    """Argparse type, raising an error if given file does not exists"""
    if not os.path.exists(filepath):
        raise argparse.ArgumentTypeError("file {} doesn't exists".format(filepath))
    return filepath


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    # positionnals
    parser.add_argument('infile', type=existant_file,
                        help="SBML or ASP file containing the graph data")
    # arguments
    parser.add_argument('--targets-file', type=str, default=None,
                        help="file containing one target per line")
    parser.add_argument('--targets', nargs='*', type=set, default=[],
                        help="targets to activate in the graph")
    # flags
    parser.add_argument('--union', action='store_true',
                        help="Print the union of all solutions")
    parser.add_argument('--intersection', action='store_true',
                        help="Print the union of all solutions")

    return parser


def parse_args(args: iter = None) -> dict:
    return cli_parser().parse_args(args)


if __name__ == '__main__':
    args = parse_args()
    graph = graph_module.graph_from_file(args.infile)
    # compute available targets (union of --targets and --targets-file)
    targets = frozenset(args.targets)
    if args.targets_file:
        with open(args.targets_file) as fd:
            targets = targets | frozenset(line for line in map(str.strip, fd) if line)
    # main work
    union_over_seeds, intersection_over_seeds = set(), set()
    for idx, seeds in enumerate(predator.search_seeds(graph, targets), start=1):
        union_over_seeds |= seeds
        intersection_over_seeds &= seeds
        print(f"Solution {idx}:\n{solution}\n")
    print('end of solutions.')
    print('\nUnion:', ', '.join(union_over_seeds))
    print('\nIntersection:', ', '.join(intersection_over_seeds))
