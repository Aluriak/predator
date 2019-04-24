"""entry point for predator.

"""

import os
import argparse
from . import graph as graph_module
from . import predator, utils, __version__, print_info


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
    # meta
    parser.add_argument('--version', '-V', action='version', version=__version__)
    parser.add_argument('--info', action='store_true',
                        help="Print general info about the input graph, and exit")
    # arguments
    parser.add_argument('--targets-file', '-tf', type=str, default=None,
                        help="file containing one target per line")
    parser.add_argument('--targets', '-t', nargs='*', type=str, default=[],
                        help="targets to activate in the graph")
    parser.add_argument('--seeds-file', '-sf', type=str, default=None,
                        help="file containing one seed per line")
    parser.add_argument('--seeds', '-s', nargs='*', type=str, default=[],
                        help="seeds already activated in the graph")
    parser.add_argument('--forbidden-seeds-file', '-fsf', type=str, default=None,
                        help="file containing one forbidden seed per line")
    parser.add_argument('--forbidden-seeds', '-fs', nargs='*', type=str, default=[],
                        help="metabolites that cannot be seed in the graph")
    parser.add_argument('--export', '-e', type=str, default=None,
                        help="file export of the ASP instance")
    parser.add_argument('--visualize', '-v', type=str, default=None,
                        help="png file to render the input graph in (default: don't render)")
    parser.add_argument('--visualize-without-reactions', '-vr', type=str, default=None,
                        help="png file to render, without reactions, the input graph in (default: don't render)")
    # flags
    parser.add_argument('--union', action='store_true',
                        help="Print the union of all solutions")
    parser.add_argument('--intersection', action='store_true',
                        help="Print the intersection of all solutions")
    parser.add_argument('--no-topological-injection', action='store_true',
                        help="Do not use topological injection found in sbml data")
    parser.add_argument('--semantic-injection', action='store_true',
                        help="Use semantic injection found in sbml data")
    parser.add_argument('--targets-are-forbidden', '-taf', action='store_true',
                        help="Targets are added to forbidden seeds")

    return parser


def parse_args(args: iter = None) -> dict:
    return cli_parser().parse_args(args)

def get_all_ids(inline, fname) -> frozenset:
    "Return the set of all seeds, targets or forbidden seeds"
    elements = frozenset(inline)
    if fname:
        with open(fname) as fd:
            elements = elements | frozenset(line for line in map(str.strip, fd) if line)
    return elements


if __name__ == '__main__':
    args = parse_args()
    print('DEBUG args =', args)
    if args.info:
        print_info(args.infile, args.visualize, args.visualize_without_reactions)
        exit()
    graph = graph_module.graph_from_file(
        args.infile,
        use_topological_injections=not args.no_topological_injection,
        use_semantic_injection=args.semantic_injection
    )
    if args.export:
        sccs, _ = predator.compute_sccs(graph)
        scc_repr = '\n'.join(f'scc({scc_name},{node}).'
                            for scc_name, nodes in sccs.items()
                            for node in nodes)
        with open(args.export, "w") as f:
            f.write(graph + '\n' + scc_repr)
    # compute available targets (union of --targets and --targets-file)
    #  (and do the same for (forbidden) seeds)
    targets = get_all_ids(args.targets, args.targets_file)
    seeds = get_all_ids(args.seeds, args.seeds_file)
    forbidden_seeds = get_all_ids(args.forbidden_seeds, args.forbidden_seeds_file)
    if args.targets_are_forbidden:  forbidden_seeds |= targets
    # main work
    union_over_seeds, intersection_over_seeds = set(), set()
    first_loop = True  # only true during the first iteration of the next loop
    for idx, seeds in enumerate(predator.search_seeds(graph, seeds, forbidden_seeds, targets), start=1):
        if args.union:  union_over_seeds |= seeds
        if args.intersection:
            if first_loop:  intersection_over_seeds = seeds
            else:           intersection_over_seeds &= seeds
        repr_seeds = ', '.join(map(str, seeds))
        print(f"Solution {idx}:\n{repr_seeds}\n")
        first_loop = False
    print('end of solutions.')
    if args.union:  print('\nUnion:', ', '.join(union_over_seeds))
    if args.intersection:  print('\nIntersection:', ', '.join(intersection_over_seeds))
    if args.visualize:
        print('Input graph rendered in', utils.render_network(graph, args.visualize, with_reactions=True))
    if args.visualize_without_reactions:
        print('Input graph rendered in', utils.render_network(graph, args.visualize_without_reactions, with_reactions=False))
