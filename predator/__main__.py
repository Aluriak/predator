"""entry point for predator.

"""

import os
import time
import argparse
from . import prerun
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
    parser.add_argument('--visualize-dag', '-vd', type=str, default=None,
                        help="png file to render the DAG of SCC (default: don't render)")

    # flags groups
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--union', action='store_true',
                       help="Print the union of all solutions")
    group.add_argument('--intersection', action='store_true',
                       help="Print the intersection of all solutions")

    # flags
    parser.add_argument('--greedy', action='store_true',
                        help="Use greedy implementation for target search")
    parser.add_argument('--no-prerun', action='store_true',
                        help="Don't perform the preruning step")
    parser.add_argument('--no-topological-injection', action='store_true',
                        help="Do not use topological injection found in sbml data")
    parser.add_argument('--semantic-injection', action='store_true',
                        help="Use semantic injection found in sbml data")
    parser.add_argument('--targets-are-forbidden', '-taf', action='store_true',
                        help="Targets are added to forbidden seeds")
    parser.add_argument('--pareto', '-p', action='store_true',
                        help="Explore the pareto front of targets/seeds ratio")
    parser.add_argument('--pareto-full', '-pf', action='store_true',
                        help="Explore the pareto front of targets/seeds/targets-as-seeds ratios")
    parser.add_argument('--scc-with-asp', action='store_true',
                        help="Use ASP search for SCC mining")
    parser.add_argument('--verbose', action='store_true',
                        help="Print everything there is to know about the search ; beware the flood")



    return parser


def parse_args(args: iter = None) -> dict:
    return cli_parser().parse_args(args)

def get_all_ids(inline, fname) -> frozenset:
    "Return the set of all seeds, targets or forbidden seeds"
    elements = frozenset(inline)
    if fname:
        elements |= frozenset(utils.get_ids_from_file(fname))
    return elements


if __name__ == '__main__':
    args = parse_args()
    print('DEBUG args =', args)
    if args.info:
        print_info(args.infile, args.visualize, args.visualize_without_reactions)
        exit()
    supp_args = {}  # more arguments for search_seeds function
    if args.verbose:  supp_args['verbose'] = True

    # Extract the ASP representation from the graph file
    time_data_extraction = time.time()
    initial_graph = graph = graph_module.graph_from_file(
        args.infile,
        use_topological_injections=not args.no_topological_injection,
        use_semantic_injection=args.semantic_injection
    )
    time_data_extraction = time.time() - time_data_extraction

    # Define targets, seeds, enumeration mode,…
    targets = get_all_ids(args.targets, args.targets_file)
    seeds = get_all_ids(args.seeds, args.seeds_file)
    forbidden_seeds = get_all_ids(args.forbidden_seeds, args.forbidden_seeds_file)
    enum_mode = 'enumeration'
    if args.union:  enum_mode = 'union'
    if args.intersection:  enum_mode = 'intersection'
    if args.targets_are_forbidden:  forbidden_seeds |= targets
    if forbidden_seeds:
        initial_graph += '\n'.join(f'annot(lower,"{fs}","F").' for fs in forbidden_seeds)
    if targets:
        initial_graph += '\n'.join(f'annot(upper,"{tg}","T").' for tg in targets)

    # Prerun
    time_prerun = 'unperformed'
    if not args.no_prerun:
        time_prerun = time.time()
        graph, targets, unreachables, dead_viz = prerun.on_graph(graph, forbidden_seeds, targets)
        initial_graph += dead_viz  # indicate what are the dead parts in the vizualizations
        time_prerun = time.time() - time_prerun
        if unreachables and targets:
            print('Some targets cannot be reached:',
                  ', '.join(map(utils.unquoted, unreachables)))
        elif unreachables:  # all targets are unreachables
            print('No target is reachable. Abort.')

    # External computation of SCCs
    graph_filename = None if args.scc_with_asp else args.infile
    if args.export:
        sccs, _ = predator.compute_sccs(graph)
        scc_repr = '\n'.join(f'scc({scc_name},{node}).'
                            for scc_name, nodes in sccs.items()
                            for node in nodes)
        with open(args.export, "w") as f:
            f.write(graph + '\n' + scc_repr)

    # If needed, extract the sccs here, in order to time it
    time_sccs_extraction = 'unperformed'
    if not args.greedy:
        if args.visualize_dag:
            time_sccs_extraction = time.time()
            supp_args['sccs'], supp_args['scc_dag'] = predator.compute_sccs(graph, graph_filename=graph_filename, verbose=args.verbose)
            time_sccs_extraction = time.time() - time_sccs_extraction
            if args.targets:
                args.visualize_dag = predator.render_scc_dag(args.visualize_dag, supp_args['sccs'], supp_args['scc_dag'], targets)
            else:
                args.visualize_dag = predator.render_scc_dag(args.visualize_dag, supp_args['sccs'], supp_args['scc_dag'], targets = set())
                del supp_args['sccs']
                del supp_args['scc_dag']

    # Compute available targets (union of --targets and --targets-file)
    #  (and do the same for (forbidden) seeds).
    time_seed_search = time.time()
    for idx, seeds in enumerate(predator.search_seeds(graph, seeds, forbidden_seeds, targets, enum_mode=enum_mode, graph_filename=graph_filename, explore_pareto=args.pareto, pareto_no_target_as_seeds=args.pareto_full, greedy=args.greedy, **supp_args), start=1):
        repr_seeds = ', '.join(map(str, seeds))
        print(f"Solution {idx}:\n{repr_seeds}\n")
    print('end of solutions.')
    time_seed_search = time.time() - time_seed_search

    # Render the graph if asked to
    time_rendering = time.time()
    print('Rendering…')
    try:
        if args.visualize:
                print('-> Input graph rendered in', utils.render_network(initial_graph, args.visualize, with_reactions=True))
        if args.visualize_without_reactions:
            print('-> Input graph rendered in', utils.render_network(initial_graph, args.visualize_without_reactions, with_reactions=False))
        if args.visualize_dag:
            print('-> DAG of SCC rendered in', args.visualize_dag)
    except KeyboardInterrupt:
        print('-> Aborted!')
    else:
        print('-> Ok!')
    time_rendering = time.time() - time_rendering

    # Show the timers
    timers = {
        'DATA EXTRACTION': time_data_extraction,
        'PRE-RUN': time_prerun,
        'SCCs EXTRACTION': time_sccs_extraction,
        'SEED SEARCH': time_seed_search,
        'RENDERING': time_rendering,
        'TOTAL': time_data_extraction + time_rendering + time_seed_search,
    }
    namewidth = max(map(len, timers))
    for name, value in timers.items():
        value = value if isinstance(value, str) else f'{round(value, 3)}s'
        print(f'TIME {name.center(namewidth)}: {value}')
