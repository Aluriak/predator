"""Routines to search for seeds.

"""
import math
import clyngor
import networkx as nx
import itertools
from . import sbml as sbml_module
from collections import defaultdict
from pkg_resources import resource_filename


ASP_SRC_ENUM_CC = resource_filename(__name__, 'asp/enum-cc.lp')
ASP_SRC_SIMPLE_SEED_SOLVING = resource_filename(__name__, 'asp/simple-seed-solving.lp')
ASP_SRC_GREEDY_TARGET_SEED_SOLVING = resource_filename(__name__, 'asp/greedy-target-seed-solving.lp')
import os
assert os.path.exists(ASP_SRC_ENUM_CC), ASP_SRC_ENUM_CC
assert os.path.exists(ASP_SRC_SIMPLE_SEED_SOLVING), ASP_SRC_SIMPLE_SEED_SOLVING


def opt_models_from_clyngor_answers(answers:iter, *, smaller_is_best:bool=True):
    "Return tuple of optimal models found by clingor.solve, when option '--opt-mode=optN' is given"
    best_opt, models = math.inf, []
    if not smaller_is_best:  best_opt *= -1
    for model, opt in answers.with_optimization:
        # print('OPT, MODEL:', opt[0], model)
        if (opt[0] < best_opt) if smaller_is_best else (opt[0] > best_opt):
            best_opt, models = opt[0], []  # model will be given again as last model, so no need to include it twice
        else:
            models.append(model)
    return tuple(models)


def search_seeds(graph_data:str, start_seeds:iter=(), forbidden_seeds:iter=(),
                 targets:set=(), graph_filename:str=None) -> [str]:
    if not targets:  # no target, just activate everything
        func = search_seeds_activate_all
    else:
        func = search_seeds_activate_targets_greedy  # TODO: use the non-greedy solution
        # func = search_seeds_activate_targets_iterative  # TODO: implement that
    yield from func(graph_data, start_seeds, forbidden_seeds, targets, graph_filename)


def search_seeds_activate_targets_greedy(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None) -> [str]:
    """Yield the set of seeds for each found solution.

    This implements the activation of targets: find the minimum sets
    of seeds that activate all targets.
    This is a greedy implementation. Do not expect it to work on a large dataset.

    """
    if not targets:
        raise ValueError("search_seeds_activate_targets_greedy() requires targets. Use another function, search_seeds(), or provide targets.")
    start_seeds_repr = ' '.join(f'seed("{s}").' for s in start_seeds)
    targets_repr = ' '.join(f'target("{t}").' for t in targets)
    forb_repr = ' '.join(f'forbidden("{s}").' for s in forbidden_seeds)
    data_repr = graph_data + start_seeds_repr + forb_repr + targets_repr
    models = clyngor.solve(ASP_SRC_GREEDY_TARGET_SEED_SOLVING, inline=data_repr, options='--opt-mode=optN')
    models = opt_models_from_clyngor_answers(models.by_predicate.discard_quotes)
    for model in models:
        seeds = frozenset(args[0] for args in model['seed'] if len(args) == 1)
        yield seeds


def search_seeds_activate_all(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None) -> [str]:
    """Yield the set of seeds for each found solution.

    This implements the most simple solution: no targets, find the minimum sets
    of seeds that activate everything.

    """
    if targets:
        raise ValueError("search_seeds_activate_all() does not handle targets. Use another method that supports it, or search_seeds().")
    start_seeds_repr = ' '.join(f'seed({s}).' for s in start_seeds)
    forbidden_repr = ' '.join(f'forbidden({s}).' for s in forbidden_seeds)
    assert isinstance(graph_data, str), graph_data
    sccs, scc_dag = compute_sccs(graph_data, graph_filename=graph_filename)
    roots = scc_dag[None]
    scc_seeds = {}  # scc name: list of {optimal seeds}
    for scc_name, nodes in sccs.items():
        scc_repr = ' '.join(f'scc({scc_name},{node}).' for node in nodes)
        scc_data = f'current_scc({scc_name}). {scc_repr} {start_seeds_repr} {forbidden_repr}'
        # print('DATA:\n    ' + scc_data)
        # print('    ' + graph_data)
        # print()
        models = clyngor.solve(ASP_SRC_SIMPLE_SEED_SOLVING, inline=graph_data + scc_data, options='--opt-mode=optN').discard_quotes
        models = opt_models_from_clyngor_answers(models.by_predicate)
        scc_seeds[scc_name] = tuple(frozenset(args[0] for args in model.get('seed', ())) for model in models)
        # print('OUTPUT SEEDS:', scc_seeds[scc_name])
    # generate all possibilities
    seed_combinations = tuple(itertools.product(*scc_seeds.values()))
    for idx, seeds_sets in enumerate(seed_combinations, start=1):
        # seeds_set is a combination of possible set of seeds for each SCC
        yield frozenset.union(*seeds_sets)


def compute_sccs(graph_data:str, graph_filename:str=None) -> [{str}]:
    """Return the DAG and nodes of Strongly Connected Components found in given graph.

    graph_data -- ASP encoding of the metabolic network.
    graph_filename -- filename of SBML encoded graph, used to speed up SCC mining.
    return -- (SCC, DAG), with DAG a SCC id -> {SCC id sucessor} mapping,
              and SCC a SCC id -> {node} mapping.

    If no graph_filename given, ASP will be used to infer them from graph_data.

    """
    if not graph_filename or True:  # the networkx method is not yet implemented
        models = clyngor.solve(ASP_SRC_ENUM_CC, inline=graph_data)
        for model in models.by_predicate:
            print('SCC MODEL:', model)
            roots = frozenset(args[0] for args in model.get('noinput', ()) if len(args) == 1)
            sccs = defaultdict(set)  # SCC identifier: nodes in SCC
            for scc_name, node in model.get('scc', ()):
                sccs[scc_name].add(node)
            scc_dag = defaultdict(set)  # SCC identifier: successor SCCs
            for scc_name, scc_succ in model.get('sccedge', ()):
                scc_dag[scc_name].add(scc_succ)
            scc_dag[None] = roots
            sccs = dict(sccs)
            scc_dag = dict(scc_dag)
        print('SCC FINAL:', sccs)
        print()
        print('         :', scc_dag)
        print()
        return sccs, scc_dag
    else:  # use graph_filename
        raise NotImplementedError("Networkx based SCC extraction is not yet implemented")
        digraph = sbml_module.read_SBML_network_as_simple_graph(graph_filename)
        sccs = {min(nodes): frozenset(nodes) for nodes in nx.strongly_connected_components(digraph)}
        return sccs, scc_dag
