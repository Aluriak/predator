"""Implementation of the prerun function, to run on a graph before
feeding it to search_seeds.

"""
import tempfile
import clyngor
import networkx as nx
import itertools
from . import graph as graph_module
from .utils import quoted, unquoted, solve


def on_graph(graph:str, forbidden_seeds:iter=(), targets:set=(), verbose=False) -> (nx.DiGraph, str, frozenset, frozenset, str):
    """Return the new graph object, the new graph filename, the new set of targets, the set of unreachables targets,
    and the ASP/biseau encoding to add to the input graph to represent the dead parts of the graph"""

    # print('GRAPH:', graph)
    nxgraph = graph_module.nx_from_asp(graph)
    tosolve = graph + r"""
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).
activated(M) :- metabolite(M) ; not forbidden(M).
activated(R) :- reaction(R) ; activated(T): reactant(T,R).
activated(M) :- metabolite(M) ; product(M,R) ; activated(R).
dead(R) :-  reaction(R)  ; not activated(R).
dead(M) :- metabolite(M) ; not activated(M).
% Ensure that activated nodes that are alone still are included in the new graph.
singleton_metabolite(M) :- metabolite(M) ; dead(R): product(M,R) ; dead(R): reactant(M,R) ; not dead(M).
#show.
#show reactant(M,R): reactant(M,R), not dead(R) ; not dead(M).
#show product(M,R): product(M,R), not dead(R), not dead(M).
#show reaction(R): reaction(R), not dead(R) ; reactant(M,R) ; not dead(M).
#show target(T): target(T), not dead(T).
#show metabolite(M): singleton_metabolite(M).
#show dead/1.
""" + (
    '\n' + '\n'.join(f'forbidden("{fs}").' for fs in forbidden_seeds) +
    '\n' + '\n'.join(f'target("{target}").' for target in targets)
    )

    # tosolve = ('a(p).')
    models = solve(inline=tosolve, delete_tempfile=False).by_predicate
    # print('\nTOSOLVE:\n' + tosolve + '\n')
    # print('CMD:', '"' + models.command + '"')  # WOOT ? TODO
    for model in models:
        graph = get_atoms(model, ('reactant', 'product', 'reaction', 'metabolite'))
        if not graph:  # everything was dead
            graph = '%* nothing in the graph was reachable *%\n'
            if verbose:
                print('Nothing in the graph is reachable: all the graph was wiped.')
        # print('NEW GRAPH:', graph)
        dead_parts = get_atoms(model, ('dead',))
        # print('DEAD parts:', dead_parts)
        biseau_rules = '\n' + dead_parts + '\ncolor(X,black):- dead(X).\ndot_property(X,fontcolor,white):- dead(X).'
        # print('BISEAU RULES:', biseau_rules)
        new_targets = {unquoted(t) for t, in model.get('target', ())}
        unfullfilled_targets = set(map(quoted, set(targets) - new_targets))
        # print('NEW TARGETS:', new_targets)
        # print('    TARGETS:', targets)
        # print('UNFULLFILLS:', unfullfilled_targets)
        break
    else:
        print('Unexpected multiple models for prerun.on_graph solving')
        exit(1)
    with tempfile.NamedTemporaryFile('w', suffix='.lp', delete=False) as fd:
        fd.write(graph)
        new_graph_filename = fd.name
        if verbose:
            print(f'Pre-processed graph saved in file {new_graph_filename}.')
    return graph, new_graph_filename, new_targets, unfullfilled_targets, biseau_rules


def get_atoms(model:dict, wanted_atoms:[str]) -> str:
    "Return atoms found in given model corresponding to given signatures"
    wanted_atoms = set(wanted_atoms)
    model = {pred: args for pred, args in model.items() if pred in wanted_atoms}
    return clyngor.utils.answer_set_to_str(model, atom_end='.')
