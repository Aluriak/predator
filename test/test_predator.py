"""Testing of the main routine.

"""

from predator import graph_from_file, search_seeds


def test_simple_reaction():
    graph = 'reactant(a,r). product(p,r). reaction(r).'
    seeds_sets = [frozenset('a')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_reaction():
    graph = 'reactant((a;b),r). product(p,r). reaction(r).'
    seeds_sets = [frozenset('ab')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_with_feedback_reaction():
    graph = 'reactant((a;b),r). product(p,r). reaction(r). reactant(p,r2). product(a,r2). reaction(r2).'
    seeds_sets = [frozenset('ab'), frozenset('ap')]
    assert seeds_sets == list(search_seeds(graph))
