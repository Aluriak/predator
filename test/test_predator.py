"""Testing of the main routine.

"""

from predator import graph_from_file, search_seeds, utils


def test_simple_reaction():
    graph = 'reactant(a,r). product(p,r). reaction(r).'
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = [frozenset('a')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_reaction():
    graph = 'reactant((a;b),r). product(p,r). reaction(r).'
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = [frozenset('ab')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_with_feedback_reaction():
    graph = 'reactant((a;b),r). product(p,r). reaction(r). reactant(p,r2). product(a,r2). reaction(r2).'
    utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ab'), frozenset('pb')}
    assert expected_seeds_sets == set(search_seeds(graph))
