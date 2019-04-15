"""Testing of the seeds search with targets.

"""
from predator import graph_from_file, search_seeds, utils
from predator import quoted_data


def test_lattice_reaction():
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('bc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='ef', forbidden_seeds='e'))
    assert expected_seeds_sets == set(search_seeds(graph, targets='f', forbidden_seeds='f'))
    expected_seeds_sets = {frozenset('ab'), frozenset('de')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='de'))
    expected_seeds_sets = {frozenset('ab')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d'))
    expected_seeds_sets = {frozenset('ab'), frozenset('bc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='e'))
    expected_seeds_sets = {frozenset('abc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='def', forbidden_seeds='def'))
    expected_seeds_sets = {frozenset('ab')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='ce'))


def test_impossible_case():
    "the code should return nothing when no solution is found"
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert not list(search_seeds(graph, forbidden_seeds={'a', 'd'}, targets={'d'}))
