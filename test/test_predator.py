"""Testing of the main routine.

"""

from predator import graph_from_file, search_seeds, utils
from predator import quoted_data


def test_simple_reaction():
    graph = quoted_data('reactant(a,r). product(p,r). reaction(r).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = [frozenset('a')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_reaction():
    graph = quoted_data('reactant((a;b),r). product(p,r). reaction(r).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = [frozenset('ab')]
    assert seeds_sets == list(search_seeds(graph))


def test_double_reactant_with_feedback_reaction():
    graph = quoted_data('reactant((a;b),r). product(p,r). reaction(r). reactant(p,r2). product(a,r2). reaction(r2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ab'), frozenset('pb')}
    assert expected_seeds_sets == set(search_seeds(graph))


def test_loop():
    graph = quoted_data('reactant(a,1). product(b,1). reaction(1).  reactant(b,2). product(a,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('a'), frozenset('b')}
    assert expected_seeds_sets == set(search_seeds(graph))


def test_with_seeds_and_forbidden():
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert not set(search_seeds(graph, start_seeds='a', forbidden_seeds='c'))
    expected_seeds_sets = {frozenset('abc')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='abc', forbidden_seeds='def'))
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='c', forbidden_seeds='d', targets='d'))
