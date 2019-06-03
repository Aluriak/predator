"""Testing of the seeds search with targets.

"""
import pytest
from predator import graph_from_file, search_seeds, utils
from predator import quoted_data


def test_lattice_reaction():
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('bc'), frozenset('fab')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='ef', forbidden_seeds='e'))
    expected_seeds_sets = {frozenset('bc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='ef', forbidden_seeds='e', compute_optimal_solutions=True))
    assert expected_seeds_sets == set(search_seeds(graph, targets='f', forbidden_seeds='f'))
    expected_seeds_sets = {frozenset('ab'), frozenset('de'), frozenset('dbc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='de'))
    expected_seeds_sets = {frozenset('ab'), frozenset('de')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='de', compute_optimal_solutions=True))
    expected_seeds_sets = {frozenset('ab')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d'))
    expected_seeds_sets = {frozenset('ab'), frozenset('bc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='e'))
    assert {frozenset('abc')} == set(search_seeds(graph, targets='e', forbidden_seeds='e', enum_mode='union'))
    assert {frozenset('b')} == set(search_seeds(graph, targets='e', forbidden_seeds='e', enum_mode='intersection'))
    assert {frozenset('abc')} == set(search_seeds(graph, targets='def', forbidden_seeds='def'))
    assert {frozenset('ab')} == set(search_seeds(graph, targets='e', forbidden_seeds='ce'))


def test_impossible_case():
    "the code should return nothing when no solution is found, or raise ValueError when dubious question is asked"
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert not tuple(search_seeds(graph, forbidden_seeds='a'))  # impossible to activate a without a
    assert not tuple(search_seeds(graph, forbidden_seeds={'a', 'd'}, targets={'d'}))
    with pytest.raises(ValueError):
        list(search_seeds(graph, start_seeds='a', forbidden_seeds='a', targets='d'))


def test_hierarchical_case():
    graph = quoted_data('reactant((a;b),1). product(c,1). reaction(1).  reactant(a,2). product(e,2). reaction(2).  reactant((d;e),3). product(b,3). reaction(3).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ab'), frozenset('ad')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='c', targets='c'))
    expected_seeds_sets = {frozenset('c')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', compute_optimal_solutions=True))
    expected_seeds_sets = {frozenset('ab'), frozenset('ad'), frozenset('bc'), frozenset('cde')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc'))
    expected_seeds_sets = {frozenset('abcde')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', enum_mode='union'))
    expected_seeds_sets = {frozenset('ab'), frozenset('ad'), frozenset('bc')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', compute_optimal_solutions=True))
    expected_seeds_sets = {frozenset('abcde')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', compute_optimal_solutions=True, enum_mode='union'))
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', compute_optimal_solutions=False, enum_mode='union'))
    expected_seeds_sets = {frozenset()}
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', compute_optimal_solutions=True, enum_mode='intersection'))
    assert expected_seeds_sets == set(search_seeds(graph, targets='bc', compute_optimal_solutions=False, enum_mode='intersection'))


def test_non_optimal_local():
    graph = quoted_data('reactant((a;b),1). product(c,1). reaction(1).  reactant((c;d),2). product(e,2). reaction(2).  reactant(e,3). product(c,3). reaction(3).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('abd'), frozenset('cd')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e'))
    expected_seeds_sets = {frozenset('cd')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='d', forbidden_seeds='e', targets='e', compute_optimal_solutions=True))
    expected_seeds_sets = {frozenset('abd')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='ce', targets='e'))
    expected_seeds_sets = {frozenset('abd'), frozenset('adc')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='a', forbidden_seeds='e', targets='e'))
    expected_seeds_sets = {frozenset('abcd')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='a', forbidden_seeds='e', targets='e', enum_mode='union'))
    expected_seeds_sets = {frozenset('da')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='a', forbidden_seeds='e', targets='e', enum_mode='intersection'))
    expected_seeds_sets = {frozenset('d')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e', enum_mode='intersection', compute_optimal_solutions=False))
    expected_seeds_sets = {frozenset('dc')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e', enum_mode='intersection', compute_optimal_solutions=True))
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e', enum_mode='intersection', compute_optimal_solutions=True, greedy=True))


def test_combinations_of_deps():
    graph = quoted_data('reactant(a,1;c,2;b,3;d,3). reaction(1;2;3). product(b,1;d,2;e,3).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ac'), frozenset('ad'), frozenset('bc'), frozenset('bd')}
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e'))
    assert expected_seeds_sets == set(search_seeds(graph, forbidden_seeds='e', targets='e', compute_optimal_solutions=True))
    assert {frozenset('abcd')} == set(search_seeds(graph, forbidden_seeds='e', targets='e', enum_mode='union'))
    assert {frozenset()} == set(search_seeds(graph, forbidden_seeds='e', targets='e', enum_mode='intersection'))
