"""Testing of the seeds search with multiple targets.

"""

from predator import search_seeds, quoted_data, utils

def test_simple_division():
    graph = quoted_data('reactant(a,1). product(d,1). reaction(1).'
                        'reactant(c,2). product(e,2). reaction(2).'
                        'reactant(b,3). product((d;e),3). reaction(3).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ac'), frozenset('b')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='de', forbidden_seeds='de'))
    expected_seeds_sets = {frozenset('b')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='de', forbidden_seeds='de', compute_optimal_solutions=True))

def test_infinite_loop_of_hypothesis():
    graph = quoted_data('reactant((a;c),1). product(b,1). reaction(1).'
                        'reactant((b;d),2). product(c,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging

    expected_seeds_sets = {frozenset('bd')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', explore_pareto=True, forbidden_seeds='c'))

    expected_seeds_sets = {frozenset('bd'), frozenset('c')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', explore_pareto=True))

    expected_seeds_sets = {frozenset('bd'), frozenset('c')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', explore_pareto=True))

def test_double_yielder():
    graph = quoted_data('reactant(a,1). product(c,1). reaction(1).'
                        'reactant(a,2). product(d,2). reaction(2).'
                        'reactant(b,3). product(c,3). reaction(3).'
                        'reactant(c,r3). product(b,r3). reaction(r3).'
                        'reactant(d,4). product(e,4). reaction(4).'
                        'reactant(e,r4). product(d,r4). reaction(r4).'
                        'reactant(e,5). product(g,5). reaction(5).'
                        'reactant(b,6). product(f,6). reaction(6).')
    graph_7 = graph + quoted_data('reactant(d,7). product(c,7). reaction(7).')
    # utils.render_network(graph_7 + 'dot_property("7",style,dashed).', 'todel.png')  # uncomment to help debugging

    expected_seeds_sets = {frozenset('a'), frozenset('cd'), frozenset('ce'), frozenset('bd'), frozenset('be')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='fg', forbidden_seeds='fg'))

    expected_seeds_sets = {frozenset('a'), frozenset('d'), frozenset('e')}
    assert expected_seeds_sets == set(search_seeds(graph_7, targets='fg', forbidden_seeds='fg'))

    # optimal search
    expected_seeds_sets = {frozenset('a')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='fg', forbidden_seeds='de', compute_optimal_solutions=True))

    expected_seeds_sets = {frozenset('a'), frozenset('d'), frozenset('e')}
    assert expected_seeds_sets == set(search_seeds(graph_7, targets='fg', forbidden_seeds='fg', compute_optimal_solutions=True))
