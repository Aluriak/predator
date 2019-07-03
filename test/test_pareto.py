"""Testing of the pareto exploration routine."""

from predator import graph_from_file, search_seeds, utils
from predator import quoted_data


def test_simple_reaction_chain():
    "A -> C -> T: C cannot be chosen, because it does not maximize the scope as much as A"
    graph = quoted_data('reactant(a,1). product(c,1). reaction(1). reactant(c,2). product(t,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('a')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True))


def test_simple_reaction_chain_with_double_root():
    "Same as previous, but with A+B -> C, allowing pareto front to choose either A+B or C."
    graph = quoted_data('reactant((a;b),1). product(c,1). reaction(1). reactant(c,2). product(t,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('c'), frozenset('ab')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True))


def test_one_scc():
    """This show the local limits: whenever solutions are in the same SCC,
    only the very best are yielded. As a consequence, both pareto and without
    pareto loose some interesting (though non-optimal) solutions"""
    graph = quoted_data("""reactant(a,1). product(b,1). reaction(1).
                           reactant((b;c),2). product(d,2). reaction(2).
                           reactant(d,3). product(a,3). reaction(3).
                           reactant((d;e),4). product(f,4). reaction(4).
                           reactant(f,5). product(c,5). reaction(5).""")
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('d'), frozenset('de')}  # 'af' would have been welcome, but does not cover as much as 'de'
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=True, verbose=True, avoid_targets_as_seeds=True))
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=True, verbose=True, avoid_targets_as_seeds=False))  # invariant
    expected_seeds_sets = {frozenset('d')}  # this is… optimized
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=False, verbose=True))

    # to get a real sens of diversity, we need to forbids the obvious solutions
    #  (here, just 'd', but it could be much more subtle).
    expected_seeds_sets = {frozenset('ac'), frozenset('af'), frozenset('bc'), frozenset('bf')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=False, verbose=True))
    expected_seeds_sets = {frozenset('af'), frozenset('bf'), frozenset('ace'), frozenset('aef'), frozenset('bce'), frozenset('bef')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, verbose=True))
    # hence, this test ask a question: do another metric for the greedy search
    #  in an SCC is needed ? Probably.
    # Something like a constraint providing *a maximal number/ratio of seed*.
