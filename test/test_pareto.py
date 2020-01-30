"""Testing of the pareto exploration routine."""

from predator import graph_from_file, search_seeds, utils
from predator import quoted_data


def test_simple_reaction_chain():
    "A -> C -> T: C cannot be chosen, because it does not maximize the scope as much as A"
    graph = quoted_data('reactant(a,1). product(c,1). reaction(1). reactant(c,2). product(t,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('a')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True, greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='t', explore_pareto=True, greedy=True))
    seeds_sets = {frozenset('a')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True, greedy=False, verbose=True))
    assert seeds_sets == set(search_seeds(graph, targets='t', explore_pareto=True, greedy=False))

def test_long_reaction_chain():
    "L -> O -> N -> G -> E -> R: only L can be chosen, because it maximizes the scope"
    graph = quoted_data('reactant(l,1). product(o,1). reaction(1).'
                        'reactant(o,2). product(n,2). reaction(2).'
                        'reactant(n,3). product(g,3). reaction(3).'
                        'reactant(g,4). product(e,4). reaction(4).'
                        'reactant(e,5). product(r,5). reaction(5).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('l')}
    assert seeds_sets == set(search_seeds(graph, prerun=False))
    assert seeds_sets == set(search_seeds(graph, targets='r', explore_pareto=True, greedy=True, prerun=False))
    assert seeds_sets == set(search_seeds(graph, targets='r', explore_pareto=True, greedy=False, prerun=False))
    # Proof that the forbidden_seeds support has hard-to-manage side-effects:
    seeds_sets = {frozenset('o')}
    assert seeds_sets == set(search_seeds(graph, targets='r', forbidden_seeds='l', explore_pareto=True, greedy=True, prerun=False))
    # If prerun is forgot, this one will fail, because no model is generated:
    assert seeds_sets == set(search_seeds(graph, targets='r', forbidden_seeds='l', explore_pareto=True, greedy=False, prerun=True))
    # Because once arrived at SCC l, it is impossible to make use of it.
    # Prerun is fixing that.

def test_forbidden_intermediate_nodes():
    "Test the handling of middle SCC that are completely forbidden"
    graph = quoted_data('reactant((a;b;c),1). product(d,1). reaction(1).'
                        'reactant(d,2). product(e,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('e'), frozenset('abc'), frozenset('ae'), frozenset('be'), frozenset('ce')}
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', explore_pareto=True, greedy=True, prerun=False))
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', explore_pareto=True, greedy=True, prerun=False))
    seeds_sets = {frozenset('e')}
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', prerun=False, compute_optimal_solutions=True))
    seeds_sets = {frozenset('d'), frozenset('abc')}
    assert seeds_sets == set(search_seeds(graph, targets='e', explore_pareto=True, prerun=True))  # d as seeds covers more metabolites + less targets as seeds
    assert seeds_sets == set(search_seeds(graph, targets='e', explore_pareto=True, prerun=False))
    seeds_sets = {frozenset('e'), frozenset('abc')}
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', prerun=False))
    seeds_sets = {frozenset('e'), frozenset('abc')}
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', explore_pareto=True, greedy=True, prerun=False, verbose=True))
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', explore_pareto=True, prerun=True, verbose=True))
    # cf test_forbidden_sources() for explanations of why the previous line would fail.


def test_forbidden_intermediate_nodes_parallelized():
    graph = quoted_data('reactant((a;b;c),1). product(d,1). reaction(1).'
                        'reactant(d,2). product(e,2). reaction(2).'
                        'reactant((m;n;o),3). product(e,3). reaction(3).')
    utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('d'), frozenset('abc'), frozenset('mno')}
    assert seeds_sets == set(search_seeds(graph, targets='e', explore_pareto=True, prerun=True))
    seeds_sets = {frozenset('mno'), frozenset('abc'), frozenset('e')}
    assert seeds_sets == set(search_seeds(graph, targets='e', forbidden_seeds='d', prerun=False))


def test_forbidden_sources():
    "Test the handling of root SCC that are completely forbidden"
    graph = quoted_data('reactant(a,1). product(b,1). reaction(1).'
                        'reactant(b,2). product(a,2). reaction(2).'
                        'reactant(b,3). product(c,3). reaction(3).'
                        'reactant(c,4). product(d,4). reaction(4).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('a'), frozenset('b')}
    assert seeds_sets == set(search_seeds(graph))
    assert seeds_sets == set(search_seeds(graph, greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=True))
    seeds_sets = {frozenset('a'), frozenset('b'), frozenset('d')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='c', greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='c'))
    seeds_sets = {frozenset('c'), frozenset('d')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='ab', greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='ab'))
    seeds_sets = {frozenset('c')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abd', greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abd'))
    seeds_sets = {frozenset('d')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abc', greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abc'))
    seeds_sets = set()
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abcd', greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='abcd'))

    # But with forbidden seeds, it's not that obvious:
    # once arrived at SCC a (containing a and b), it is impossible to do anything.
    seeds_sets = {frozenset('a'), frozenset('b')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True))
    seeds_sets = {frozenset('a')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='b', explore_pareto=True, greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='b', explore_pareto=True, compute_optimal_solutions=True))
    # seeds_sets = {frozenset('c'), frozenset('a')}  # pareto exploration find optimals before compute_optimal_solutions post-process
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='b', explore_pareto=True, compute_optimal_solutions=False))
    seeds_sets = {frozenset('b')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='a', explore_pareto=True, greedy=True))
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='a', explore_pareto=True, compute_optimal_solutions=True))
    # seeds_sets = {frozenset('c'), frozenset('b')}  # pareto exploration find optimals before compute_optimal_solutions post-process
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='a', explore_pareto=True, compute_optimal_solutions=False))
    seeds_sets = {frozenset('c')}
    assert seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='ab', explore_pareto=True))

def test_simple_reaction_chain_with_double_root():
    "Same as previous, but with A+B -> C, allowing pareto front to choose either A+B or C."
    graph = quoted_data('reactant((a;b),1). product(c,1). reaction(1). reactant(c,2). product(t,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('c'), frozenset('ab')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True, greedy=True))
    seeds_sets = {frozenset('c'), frozenset('ab')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True, greedy=False))
    seeds_sets = {frozenset('c')}
    assert seeds_sets == set(search_seeds(graph, targets='t', forbidden_seeds='t', explore_pareto=True, greedy=False, verbose=True, compute_optimal_solutions=True))


def test_one_scc():
    """This show the local limits: whenever solutions are in the same SCC,
    only the very best are yielded. As a consequence, both pareto and without
    pareto loose some interesting (though non-optimal) solutions

    This problem has been fixed by implementing iterative pareto.
    (see last tests)

    """
    graph = quoted_data("""reactant(a,1). product(b,1). reaction(1).
                           reactant((b;c),2). product(d,2). reaction(2).
                           reactant(d,3). product(a,3). reaction(3).
                           reactant((d;e),4). product(f,4). reaction(4).
                           reactant(f,5). product(c,5). reaction(5).""")
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('d'), frozenset('de')}  # 'af' would have been welcome, but does not cover as much as 'de'
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=True, greedy=True, avoid_targets_as_seeds=True))
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=True, greedy=True, avoid_targets_as_seeds=False))  # invariant
    expected_seeds_sets = {frozenset('d')}  # this is… optimized
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', explore_pareto=False))

    # to get a real sens of diversity, we need to forbids the obvious solutions
    #  (here, just 'd', but it could be much more subtle).
    expected_seeds_sets = {frozenset('ac'), frozenset('af'), frozenset('bc'), frozenset('bf')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=False))
    expected_seeds_sets = {frozenset('af'), frozenset('bf'), frozenset('ace'), frozenset('aef'), frozenset('bce'), frozenset('bef')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, greedy=True))
    # hence, this test ask a question: do another metric for the greedy search
    #  in an SCC is needed ? Probably.
    # Something like a constraint providing *a maximal number/ratio of seed*.

    # This is approached using a mix of iterative and pareto solution:
    #  while using greedy and pareto.
    expected_seeds_sets = {frozenset('af'), frozenset('bf'), frozenset('ace'), frozenset('aef'), frozenset('bce'), frozenset('bef')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, greedy=False, verbose=True,
                                                   compute_optimal_solutions=False, filter_included_solutions=False))
    expected_seeds_sets = {frozenset('af'), frozenset('bf'), frozenset('ace'), frozenset('bce')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, greedy=False, verbose=True,
                                                   compute_optimal_solutions=False, filter_included_solutions=True))
    expected_seeds_sets = {frozenset('af'), frozenset('bf')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='d', forbidden_seeds='d', explore_pareto=True, greedy=False, verbose=True,
                                                   compute_optimal_solutions=True, filter_included_solutions=True))
