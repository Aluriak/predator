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
