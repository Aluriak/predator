"""Some test over the sbml support, using toy examples"""

import predator

def test_toy1():
    "Usage of only graph_filename as parameter, enabling the use of networkx to search for SCCs"
    seeds_sets = set(predator.search_seeds(graph_filename='./networks/toy_1.sbml'))
    assert seeds_sets == set(map(frozenset, ({'s', 'b'}, {'s', 'c'}, {'s', 'd'}, {'s', 'e'})))
