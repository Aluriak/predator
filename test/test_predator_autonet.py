"""Testing of the main routine using the automatic discovery of example networks
and their .test (if available).

"""

import os
import glob
from predator import graph_from_file, search_seeds, utils


def seeds_from_file(fname:str) -> [str]:
    """Yield sets of seeds found in given file

    In the file, an empty line describes the beginning of a new set of seed
    (solution expected).
    A line beginning by `#` is ignored.

    """
    with open(fname) as fd:
        seeds = set()  # current set of current of seeds, next to be yield
        for line in map(str.strip, fd):
            if not line:  # a new set is beginning
                yield frozenset(seeds)
                seeds = set()
            elif line.startswith('#'):  # comment
                pass
            else:  # add to current set
                seeds.add(line)
    if seeds:  # last set
        yield frozenset(seeds)


def template_test_network_example(graph_file, seeds_file):
    expected_seeds_sets = frozenset(seeds_from_file(seeds_file))
    graph = graph_from_file(graph_file)
    def test_network_example():
        seeds_sets = set(search_seeds(graph))
        # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
        assert seeds_sets == expected_seeds_sets
    return test_network_example


# generation of test functions
for testfile in glob.glob('networks/*.test'):
    basename = os.path.splitext(testfile)[0]
    name = os.path.split(basename)[1]
    for ext in 'xml,sbml,lp'.split(','):
        graphfile = basename + '.' + ext
        # assert False, graphfile
        if os.path.exists(graphfile):
            globals()[f'test_network_{name}_{ext}'] = template_test_network_example(graphfile, testfile)
