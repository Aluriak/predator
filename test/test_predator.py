"""Testing of the main routine.

"""

from predator import graph_from_file, search_seeds, utils, EnumMode
from predator import quoted_data


def test_basic_EnumMode_API():
    for obj in EnumMode:
        assert getattr(EnumMode, obj.value.title()) is EnumMode(obj.value)
    assert EnumMode('enumeration').clingo_option == ''
    assert EnumMode('union').clingo_option == '--enum-mode brave'
    assert EnumMode('intersection').clingo_option == '--enum-mode cautious'

def test_simple_reaction():
    graph = quoted_data('reactant(a,r). product(p,r). reaction(r).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('a')}
    assert seeds_sets == set(search_seeds(graph))
    seeds_sets = {frozenset('a'), frozenset('p')}
    assert seeds_sets == set(search_seeds(graph, targets='p'))


def test_double_reactant_reaction():
    graph = quoted_data('reactant((a;b),r). product(p,r). reaction(r).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    seeds_sets = {frozenset('ab')}
    assert seeds_sets == set(search_seeds(graph))
    seeds_sets = {frozenset('ab'), frozenset('p')}
    assert seeds_sets == set(search_seeds(graph, targets='p', verbose=True))


def test_loss_because_of_optimal_locality():
    DATA = """
        reactant(X,R):- reaction(X,R,_). product(X,R) :- reaction(_,R,X). reaction(R) :- reaction(_,R,_).
        reaction((a;b;c),r1,(d;e;f)).
        reaction((d;e;f),r1r,(a;b;c)).
        reaction((d;e;f),r2,(g;h)).
        reaction(g,r3,i).
        reaction(i,r4,j).
        reaction(j,r5,g).
        reaction(h,r6,k).
        reaction(k,r7,l).
        reaction(l,r8,h).
        reaction((j;k),r9,m).
    """
    # two versions of the graph: with or without z
    graph = quoted_data(DATA)
    zgraph = quoted_data(DATA + 'reaction(z,r0,(a;b;c)).')
    # utils.render_network(zgraph+'dot_property("z",style,dashed).', 'todel.png')  # uncomment to help debugging

    # fully expected results
    seeds_sets = {frozenset('abc'), frozenset('def')}
    assert seeds_sets == set(search_seeds(graph))
    seeds_sets = {frozenset('z')}
    assert seeds_sets == set(search_seeds(zgraph))

    seeds_sets = {
        # frozenset('z'),  # not in the zgraph
        frozenset('abc'), frozenset('def'),
        frozenset('hi'), frozenset('hj'), frozenset('hg'),
        frozenset('ki'), frozenset('kj'), frozenset('kg'),
        frozenset('li'), frozenset('lj'), frozenset('lg')
    }
    assert seeds_sets == set(search_seeds(graph, targets='m', forbidden_seeds='m', verbose=True))
    seeds_sets = {
        frozenset('z'),
        # frozenset('abc'), frozenset('def'),  # for a subtil reason, these does not belong to the solutions
        frozenset('hi'), frozenset('hj'), frozenset('hg'),
        frozenset('ki'), frozenset('kj'), frozenset('kg'),
        frozenset('li'), frozenset('lj'), frozenset('lg')
    }
    assert seeds_sets == set(search_seeds(zgraph, targets='m', forbidden_seeds='m', verbose=True))
    # NOTE: In that case, abc and def are no longer possible seeds.
    # This is quite unnatural, but it's logical with current implementation:
    #  when clingo is studying the SCC 'a', containing abcdef, and child of SCC containing z,
    #  it faces two main choices: make abc (or def) seeds, or ask z to be activated.
    #  3 seeds vs 1 input -> only the z solution is kept.

    seeds_sets = {frozenset('hi'), frozenset('hj'), frozenset('hg'),
                  frozenset('ki'), frozenset('kj'), frozenset('kg'),
                  frozenset('li'), frozenset('lj'), frozenset('lg')}
    assert seeds_sets == set(search_seeds(graph, targets='m', forbidden_seeds='m', verbose=True, compute_optimal_solutions=True))
    seeds_sets = {frozenset('z')}
    assert seeds_sets == set(search_seeds(zgraph, targets='m', forbidden_seeds='m', verbose=True, compute_optimal_solutions=True))


def test_double_reactant_with_feedback_reaction():
    graph = quoted_data('reactant((a;b),r). product(p,r). reaction(r). reactant(p,r2). product(a,r2). reaction(r2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('ab'), frozenset('pb')}
    assert expected_seeds_sets == set(search_seeds(graph))
    expected_seeds_sets = {frozenset('abp')}
    assert expected_seeds_sets == set(search_seeds(graph, enum_mode='union'))
    expected_seeds_sets = {frozenset('b')}
    assert expected_seeds_sets == set(search_seeds(graph, enum_mode='intersection'))


def test_local_minimization():
    "Is it sure that minimizing the amount of seeds/ingoings is the right way ?"
    graph = quoted_data('reactant(a,1). product((b;c),1). reaction(1).  reactant((b;c),2). product(a,2). reaction(2).  reactant(d,3). product(b,3). reaction(3).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    expected_seeds_sets = {frozenset('a')}
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', forbidden_seeds='c'))
    assert expected_seeds_sets == set(search_seeds(graph, targets='c', forbidden_seeds='c', compute_optimal_solutions=True))


def test_loop():
    graph = quoted_data('reactant(a,1). product(b,1). reaction(1).  reactant(b,2). product(a,2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert {frozenset('a'), frozenset('b')} == set(search_seeds(graph))
    assert {frozenset('ab')} == set(search_seeds(graph, enum_mode='union'))
    assert {frozenset()} == set(search_seeds(graph, enum_mode='intersection'))


def test_and_or():
    graph = quoted_data('reactant("b","r2"). reactant("c","r3"). reactant("a","r1"). reactant("b","r1"). reactant("c","r1"). product("c","r2"). product("b","r3"). product("d","r1"). reaction("r2"). reaction("r3"). reaction("r1").')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert {frozenset('ab'), frozenset('ac')} == set(search_seeds(graph))
    assert {frozenset('ab'), frozenset('ac'), frozenset('d')} == set(search_seeds(graph, targets='d'))
    assert {frozenset('ab'), frozenset('ac')} == set(search_seeds(graph, targets='d', forbidden_seeds='d'))
    assert {frozenset('abc')} == set(search_seeds(graph, enum_mode='union'))
    assert {frozenset('a')} == set(search_seeds(graph, enum_mode='intersection'))


def test_with_seeds_and_forbidden():
    graph = quoted_data('reactant((a;b),1). product((d;e),1). reaction(1).  reactant((b;c),2). product((e;f),2). reaction(2).')
    # utils.render_network(graph, 'todel.png')  # uncomment to help debugging
    assert not set(search_seeds(graph, start_seeds='a', forbidden_seeds='c'))
    expected_seeds_sets = {frozenset('abc')}
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='abc', forbidden_seeds='def'))
    assert expected_seeds_sets == set(search_seeds(graph, start_seeds='c', forbidden_seeds='d', targets='d'))
