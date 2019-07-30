"""Testing of graph functions"""


from predator import predator
from predator.graph import nx_from_asp, sccs_dag_from_nxdigraph


def test_sccs_dag_simple():
    dag = sccs_dag_from_nxdigraph(nx_from_asp('edge(a,b).edge(b,c).'))
    assert dag == {None: {'a'}, 'a': {'b'}, 'b': {'c'}}


def test_sccs_dag_complex():
    dag = sccs_dag_from_nxdigraph(nx_from_asp('edge(a,b).edge(b,c).edge(c,a).'
                                              'edge(d,e).edge(e,d).'
                                              'edge(c,f).edge(d,f).'
                                              'edge(f,g). edge(g,f).'
                                              'edge(f,h).'
                                              'edge(h,i).edge(i,h).'))
    assert dag == {None: {'a', 'd'}, 'a': {'f'}, 'd': {'f'}, 'f': {'h'}}


def test_sccs_dag_cycle_4():
    dag = sccs_dag_from_nxdigraph(nx_from_asp('edge(a,b).edge(b,c).edge(c,d).edge(d,a).'))
    assert dag == {None: {'a'}}


def test_sccs_dag_cycle_4_alt():
    sccs, dag = predator.compute_sccs('reactant(a,r1).reaction(r1).product(b,r1).'
                                      'reactant(b,r2).reaction(r2).product(c,r2).'
                                      'reactant(c,r3).reaction(r3).product(d,r3).'
                                      'reactant(d,r4).reaction(r4).product(a,r4).'
                                      , verbose=True)
    assert dag == {None: {'a'}}
