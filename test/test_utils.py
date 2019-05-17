
from predator.utils import inverted_dag, get_terminal_nodes, remove_terminal


def test_inverted_dag():
    dag = {1: {3, 4}, 2: {4}, 3: {5}, 4: {5}, 5: set(), None: {1, 2}}
    found = inverted_dag(dag)
    expected = {1: set(), 2: set(), 3: {1}, 4: {1, 2}, 5: {3, 4}, None: {5}}
    assert found == expected

def test_inverted_dag_from_bug1():
    dag = {'c': {'e', 'f'}, 'b': {'e', 'f', 'd'}, 'a': {'e', 'd'}, None: {'b', 'a', 'c'}}
    found = inverted_dag(dag)
    expected = {'e': {'c', 'b', 'a'}, 'f': {'b', 'c'}, 'd': {'a', 'b'},
                'a': set(), 'b': set(), 'c': set(), None: {'d', 'e', 'f'}}
    assert found == expected

def test_inverted_dag_from_bug1_without_none():
    # same as previous, but with None key removed
    dag = {'c': {'e', 'f'}, 'b': {'e', 'f', 'd'}, 'a': {'e', 'd'}}
    found = inverted_dag(dag)
    expected = {'e': {'c', 'b', 'a'}, 'f': {'b', 'c'}, 'd': {'a', 'b'},
                'a': set(), 'b': set(), 'c': set(), None: {'d', 'e', 'f'}}
    assert found == expected


def test_get_terminal_nodes():
    dag = {1: {3, 4}, 2: {4}, 3: {5}, 4: {5}, 5: {}, None: {1, 2}}
    found = set(get_terminal_nodes(dag))
    expected = {5}
    assert found == expected


def test_remove_terminal():
    dag = {1: {3, 4}, 2: {4}, 3: {5}, 4: {5}, 5: set(), None: {1, 2}}
    remove_terminal(5, dag)
    assert dag == {1: {3, 4}, 2: {4}, 3: set(), 4: set(), None: {1, 2}}


def test_remove_terminal_horizontal():
    dag = {'c': {'e', 'f'}, 'b': {'e', 'f', 'd'}, 'a': {'e', 'd'}, None: {'b', 'a', 'c'}}
    remove_terminal('d', dag)
    assert dag == {'c': {'e', 'f'}, 'b': {'e', 'f'}, 'a': {'e'}, None: {'b', 'a', 'c'}}
    remove_terminal('e', dag)
    assert dag == {'c': {'f'}, 'b': {'f'}, 'a': set(), None: {'b', 'a', 'c'}}
    remove_terminal('f', dag)
    assert dag == {'c': set(), 'b': set(), 'a': set(), None: {'b', 'a', 'c'}}
    remove_terminal('a', dag)
    assert dag == {'c': set(), 'b': set(), None: {'b', 'c'}}
    remove_terminal('b', dag)
    assert dag == {'c': set(), None: {'c'}}
    remove_terminal('c', dag)
    assert dag == {None: set()}
