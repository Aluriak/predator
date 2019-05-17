"""Utilitaries"""

import clyngor
from collections import defaultdict


BISEAU_VIZ = """
link(T,R) :- reactant(T,R).
link(R,P) :- product(P,R).
shape(R,rectangle) :- reaction(R).
obj_property(edge,arrowhead,vee).
"""
BISEAU_VIZ_NOREACTION = """
link(M,P) :- product(P,R) ; reactant(M,R).
obj_property(edge,arrowhead,vee).
"""


def solve(*args, **kwargs):
    "Wrapper around clyngor.solve"
    try:
        return clyngor.solve(*args, **kwargs)
    except FileNotFoundError as err:
        if 'clingo' in err.args[0]:
            print('ERROR binary file clingo is not accessible in the PATH.')
            exit(1)
        else:  raise err


def inverted_dag(dag:dict) -> [str]:
    """Return the same dag, but with edges reversed. Key None, if exists,
    must be associated with roots.

    """
    revdag = defaultdict(set)
    for pred, succs in dag.items():
        for succ in succs:
            if pred is None:
                revdag[succ]  # succ is now terminal
            else:
                revdag[succ].add(pred)  # build the reverse link
                revdag[pred]  # ensure that predecessor exists
            if succ not in dag:
                revdag[None].add(succ)
        if not succs:
            revdag[None].add(pred)
    return dict(revdag)

def get_terminal_nodes(dag:dict) -> [str]:
    """Yield terminal nodes of given dag, with doublons.
    Key None must be associated with roots
    """
    roots = list(dag[None])
    while roots:
        root = roots.pop()
        if root in dag and dag[root]:
            for sub in dag[root]:
                roots.append(sub)
        else:
            yield root


def remove_terminal(node, dag:dict, parents:set=None):
    "Modify dag in-place, removing given node as if it was a terminal"
    if not parents:
        parents = set()
        for pred, succs in dag.items():
            if node in succs:
                parents.add(pred)
    for parent in parents:
        dag[parent].remove(node)
    if node in dag:
        del dag[node]


def render_network(graph:str, outfile:str, with_reactions=True):
    "Render given graph using biseau, as a single image saved with given filename."
    import biseau
    viz = BISEAU_VIZ if with_reactions else BISEAU_VIZ_NOREACTION
    biseau.compile_to_single_image(graph + viz, outfile=outfile)
    return outfile


def quoted(string:str) -> str:
    r"""Return string, double quoted

    >>> quoted('"a').replace('\\', '$')
    '"$"a"'
    >>> quoted('"a b"').replace('\\', '$')
    '"a b"'
    >>> quoted('a b').replace('\\', '$')
    '"a b"'
    >>> quoted('a\\"').replace('\\', '$')
    '"a$""'
    >>> quoted('a"').replace('\\', '$')
    '"a$""'
    >>> quoted('\\"a"').replace('\\', '$')
    '"$"a$""'
    >>> quoted('"').replace('\\', '$')
    '"$""'

    """
    if len(string) > 1 and string[0] == '"' and string[-2] != '\\' and string[-1] == '"':
        return string
    else:
        return '"' + string.replace('\\"', '"').replace('"', '\\"') + '"'


def unquoted(string:str) -> str:
    r"""Remove surrounding double quotes if they are acting as such

    >>> unquoted('"a').replace('\\', '$')
    '$"a'
    >>> unquoted('"a b"')
    'a b'
    >>> unquoted('b"').replace('\\', '$')
    'b$"'
    >>> unquoted('"b\\"').replace('\\', '$')
    '$"b$"'

    """
    if string[0] == '"' and string[-2] != '\\' and string[-1] == '"':
        return string[1:-1]
    else:
        return string.replace('\\"', '"').replace('"', '\\"')


def quoted_data(asp:str) -> str:
    "Return the same atoms as found in given asp code, but with all arguments quoted"
    def gen():
        for model in clyngor.solve(inline=asp):
            for pred, args in model:
                yield f'{pred}(' + ','.join(quoted(str(arg)) for arg in args) + ').'
    return ' '.join(gen())
