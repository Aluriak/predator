"""Utilitaries"""

import clyngor


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

    """
    if string[0] == '"' and string[-2] != '\\' and string[-1] == '"':
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
