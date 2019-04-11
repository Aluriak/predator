"""Utilitaries"""


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


def render_network(graph:str, outfile:str, with_reactions=False):
    "Render given graph using biseau, as a single image saved with given filename."
    import biseau
    viz = BISEAU_VIZ if with_reactions else BISEAU_VIZ_NOREACTION
    biseau.compile_to_single_image(graph + viz, outfile=outfile)
    return outfile
