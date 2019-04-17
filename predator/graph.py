"""Routines for graph data extractions.

"""


import os
import clyngor
import networkx as nx
from . import sbml, utils
from collections import defaultdict


def graph_from_file(fname:str, **kwargs) -> str:
    """Return ASP data encoding graph in given file."""
    ext = os.path.splitext(fname)[1]
    if ext in {'.xml', '.sbml'}:
        return '\n'.join(sbml.readSBMLnetwork(fname, **kwargs))
    elif ext == '.lp':
        with open(fname) as fd:
            return fd.read()
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")


def nxgraph_from_file(fname:str) -> str:
    """Return nx.Digraph encoding graph in given file."""
    ext = os.path.splitext(fname)[1]
    if ext in {'.xml', '.sbml'}:
        return sbml.read_SBML_network_as_simple_graph(fname)
    elif ext == '.lp':
        return nx_from_asp(graph_from_file(fname))
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")


def nx_from_asp(asp:str, directed:bool=True):
    "Return a nx.Digraph describing the graph given in ASP format"
    graph = (nx.Digraph if directed else nx.Graph)()
    models = clyngor.solve(inline=asp).by_predicate.discard_quotes
    for model in models:
        for args in model.get('node', ()):
            if len(args) == 1:
                graph.add_node(*args)
        for args in model.get('edge', ()):
            if len(args) == 2:
                graph.add_edge(*args)
    return graph


def print_info(fname:str, render_in:str=None, render_in_without_reactions:str=None):
    graph = nxgraph_from_file(fname)
    size = defaultdict(int)  # nb node: nb of SCC having this number of node
    for nodes in nx.strongly_connected_components(graph):
        nodes = frozenset(nodes)
        size[len(nodes)] += 1
    nb_node_computed = sum(nbnode * count for nbnode, count in size.items())
    assert nb_node_computed == graph.number_of_nodes(), (nb_node_computed, graph.number_of_nodes())
    print(f"- {graph.number_of_nodes()} nodes")
    print(f"- {graph.number_of_edges()} edges")
    print(f"- {sum(size.values())} Strongly Connected Components")
    nb_nontrivial_scc = sum(count for nbnode, count in size.items() if nbnode > 1)
    print(f"- {nb_nontrivial_scc} Strongly Connected Components having more than 1 node")
    headers = '#node', '#scc'
    colwidth = len(headers[0])
    print(' | '.join(header.center(len(header)) for header in headers))
    for nbnode in sorted(tuple(size.keys())):
        print(str(nbnode).rjust(len(headers[0])), '|', str(size[nbnode]).rjust(len(headers[1])))

    if render_in or render_in_without_reactions:
        asp_graph = graph_from_file(fname)
    if render_in:
        print('Input graph rendered in', utils.render_network(asp_graph, render_in, with_reactions=True))
    if render_in_without_reactions:
        print('Input graph rendered in', utils.render_network(asp_graph, render_in_without_reactions, with_reactions=False))
