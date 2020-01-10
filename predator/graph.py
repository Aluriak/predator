"""Routines for graph data extractions.

"""


import os
import itertools
import networkx as nx
from . import sbml, utils
from collections import defaultdict


def graph_from_file(fname:str, **kwargs) -> str:
    """Return ASP data encoding graph in given file."""
    ext = os.path.splitext(fname)[1]
    if ext in {'.xml', '.sbml'}:
        return '\n'.join(sbml.read_SBML_network(fname, **kwargs))
    elif ext == '.lp':
        with open(fname) as fd:
            return fd.read()
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")


def nxgraph_from_file(fname:str, asp_data:str=None, **kwargs) -> str:
    """Return nx.Digraph encoding graph in given file."""
    ext = os.path.splitext(fname)[1]
    if ext in {'.xml', '.sbml'}:
        return sbml.read_SBML_network_as_simple_graph(fname, **kwargs)
    elif ext == '.lp':
        return nx_from_asp(asp_data if asp_data else graph_from_file(fname), **kwargs)
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")


def nx_from_asp(asp:str, directed:bool=True, with_reaction_nodes:bool=False, quoted_names:bool=False) -> nx.Graph or nx.Digraph:
    """Return a nx.Digraph describing the graph given in ASP format,
    using node/1, edge/2, metabolite/1, reaction/1, product/2 and reactant/2 atoms."""
    quoted = utils.quoted if quoted_names else lambda x: x
    graph = (nx.DiGraph if directed else nx.Graph)()
    models = utils.solve(inline=asp).by_arity.discard_quotes
    for model in models:
        for node, in model.get('node/1', ()):
            graph.add_node(quoted(node))
        for src, trg in model.get('edge/2', ()):
            graph.add_edge(quoted(src), quoted(trg))
        for node, in model.get('metabolite/1', ()):
            graph.add_node(quoted(node))
        if with_reaction_nodes:
            for node, in model.get('reaction/1', ()):
                graph.add_node(quoted(node))
            for prd, rct in model.get('product/2', ()):
                graph.add_edge(quoted(rct), quoted(prd))
            for met, rct in model.get('reactant/2', ()):
                graph.add_edge(quoted(met), quoted(rct))
        else:  # no reaction node: directly link reactant to products.
            reactions = defaultdict(lambda: (set(), set()))  # name: ({reactants}, {products})
            for prd, rct in model.get('product/2', ()):
                reactions[rct][1].add(quoted(prd))
            for met, rct in model.get('reactant/2', ()):
                reactions[rct][0].add(quoted(met))
            for reactants, products in reactions.values():
                for reactant in reactants:
                    for product in products:
                        graph.add_edge(reactant, product)
    return graph


def sccs_dag_from_nxdigraph(graph:nx.DiGraph, sccs:dict=None) -> dict:
    "Return the mapping {scc id: {successor scc id}} from input graph"
    if sccs is None:
        sccs = {min(nodes): frozenset(nodes) for nodes in nx.strongly_connected_components(graph)}
    scc_dag = defaultdict(set)  # SCC identifier: successor SCCs
    edges, no_successors, non_root = frozenset(graph.edges), set(), set()
    for scc_a, nodes_a in sccs.items():
        has_successor = False
        for scc_b, nodes_b in sccs.items():
            if nodes_a is nodes_b: continue
            if any((noda, nodb) in edges for noda in nodes_a for nodb in nodes_b):
                scc_dag[scc_a].add(scc_b)
                has_successor = True
                non_root.add(scc_b)
        if not has_successor:
            no_successors.add(scc_a)
    roots = {scc for scc in scc_dag.keys() if scc not in non_root}
    # The SCCs that are in `sccs` but not linked to any other SCC are not yet handled,
    #  because they are not yet added to the DAG. All of them has no successor,
    #  hence are in `no_successors`, but careful not to add the DAG terminals
    all_dag_leafs = set(itertools.chain.from_iterable(scc_dag.values()))
    roots |= {scc for scc in sccs.keys() if scc in no_successors and scc not in all_dag_leafs}
    # Ensure the treatment of roots without childs:
    #  if a root is only in scc_dag[None], it will not be treated, so we
    #  add them with empty set of child in the scc_dag.
    for root in roots:
        if root not in scc_dag:
            scc_dag[root] = set()
    # last verifications before returning the DAG
    assert not non_root & roots, non_root & roots
    scc_dag[None] = roots
    dag_values = set(itertools.chain.from_iterable(scc_dag.values())) | roots
    assert dag_values == set(sccs), (dag_values, set(sccs))
    return dict(scc_dag)


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
