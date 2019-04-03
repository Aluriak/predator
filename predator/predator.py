"""Routines to search for seeds.

"""
import networkx as nx
from . import sbml as sbml_module


def search_seeds(graph_data:str, targets:set=None, graph_filename:str=None) -> [str]:
    """Yield the set of seeds for each found solution.

    """
    # compute SCC DAG
    # find source SCCs
    # for each SCC in width-first traversal:
    #     find the minimal set of seed in the SCC to fully 
    yield from ()  # TODO: implement that



def compute_sccs(graph_data:str, graph_filename:str=None) -> [{str}]:
    """Yield each Strongly Connected Component found in given graph.

    If no graph_filename given, ASP will be used to infer them from graph_data.

    """
    if not graph_filename:
        raise NotImplementedError("ASP search for connected components not implemented")

    else:  # use graph_filename
        digraph = sbml_module.read_SBML_network_as_simple_graph(graph_filename)
        for scc in nx.strongly_connected_components(digraph):
            yield frozenset(scc)
