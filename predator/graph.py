"""Routines for graph data extractions.

"""


import os
from . import sbml


def graph_from_file(fname:str, **kwargs) -> str:
    """Return ASP data encoding graph in given file."""
    ext = os.path.splitext(fname)[1]
    if ext in {'.xml', '.sbml'}:
        return sbml.readSBMLnetwork(fname, **kwargs)
    elif ext == '.lp':
        with open(fname) as fd:
            return fd.read()
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")
