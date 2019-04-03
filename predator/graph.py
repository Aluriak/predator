"""Routines for graph data extractions.

"""


import os
from . import sbml


def graph_from_file(fname:str) -> str:
    """Return ASP data encoding graph in given file."""
    basename, ext = os.path.splitext(fname)[1]
    if ext == '.sbml':
        return sbml.readSBMLnetwork(fname, os.path.split(basename)[1])
    elif ext == '.lp':
        with open(fname) as fd:
            return fd.read()
    else:
        raise NotImplementedError(f"File of extension {ext} cannot be treated")
