"""Routines to extract information from SBML files.

"""

import itertools
import networkx as nx
import xml.etree.ElementTree as etree
from xml.etree.ElementTree import XML, fromstring, tostring
from . import utils


def get_sbml_tag(element) -> str:
    "Return tag associated with given SBML element"
    if element.tag[0] == "{":
        _, tag = element.tag[1:].split("}")  # uri is not used
    else:
        tag = element.tag
    return tag


def get_model(sbml):
    """
    return the model of a SBML
    """
    model_element = None
    for e in sbml:
        tag = get_sbml_tag(e)
        if tag == "model":
            model_element = e
            break
    return model_element


def get_listOfSpecies(model):
    """
    return list of species of a SBML model
    """
    listOfSpecies = None
    for e in model:
        tag = get_sbml_tag(e)
        if tag == "listOfSpecies":
            listOfSpecies = e
            break
    return listOfSpecies


def get_listOfReactions(model) -> list:
    """return list of reactions of a SBML model"""
    listOfReactions = []
    for e in model:
        tag = get_sbml_tag(e)
        if tag == "listOfReactions":
            listOfReactions = e
            break
    return listOfReactions


def get_listOfReactants(reaction) -> list:
    """return list of reactants of a reaction"""
    listOfReactants = []
    for e in reaction:
        tag = get_sbml_tag(e)
        if tag == "listOfReactants":
            listOfReactants = e
            break
    return listOfReactants


def get_listOfProducts(reaction) -> list:
    """return list of products of a reaction"""
    listOfProducts = []
    for e in reaction:
        tag = get_sbml_tag(e)
        if tag == "listOfProducts":
            listOfProducts = e
            break
    return listOfProducts


def read_SBML_network(filename,
                      use_topological_injections:bool=True,
                      use_semantic_injection:bool=False):
    """Yield ASP atoms describing the network found in given SBML network

    filename -- sbml/xml file name.
    graph -- if given, nx.Digraph instance to populate with the network
    use_topological_injections -- consider as seed any species produced by reactant-free reactions.
    use_semantic_injection -- consider as seed any species marked as such in the SBML data.

    """
    tree = etree.parse(filename)
    sbml = tree.getroot()
    model = get_model(sbml)

    def make_reaction(rxn_name, reactants, products, reversible:bool):
        """Add reverse reaction"""
        yield f'reaction("{rxn_name}").'

        for reactant in reactants:
            rname = reactant.attrib.get('species')
            yield f'reactant("{rname}","{rxn_name}").'
            if use_topological_injections and reversible and not products:  # this reaction is a generator of seed
                yield f'seed("{rname}").'

        for product in products:
            pname = product.attrib.get('species')
            yield f'product("{pname}","{rxn_name}").'
            if use_topological_injections and not reactants:  # this reaction is a generator of seed
                yield f'seed("{pname}").'

    reactions = get_listOfReactions(model)
    for e in reactions:
        tag = get_sbml_tag(e)
        if tag == "reaction":
            reactionId = e.attrib.get("id")
            reactants = get_listOfReactants(e)
            products = get_listOfProducts(e)
            reversible = e.attrib.get("reversible") == 'true'
            yield from make_reaction(reactionId, reactants, products, reversible)

            if reversible:
                yield from make_reaction('rev_' + reactionId, products, reactants, reversible)
                yield f'reversible("{reactionId}").'


def read_SBML_network_as_simple_graph(filename:str, with_reaction_nodes:bool=True, quoted_names:bool=True) -> nx.DiGraph:
    """Return nx.Digraph describing the network found in given SBML network

    filename -- sbml/xml file name.
    with_reaction -- if set, reactions will be nodes. Else, reactants are directly linked to products
    quoted_names -- if set, all nodes will be utils.quoted

    """
    tree = etree.parse(filename)
    sbml = tree.getroot()
    model = get_model(sbml)
    digraph = nx.DiGraph()
    quoted = utils.quoted if quoted_names else lambda x: x

    if with_reaction_nodes:
        def make_reaction(rxn_name, reactants, products):
            """Add reverse reaction"""
            digraph.add_node(quoted(rxn_name), isreaction=True)
            for reactant in reactants:
                rname = reactant.attrib.get('species')
                digraph.add_edge(quoted(rname), quoted(rxn_name))
            for product in products:
                pname = product.attrib.get('species')
                digraph.add_edge(quoted(rxn_name), quoted(pname))
    else:  # don't include reaction nodes in the graph
        def make_reaction(rxn_name, reactants, products):
            """Add reverse reaction"""
            reactants = tuple(reactant.attrib.get('species') for reactant in reactants)
            products = tuple(product.attrib.get('species') for product in products)
            for reactant, product in itertools.product(reactants, products):
                digraph.add_edge(quoted(reactant), quoted(product))

    reactions = get_listOfReactions(model)
    for e in reactions:
        tag = get_sbml_tag(e)
        if tag == "reaction":
            reactionId = e.attrib.get("id")
            reactants = get_listOfReactants(e)
            products = get_listOfProducts(e)
            reversible = e.attrib.get("reversible") == 'true'
            make_reaction(reactionId, reactants, products)

            if reversible:
                make_reaction('rev_' + reactionId, products, reactants)

    return digraph


def readSBMLnetwork_irrev(filename, name):
    """
    Read a SBML network and turn it into ASP-friendly data
    """
    lpfacts = TermSet()

    tree = etree.parse(filename)
    sbml = tree.getroot()
    model = get_model(sbml)

    reactions = get_listOfReactions(model)
    for e in reactions:
        reversibool = False
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "reaction":
            reactionId = e.attrib.get("id")
            lpfacts.add(Term("reaction", ['"' + reactionId + '"']))  # , "\""+name+"\""
            if e.attrib.get("reversible") == "true":
                reversibool = True
                lpfacts.add(Term("reaction", ['"' + reactionId + '_rev"']))
                # lpfacts.add(Term('reversible', ["\""+reactionId+"\""]))

            listOfReactants = get_listOfReactants(e)
            if listOfReactants == None:
                print("\n Warning:", reactionId, "listOfReactants=None")
            else:
                for r in listOfReactants:
                    lpfacts.add(Term('reactant', ["\""+r.attrib.get("species")+"\"", "\""+reactionId+"\""])) #,"\""+name+"\""
                    if reversibool:
                        lpfacts.add(Term('product', ["\""+r.attrib.get("species")+"\"", "\""+reactionId+"_rev\""])) #,"\""+name+"\""

            listOfProducts = get_listOfProducts(e)
            if listOfProducts == None:
                print("\n Warning:", reactionId, "listOfProducts=None")
            else:
                for p in listOfProducts:
                    lpfacts.add(Term('product', ["\""+p.attrib.get("species")+"\"", "\""+reactionId+"\""])) #,"\""+name+"\""
                    if reversibool:
                        lpfacts.add(Term('reactant', ["\""+p.attrib.get("species")+"\"", "\""+reactionId+"_rev\""]))
    #print(lpfacts)
    return lpfacts
