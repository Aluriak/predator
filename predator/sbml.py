"""Routines to extract information from SBML files.

"""

import re
import xml.etree.ElementTree as etree
from xml.etree.ElementTree import XML, fromstring, tostring


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


def readSBMLnetwork(filename,
                    use_topological_injections:bool=True,
                    use_semantic_injection:bool=False):
    """Yield ASP atoms describing the network found in given SBML network

    filename -- sbml/xml file name.
    use_topological_injections -- consider as seed any species produced by reactant-free reactions.
    use_semantic_injection -- consider as seed any species marked as such in the SBML data.

    """
    tree = etree.parse(filename)
    sbml = tree.getroot()
    model = get_model(sbml)

    reactions = get_listOfReactions(model)
    for e in reactions:
        tag = get_sbml_tag(e)
        if tag == "reaction":
            reactionId = e.attrib.get("id")
            yield f'dreaction("{reactionId}").'
            reversible = e.attrib.get("reversible") == 'true'
            if reversible:
                yield f'reversible("{reactionId}").'

            reactants = get_listOfReactants(e)
            products = get_listOfProducts(e)

            for reactant in reactants:
                name = reactant.attrib.get('species')
                yield f'reactant("{name}","{reactionId}").'
                if use_topological_injections and reversible and not products:  # this reaction is a generator of seed
                    yield f'seed("{name}").'

            for product in products:
                name = product.attrib.get('species')
                yield f'product("{name}","{reactionId}").'
                if use_topological_injections and not reactants:  # this reaction is a generator of seed
                    yield f'seed("{name}").'


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
            lpfacts.add(Term("dreaction", ['"' + reactionId + '"']))  # , "\""+name+"\""
            if e.attrib.get("reversible") == "true":
                reversibool = True
                lpfacts.add(Term("dreaction", ['"' + reactionId + '_rev"']))
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
