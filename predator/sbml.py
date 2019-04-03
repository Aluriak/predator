"""Routines to extract information from SBML files.

"""

import re
import xml.etree.ElementTree as etree
from xml.etree.ElementTree import XML, fromstring, tostring


def get_model(sbml):
    """
    return the model of a SBML
    """
    model_element = None
    for e in sbml:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
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
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "listOfSpecies":
            listOfSpecies = e
            break
    return listOfSpecies


def get_listOfReactions(model):
    """
    return list of reactions of a SBML model
    """
    listOfReactions = None
    for e in model:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "listOfReactions":
            listOfReactions = e
            break
    return listOfReactions


def get_listOfReactants(reaction):
    """
    return list of reactants of a reaction
    """
    listOfReactants = None
    for e in reaction:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "listOfReactants":
            listOfReactants = e
            break
    return listOfReactants


def get_listOfProducts(reaction):
    """
    return list of products of a reaction
    """
    listOfProducts = None
    for e in reaction:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "listOfProducts":
            listOfProducts = e
            break
    return listOfProducts


def readSBMLnetwork(filename, name):
    """Yield ASP atoms describing the network found in given SBML network"""
    tree = etree.parse(filename)
    sbml = tree.getroot()
    model = get_model(sbml)

    reactions = get_listOfReactions(model)
    for e in reactions:
        if e.tag[0] == "{":
            uri, tag = e.tag[1:].split("}")
        else:
            tag = e.tag
        if tag == "reaction":
            reactionId = e.attrib.get("id")
            yield f'dreaction("{reactionId}").'
            if e.attrib.get("reversible") == "true":
                yield f'reversible("{reactionId}").'

            reactants = get_listOfReactants(e)
            if reactants == None:
                print("\n Warning:", reactionId, "reactants=None")
            else:
                for r in reactants:
                    yield f"reactant(\"{r.attrib.get('species')}\",\"{reactionId}\")."

            products = get_listOfProducts(e)
            if products == None:
                print("\n Warning:", reactionId, "products=None")
            else:
                for p in products:
                    yield f"product(\{p.attrib.get('species')}\",\"{reactionId}\")."


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
