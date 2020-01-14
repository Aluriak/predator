"""Routines to search for seeds.

The main API is search_seeds() function, that dispatch call
to the suited routine depending of overall parameters.

A routine should have the following signature:

def search_seed_method_name(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration) -> [{str}]

It can have other parameters, that search_seeds() must take care of.
The yielded values are sets of strings, each set describing a solution,
i.e. a set of seed names that fullfill the asked conditions.

"""
import math
import clyngor
import networkx as nx
import itertools
from clyngor import opt_models_from_clyngor_answers
from . import graph as graph_module
from .utils import quoted, unquoted, solve, get_terminal_nodes, inverted_dag, remove_terminal, render_network
from functools import partial
from collections import defaultdict
from pkg_resources import resource_filename
from enum import Enum


class EnumMode(Enum):
    "The mode of solution enumeration"
    Enumeration, Union, Intersection = 'enumeration', 'union', 'intersection'

    @property
    def clingo_option(self) -> str:
        return {
            'enumeration': '',
            'union': '--enum-mode brave',
            'intersection': '--enum-mode cautious',
        }[self.value]

    @property
    def clingo_option_for_iterative_search(self) -> str:
        return {
            'enumeration': '',
            'union': '--enum-mode brave',
            'intersection': '',
        }[self.value]

    @property
    def frozenset_operation(self) -> str:
        return {
            'enumeration': lambda sets: sets,
            'union': lambda sets: {frozenset.union(*sets)},
            'intersection': lambda sets: {frozenset.intersection(*sets)},
        }[self.value]


ASP_SRC_ENUM_CC = resource_filename(__name__, 'asp/enum-cc.lp')
ASP_SRC_SIMPLE_SEED_SOLVING = resource_filename(__name__, 'asp/simple-seed-solving.lp')
ASP_SRC_GREEDY_TARGET_SEED_SOLVING = resource_filename(__name__, 'asp/greedy-target-seed-solving.lp')
ASP_SRC_GREEDY_TARGET_SEED_SOLVING_SEED_MINIMALITY = resource_filename(__name__, 'asp/greedy-target-seed-solving-seed-minimality-constraint.lp')
ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM = resource_filename(__name__, 'asp/iterative-target-seed-solving--aim.lp')
ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM__MINIMALITY_CONSTRAINT = resource_filename(__name__, 'asp/iterative-target-seed-solving--aim--minimality-constraint.lp')
ASP_SRC_PARETO_OPTIMIZATIONS = resource_filename(__name__, 'asp/pareto_utnu_scope_seeds.lp')
ASP_SRC_PARETO_OPTIMIZATIONS_ALL = resource_filename(__name__, 'asp/pareto_utnu_avoid_targets_as_seeds.lp')
ASP_SRC_ITERATIVE_PARETO_OPTIMIZATIONS = resource_filename(__name__, 'asp/pareto_iterative_scope_seeds.lp')
ASP_SRC_ITERATIVE_PARETO_OPTIMIZATIONS_ALL = resource_filename(__name__, 'asp/pareto_iterative_avoid_targets_as_seeds.lp')

import os
assert os.path.exists(ASP_SRC_ENUM_CC), ASP_SRC_ENUM_CC
assert os.path.exists(ASP_SRC_SIMPLE_SEED_SOLVING), ASP_SRC_SIMPLE_SEED_SOLVING


def search_seeds(graph_data:str='', start_seeds:iter=(), forbidden_seeds:iter=(),
                 targets:set=(), graph_filename:str=None, enum_mode:str=EnumMode.Enumeration,
                 explore_pareto:bool=False, pareto_no_target_as_seeds:bool=False,
                 greedy:bool=False, **kwargs) -> [{set}]:
    "Wrapper around all seeds search methods. The used methods depends of given parameters."
    if not graph_data and graph_filename:
        graph_data = graph_module.graph_from_file(graph_filename)
    if not graph_data and not graph_filename:
        raise ValueError("No input data provided: expecting graph_data or graph_filename argument")
    enum_mode = EnumMode(enum_mode)
    if not targets:  # no target, just activate everything
        func = search_seeds_activate_all
        if 'compute_optimal_solutions' in kwargs:  del kwargs['compute_optimal_solutions']
    elif explore_pareto or pareto_no_target_as_seeds:  # explore the pareto front
        func = search_pareto_front if greedy else search_iterative_pareto_front
        kwargs.setdefault('avoid_targets_as_seeds', pareto_no_target_as_seeds)
    elif greedy:  # non efficient search of targets
        func = search_seeds_activate_targets_greedy
        if 'compute_optimal_solutions' in kwargs:  del kwargs['compute_optimal_solutions']
    else:  # efficient search of targets
        func = search_seeds_activate_targets_iterative
    if kwargs.get('verbose'):
        print('FUNC:', func)  # search_seeds_activate_targets_iterative
        print('KWARGS:', kwargs)
    yield from func(graph_data, frozenset(start_seeds), frozenset(forbidden_seeds), frozenset(targets), graph_filename, enum_mode, **kwargs)


def search_iterative_pareto_front(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                                  graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration,
                                  avoid_targets_as_seeds:bool=False, verbose:bool=False, **kwargs) -> [set]:
    """Yield the set of seeds/optimal solutions.

    This is a wrapper around search_seeds_activate_targets_iterative() function,
    where a pareto exploration is used instead of the standard
    minimize-seed-number approach for each SCC.

    This function basically prepare the `compute_hypothesis_from_scc` argument
    and yield from the decorated function call.

    """
    # the scc->hypothesis function relies on the pareto approach if asked
    func_scc2hyp = partial(_compute_hypothesis_from_scc__pareto, avoid_targets_as_seeds=avoid_targets_as_seeds)
    if 'avoid_targets_as_seeds' in kwargs:  # already passed to func_scc2hyp
        del kwargs['avoid_targets_as_seeds']
    yield from search_seeds_activate_targets_iterative(
        graph_data, start_seeds, forbidden_seeds, targets, graph_filename,
        enum_mode, compute_hypothesis_from_scc=func_scc2hyp, verbose=verbose,
        **kwargs
    )

def search_pareto_front(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                        graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration,
                        avoid_targets_as_seeds:bool=False, verbose:bool=False) -> [set]:
    """Yield the set of seeds for each found solution on the pareto front."""
    if enum_mode in {EnumMode.Union, EnumMode.Intersection}:
        raise ValueError(f"Mode {str(enum_mode)} is not supported for this routine.")
    # treat parameters
    if start_seeds and forbidden_seeds and start_seeds & forbidden_seeds:
        raise ValueError(f"start_seeds and forbidden_seeds shares some seeds: {start_seeds & forbidden_seeds}.")
    targets = frozenset(map(quoted, targets))
    start_seeds = frozenset(map(quoted, start_seeds))
    forbidden_seeds = frozenset(map(quoted, forbidden_seeds))
    # data representation
    start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
    targets_repr = ' '.join(f'target({quoted(t)}).' for t in targets)
    forb_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
    data_repr = graph_data + start_seeds_repr + targets_repr + forb_repr
    # solving
    models = _pareto_with_asprin(data_repr, iterative=False, verbose=verbose)
    for model, opts, isoptimum in models:
        if isoptimum:
            yield frozenset(args[0] for args in model.get('seed', ()))

def _pareto_with_asprin(data:str, avoid_targets_as_seeds:bool=False, discard_quotes:bool=True, iterative:bool=False, verbose:bool=False) -> clyngor.Answers:
    """Return ASP answers found by calling asprin on given data for (iterative)
    seed search.
    """
    if iterative:
        pareto_constraints = ASP_SRC_ITERATIVE_PARETO_OPTIMIZATIONS_ALL if avoid_targets_as_seeds else ASP_SRC_ITERATIVE_PARETO_OPTIMIZATIONS
        files = pareto_constraints, ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM
    else:
        pareto_constraints = ASP_SRC_PARETO_OPTIMIZATIONS_ALL if avoid_targets_as_seeds else ASP_SRC_PARETO_OPTIMIZATIONS
        files = pareto_constraints, ASP_SRC_GREEDY_TARGET_SEED_SOLVING
    models = solve(files, inline=data, clingo_bin_path='asprin').with_optimality.by_predicate
    if discard_quotes: models.discard_quotes
    if verbose:
        print('_pareto_with_asprin:')
        print('\tDATA:', data)
        print('\tCMD:', models.command)
    return models

def search_seeds_activate_targets_iterative(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                                            graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration,
                                            compute_optimal_solutions:bool=False, filter_included_solutions:bool=True,
                                            compute_hypothesis_from_scc:callable=None, verbose:bool=False,
                                            sccs:dict=None, scc_dag:dict=None) -> [set]:
    """Yield the set of seeds for each found solution.

    compute_optimal_solutions -- if True, will post-process the solutions to get
                                 the optimal solutions.
    filter_included_solutions -- post process solution list to discard solutions
                                 that are subset of others.

    This implements the activation of targets: find the minimum sets
    of seeds that activate all targets.
    This is an iterative implementation, starting from SCC DAG leaves,
    finishing on roots.

    Definition: *aim*
        in a SCC, an *aim* is the union of the sets of (1) targets inside the SCC,
        and (2) outgoing reactions that must be activated.
        It is therefore dependent of the context (named later Hypothesis)
        in which the SCC is solved.

    Proposition: an SCC with no *aim* is either a terminal SCC without target,
        or any SCC with subsequent SCC carrying no targets or that need no input to activate it
        (because seeds are producing the subsequent targets).
        Such SCC can be ignored during the computation.

    Algorithm:
        While there is non treated SCC:
            for each terminal SCC named T:
                for each hypothesis H involving T:
                    compute the *aim* of T in given H
                    decide S, the set of hypothesis that fullfill the *aim*
                    for each hypothesis of S:
                        add hypothesis to the set of all hypothesis (parents will use it)
                remove T from the set of SCC, recompute set of terminal SCC accordingly.

    An *Hypothesis* is an object created by childs and propagated toward parents.
    The Hypothesis attached to root are the final description of solutions.

    An Hypothesis is:

        set of seed, {scc -> {reaction name -> reactants in scc})}, set of fullfilled targets

    The *set of seed* contains *all* seeds accross the SCCs that need to be activated according to the hypothesis,
    and constitute a whole solution.

    The *set of fullfilled targets* contains the targets that have been activated in the hypothesis.

    The *scc* found in an Hypothesis are always the parent of the last treated
    SCCs: an SCC cannot propagate an Hypothesis with self in the Hypothesis.

    The set of all Hypothesis H is maintained during the propagation, so that
    a SCC will extract (and remove) from H the hypothesis that are implying her.
    For each hypothesis extracted, the SCC will add to H one or more hypothesis.

    Other properties:

    1. When solving an scc, all listed reactions must be activated. The reactants in the scc are therefore targets.
    2. The reactants outside the scc will necessarily be activated: therefore the reactants outside the scc are considered activated.
    3. The current SCC will be removed from hypothesis during propagation, replaced by a seed or another reaction in a parent SCC.
    4. However, another SCC cannot be modified: an SCC removes itself from the Hypothesis, and add its parents.
    5. All Hypothesis objects are stored globally, and an SCC needs first to extract those that concern itself.

    Optimizations to think about:

    1. do not add an hypothesis to H when a better one is already in H
    2. EnumMode.Intersection is currently naïve. Maybe it could be treated specifically like Union.
    3. when an SCC is composed of one node, with one ingoing and one outgoing reaction, ASP is unneeded.

    """
    compute_hypothesis_from_scc = compute_hypothesis_from_scc or _compute_hypothesis_from_scc__minimization
    _print = print if verbose else lambda *_, **__: None
    _print(start_seeds, forbidden_seeds, start_seeds & forbidden_seeds)
    if start_seeds and forbidden_seeds and start_seeds & forbidden_seeds:
        raise ValueError(f"start_seeds and forbidden_seeds shares some seeds: {start_seeds & forbidden_seeds}.")

    def find_aim(scc_name:str, hypothesis:tuple) -> [str]:
        "Yield targets or reactions associated with given scc for given hypothesis"
        scc_targets = targets & sccs[scc_name]
        seeds, scc_reactions, fullfilled = hypothesis
        for scc, reactions in scc_reactions.items():
            if scc == scc_name:
                for reaction, external_reactants in reactions.items():
                    scc_targets |= external_reactants
        _print('\tFIND AIM:', scc_name, hypothesis, '->', scc_targets)
        return scc_targets

    all_hypothesis = []  # list of hypothesis: (seeds, {scc: {reaction name: reactants in scc}}, fullfilled targets)

    def get_hypothesis_of(scc_name:str) -> iter:
        """Yield hypothesis that are associated with given scc_name,
        removing them from all_hypothesis list"""
        nonlocal all_hypothesis  # TODO: maybe there is a better way to maintain all_hypothesis
        new_all_hypothesis = []  # will replace all_hypothesis
        for hypothesis in all_hypothesis:
            if scc_name in hypothesis[1]:
                yield hypothesis
            else:  # keep the hypothesis
                new_all_hypothesis.append(hypothesis)
        all_hypothesis = new_all_hypothesis

    def get_null_hypothesis(scc_name:str=None) -> tuple:
        "return a null hypothesis object"
        hyp = set(), {}, set()
        if scc_name:  hyp[1].setdefault(scc_name, dict())
        return hyp

    def get_parents_of(scc_name:str) -> iter:
        """Yield parent SCCs. If no parent, yield None"""
        yield from rev_scc_dag[scc_name] or [None]

    if not sccs or not scc_dag:
        _print('COMPUTE SCCs…', end='', flush=True)
        sccs, scc_dag = compute_sccs(graph_data, graph_filename=graph_filename, verbose=verbose)
        _print(' OK!')
    terminals = frozenset(get_terminal_nodes(scc_dag))
    rev_scc_dag = dict(inverted_dag(scc_dag))
    _print('  SCC DAG:', scc_dag)
    # _print('TERMINALS:', terminals)
    # _print('         :', frozenset(rev_scc_dag[None]))
    assert terminals == frozenset(rev_scc_dag[None])
    # hypothesis = defaultdict(list)  # {scc: hypothesis} with hypothesis == iterable of (seeds, reactions)
    # associate an empty hypothesis for each SCC having a target.
    targets = frozenset(map(quoted, targets))
    for scc, nodes in sccs.items():
        _print('SEARCHING TARGETS:', nodes, targets, nodes & targets)
        if nodes & targets:  # the terminal SCC has an aim
            all_hypothesis.append(get_null_hypothesis(scc))
    # iteratively find hypothesis
    while len(scc_dag) > 1:  # last valid key is None
        for terminal in frozenset(get_terminal_nodes(scc_dag)):
            _print('\nTERMINAL:', terminal)
            if terminal is None: continue
            self_hypothesis = tuple(get_hypothesis_of(terminal))
            _print('ALL HYPS:', 'len:', len(all_hypothesis), all_hypothesis, '\t(does not contains TERM HYP)')
            _print('TERM HYP:', 'len:', len(self_hypothesis), self_hypothesis)
            _print('    DAG :', 'len:', len(scc_dag), scc_dag)
            for current_hypothesis in self_hypothesis:
                aim = find_aim(terminal, current_hypothesis)
                if not aim:  # the hypothesis has nothing to produce
                    seeds, scc_reactions, fullfilled = current_hypothesis
                    for parent in get_parents_of(terminal):
                        new_scc_reactions = dict(scc_reactions)  # don't share among parents
                        new_scc_reactions[parent] = scc_reactions[terminal]
                        del new_scc_reactions[terminal]
                        new_hypothesis = set(seeds), new_scc_reactions, fullfilled & targets
                        all_hypothesis.append(new_hypothesis)
                        _print(f"\tTerminal {terminal} has nothing to do ; propagation of hypothesis '{new_hypothesis}' to parent {parent}.")
                        break  # don't propagate such trivial hypothesis to all parents. Only one will have the job done
                    continue  # the current hypothesis has been treated
                start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
                targets_repr = ' '.join(f'target({quoted(t)}).' for t in aim)
                forb_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
                scc_repr = ' '.join(f'scc({terminal},{node}).' for node in sccs[terminal])
                scc_data = graph_data + f'current_scc({terminal}). {scc_repr} {start_seeds_repr} {forb_repr} {targets_repr}'
                _print('\tCURRENT HYP:', current_hypothesis)
                _print('\tCURR.  AIM :', aim)
                # _print('\t  SCC DATA:', scc_data)
                for new_seeds, new_targets, new_fullfill in compute_hypothesis_from_scc(terminal, scc_data, sccs, rev_scc_dag, enum_mode, verbose=verbose):
                    new_seeds |= current_hypothesis[0]
                    new_targets = {**current_hypothesis[1], **new_targets}
                    del new_targets[terminal]  # remove self from the hypothesis
                    new_fullfill |= current_hypothesis[2]
                    new_hypothesis = new_seeds, new_targets, new_fullfill & targets
                    _print('\t\tADDED HYPOTHESIS:', new_hypothesis)
                    all_hypothesis.append(new_hypothesis)
            # now, remove terminal from the dag
            remove_terminal(terminal, scc_dag, frozenset(rev_scc_dag.get(terminal, ())))
    _print('OUTPUT HYPOTHESIS:', all_hypothesis, targets, 'compute_optimal_solutions=', compute_optimal_solutions)
    solutions = _solutions_from_hypothesis(all_hypothesis, targets, compute_optimal_solutions, filter_included_solutions, enum_mode, verbose)
    return enum_mode.frozenset_operation(solutions)

def _solutions_from_hypothesis(all_hypothesis:list, targets:set, compute_optimal_solutions:bool, filter_included_solutions:bool, enum_mode:EnumMode, verbose:bool) -> frozenset:
    """Compute the solutions from hypothesis and targets.
    If compute_optimal_solutions is given, will filter out all non-optimal solutions.

    """
    _print = print if verbose else lambda *_, **__: None
    if not all_hypothesis:  # no hypothesis
        return frozenset()
    # verify that all targets are reachable
    reachables = set()
    for _, __, reachable in all_hypothesis:
        reachables |= targets & reachable
    if reachables != targets:  # not all are reachables
        # _print('WARNING: Not all targets are reachables:', ', '.join(targets - reachables) + '. An error will be raised.')
        raise ValueError(f"Not all targets are reachables: {targets - reachables}.")
        # targets = reachables  # TODO: should raise an error, instead ?
    # quick data integrity check
    unfullfilled_targets = set(targets)
    for seeds, scc_reactions, fullfilled in all_hypothesis:
        assert len(scc_reactions) == 1, scc_reactions
        assert None in scc_reactions, scc_reactions
        assert len(scc_reactions[None]) == 0, scc_reactions[None]
        if enum_mode is not EnumMode.Intersection:
            assert seeds, seeds
        assert fullfilled, (seeds, fullfilled)
        seeds = frozenset(map(unquoted, seeds))
        unfullfilled_targets -= fullfilled
        # _print(f'SOLUTION: seeds {set(seeds)} are fullfilling {set(map(unquoted, fullfilled))}')
    assert not unfullfilled_targets, unfullfilled_targets  # some targets where not reachable ; what is this magic ?
    # collect, for each target, the set of seeds reaching it
    seeds_for_target = defaultdict(list)  # target: {seeds}
    for seeds, scc_reactions, fullfilled in all_hypothesis:
        for target in fullfilled:
            seeds_for_target[target].append(seeds)
    # generate all possible solutions
    # _print('SEEDS FOR TARGETS:', seeds_for_target)
    # _print('                 :', tuple(itertools.product(*seeds_for_target.items())))
    solutions = set()  # a solution is a set of seeds activating all targets
    for seeds_sets in itertools.product(*seeds_for_target.values()):
        # _print('SSEEDS SETS:', seeds_sets)
        seeds_set = set.union(*seeds_sets)
        solutions.add(frozenset(map(lambda x: x.strip('"'), seeds_set)))
        _print(f'OPT SOLUTION: seeds {seeds_set} are fullfilling {targets}')
    if filter_included_solutions:
        filtered_solutions = set()  # a solution is a set of seeds activating all targets
        for seeds_sets in solutions:
            if not any(other_set < seeds_sets for other_set in solutions):
                filtered_solutions.add(seeds_sets)
        solutions = filtered_solutions
        # _print('FILTERED:', solutions)
    if compute_optimal_solutions:
        # reduce solutions using seed number
        opt_seed_number = min(map(len, solutions))
        solutions = (s for s in solutions if len(s) == opt_seed_number)
    return frozenset(solutions)


def __compute_hypothesis_from_scc(models:iter, scc_name:str, scc_encoding:set, sccs:dict, rev_scc_dag:dict, enum_mode:EnumMode, verbose:bool) -> [(set, dict, set)]:
    """Yield hypothesis derived from given models, computed from either clingo or asprin.

    For a pareto front exploration, see _compute_hypothesis_from_scc__pareto.
    For a seed minimization, see _compute_hypothesis_from_scc__minimization.

    """
    # the following call will provide us a model for each hypothesis.
    _print = print if verbose else lambda *_, **__: None
    _print('\tCOMPUTE ALL HYPOTHESIS:', scc_name)
    _print('                NB MODEL:', len(models))
    _print('             SCC PARENTS:', rev_scc_dag[scc_name])
    for new_hypothesis in models:
        # new targets will need to be activated by parent SCC.
        new_targets = defaultdict(set)  # reaction: reactants
        new_seeds, new_fullfilled = set(), set()
        for args in new_hypothesis.get('new_target', ()):
            if len(args) == 2:
                target, reaction = args
                new_targets[reaction].add(target)
        for args in new_hypothesis.get('seed', ()):
            if len(args) == 1:
                new_seeds.add(args[0])
        for args in new_hypothesis.get('activated_local_target', ()):
            if len(args) == 1:
                new_fullfilled.add(args[0])
        _print('\t\tFound Hypothesis:   targets:', new_targets, '\tseeds:', new_seeds, '\tfullfilled:', new_fullfilled)
        # create hypothesis with each parent SCC that have a reactant in it (or None for roots)
        scc_reactions = {None: dict(new_targets)}  # default case: no parent
        alien_reactants = frozenset(itertools.chain.from_iterable(new_targets.values()))
        _print('\t\tSCC_REACTION BUILDING…')
        _print('\t\t\t', rev_scc_dag[scc_name], new_targets)
        _print('\t\t\t', alien_reactants)
        if rev_scc_dag[scc_name] and new_targets:
            scc_reactions = {parent: dict(new_targets) for parent in rev_scc_dag[scc_name]
                             if any(reactant in sccs[parent] for reactant in alien_reactants)}
        _print('\t\t\t', scc_reactions)
        yield new_seeds, scc_reactions, new_fullfilled

def _compute_hypothesis_from_scc__minimization(scc_name:str, scc_encoding:set, sccs:dict, rev_scc_dag:dict, enum_mode:EnumMode, verbose:bool) -> [(set, dict, set)]:
    """Yield hypothesis computed from given scc_name to consider for next SCCs,
    using a seed minimization approach.

    For a more subtle implementation using pareto fronts,
    see _compute_hypothesis_from_scc__pareto implementation.

    """
    # the following call will provide us a model for each hypothesis.
    _print = print if verbose else lambda *_, **__: None
    files = ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM, ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM__MINIMALITY_CONSTRAINT
    models = solve(files, inline=scc_encoding, options='--opt-mode=optN ' + enum_mode.clingo_option_for_iterative_search, delete_tempfile=False).by_predicate
    _print('CMD:', models.command)
    if enum_mode is EnumMode.Union:  # optimized case
        model = None
        for model in models:  pass
        models = [model] if model else []
    else:  # for intersection and enumeration, just get optimums
        models = tuple(opt_models_from_clyngor_answers(models))
    yield from __compute_hypothesis_from_scc(models, scc_name, scc_encoding, sccs, rev_scc_dag, enum_mode, verbose)


def _compute_hypothesis_from_scc__pareto(scc_name:str, scc_encoding:set, sccs:dict, rev_scc_dag:dict, enum_mode:EnumMode, verbose:bool, avoid_targets_as_seeds:bool=False) -> [(set, dict, set)]:
    """Yield hypothesis computed from given scc_name to consider for next SCCs,
    using a pareto front exploration to discover the different hypothesis.

    Make use of _pareto_with_asprin() and __compute_hypothesis_from_scc() functions.

    """
    # the following call will provide us a model for each hypothesis.
    _print = print if verbose else lambda *_, **__: None
    models = _pareto_with_asprin(scc_encoding, avoid_targets_as_seeds=avoid_targets_as_seeds, discard_quotes=False, iterative=True, verbose=verbose)
    if enum_mode is EnumMode.Union:  # optimized case
        raise NotImplementedError("Iterative Pareto search with Union as enum_mode")
        models = tuple(models)
        assert len(models) == 1
    else:  # for intersection and enumeration, just get optimums
        pass  # optimums are already provided by pareto search
        models = (model for model, opts, isoptimum in models if isoptimum)
    models = tuple(models)
    yield from __compute_hypothesis_from_scc(models, scc_name, scc_encoding, sccs, rev_scc_dag, enum_mode, verbose)


def search_seeds_activate_targets_greedy(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                                         graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration,
                                         compute_optimal_solutions:bool=False, filter_included_solutions:bool=True) -> [{set}]:
    """Yield the set of seeds for each found solution.

    This implements the activation of targets: find the minimum sets
    of seeds that activate all targets.
    This is a greedy implementation. Do not expect it to work on a large dataset.

    Both compute_optimal_solutions and filter_included_solutions are unused.
    They are here only to mimic search_seeds_activate_targets_iterative API.

    """
    if not targets:
        raise ValueError("search_seeds_activate_targets_greedy() requires targets. Use another function, search_seeds(), or provide targets.")
    start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
    targets_repr = ' '.join(f'target({quoted(t)}).' for t in targets)
    forb_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
    data_repr = graph_data + start_seeds_repr + forb_repr + targets_repr
    models = solve((ASP_SRC_GREEDY_TARGET_SEED_SOLVING, ASP_SRC_GREEDY_TARGET_SEED_SOLVING_SEED_MINIMALITY),
                   inline=data_repr, options='--opt-mode=optN ' + enum_mode.clingo_option).discard_quotes.by_predicate
    if enum_mode is EnumMode.Enumeration:
        models = opt_models_from_clyngor_answers(models)
    else:
        for model in models:
            _models = [model]
        models = _models
    for model in models:
        seeds = frozenset(args[0] for args in model['seed'] if len(args) == 1)
        yield seeds


def search_seeds_activate_all(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration) -> [{set}]:
    """Yield the set of seeds for each found solution.

    This implements the most simple solution: no targets, find the minimum sets
    of seeds that activate everything.

    Use SCC optimization, starting from SCC root, going to leaves.

    """
    if targets:
        raise ValueError("search_seeds_activate_all() does not handle targets. Use another method that supports it, or search_seeds().")
    start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
    forbidden_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
    assert isinstance(graph_data, str), graph_data
    sccs, scc_dag = compute_sccs(graph_data, graph_filename=graph_filename)
    roots = scc_dag[None]
    scc_seeds = {}  # scc name: list of {optimal seeds}
    for scc_name, nodes in sccs.items():
        scc_repr = ' '.join(f'scc({scc_name},{node}).' for node in nodes)
        scc_data = f'current_scc({scc_name}). {scc_repr} {start_seeds_repr} {forbidden_repr}'
        # print('DATA:\n    ' + scc_data)
        # print('    ' + graph_data)
        # print()
        models = solve(ASP_SRC_SIMPLE_SEED_SOLVING, inline=graph_data + scc_data, options='--opt-mode=optN ' + enum_mode.clingo_option, delete_tempfile=False)
        models = opt_models_from_clyngor_answers(models.by_predicate.discard_quotes)
        scc_seeds[scc_name] = tuple(frozenset(args[0] for args in model.get('seed', ())) for model in models)
        # print('OUTPUT SEEDS:', scc_seeds[scc_name])
    # generate all possibilities
    sets_of_seeds_sets = scc_seeds.values()
    # seeds_set is a combination of possible set of seeds for each SCC
    seed_combinations = (frozenset.union(*seeds_sets) for seeds_sets in itertools.product(*sets_of_seeds_sets))
    if enum_mode is EnumMode.Union:
        seed_combinations = [frozenset.union(*seed_combinations)]
    elif enum_mode is EnumMode.Intersection:
        seed_combinations = [frozenset.intersection(*seed_combinations)]
    yield from seed_combinations


def compute_sccs(graph_data:str, graph_filename:str=None, verbose:bool=False) -> [{str}]:
    """Return the DAG and nodes of Strongly Connected Components found in given graph.

    graph_data -- ASP encoding of the metabolic network.
    graph_filename -- filename of SBML encoded graph, used to speed up SCC mining.
    return -- (SCC, DAG), with DAG a SCC id -> {SCC id sucessor} mapping,
              and SCC a SCC id -> {node} mapping.

    If no graph_filename given, ASP will be used to infer them from graph_data.

    """
    _print = print if verbose else lambda *_, **__: None
    _print('COMPUTE SCCS:', graph_filename)
    # compute the SCCs components and the SCC DAG
    if graph_filename:  # use nx.Digraph extraction to quickly extract connected components
        graph = graph_module.nxgraph_from_file(graph_filename, graph_data, with_reaction_nodes=False, quoted_names=True)
        yield_sccs = nx.strongly_connected_components or nx.kosaraju_strongly_connected_components
        sccs = {min(nodes): frozenset(nodes) for nodes in yield_sccs(graph)}
        scc_dag = graph_module.sccs_dag_from_nxdigraph(graph, sccs)
    else:  # no graph_filename, use pure ASP computation
        models = solve(ASP_SRC_ENUM_CC, inline=graph_data)
        for model in models.by_predicate:
            _print('SCC MODEL:', model)
            roots = {args[0] for args in model.get('noinput', ()) if len(args) == 1}
            sccs = defaultdict(set)  # SCC identifier: nodes in SCC
            for scc_name, node in model.get('scc', ()):
                sccs[scc_name].add(node)
            sccs = dict(sccs)
            # Compute the DAG of SCCs
            scc_dag = defaultdict(set)  # SCC identifier: successor SCCs
            for scc_name, scc_succ in model.get('sccedge', ()):
                scc_dag[scc_name].add(scc_succ)
            scc_dag[None] = roots
            sccs = dict(sccs)
            scc_dag = dict(scc_dag)
            break
        else:
            raise RuntimeError(r"ASP SCC enumeration yielded more than one model")
    _print('SCC FINAL:', sccs)
    _print()
    _print('         :', scc_dag)
    _print()
    return sccs, scc_dag


def render_scc_dag(out_filename:str, sccs, scc_dag, targets):
    scc_dag_graph = '\n'.join(
        f'link({"root" if parent is None else parent},{child}).'
        for parent, childs in scc_dag.items()
        for child in childs
    ) + '\n' + '\n'.join(
        f'annot(upper,{scc},"S").'
        for scc, nodes in sccs.items()
        if len(nodes) == 1
    ) + '\n' + '\n'.join(
        f'annot(lower,{scc},"T").'
        for scc, nodes in sccs.items()
        if any(quoted(target) in nodes for target in targets)
    ) + (
        '\ndot_property(root,color,transparent).'
        '\nobj_property(edge,arrowhead,normal).'
    )
    return render_network(scc_dag_graph, out_filename)
