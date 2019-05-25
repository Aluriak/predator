"""Routines to search for seeds.

The main API is search_seeds() function, that dispatch call
to the suited routine depending of overall parameters.

A routine should have the following signature:

def search_seed_method_name(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration) -> [{str}]

It can have other parameters, that search_seeds() must take care of.
The yielded values are sets of strings, each set describing a solution, i.e. a set of seed names that fullfill the asked conditions.

"""
import math
import networkx as nx
import itertools
from . import sbml as sbml_module
from .utils import quoted, unquoted, solve, get_terminal_nodes, inverted_dag, remove_terminal
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


ASP_SRC_ENUM_CC = resource_filename(__name__, 'asp/enum-cc.lp')
ASP_SRC_SIMPLE_SEED_SOLVING = resource_filename(__name__, 'asp/simple-seed-solving.lp')
ASP_SRC_GREEDY_TARGET_SEED_SOLVING = resource_filename(__name__, 'asp/greedy-target-seed-solving.lp')
ASP_SRC_GREEDY_TARGET_SEED_SOLVING_SEED_MINIMALITY = resource_filename(__name__, 'asp/greedy-target-seed-solving-seed-minimality-constraint.lp')
ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM = resource_filename(__name__, 'asp/iterative-target-seed-solving--aim.lp')
ASP_SRC_PARETO_OPTIMIZATIONS = resource_filename(__name__, 'asp/pareto_utnu.lp')
import os
assert os.path.exists(ASP_SRC_ENUM_CC), ASP_SRC_ENUM_CC
assert os.path.exists(ASP_SRC_SIMPLE_SEED_SOLVING), ASP_SRC_SIMPLE_SEED_SOLVING


def opt_models_from_clyngor_answers(answers:iter, *, smaller_is_best:bool=True):
    "Return tuple of optimal models found by clingor.solve, when option '--opt-mode=optN' is given"
    best_opt, models = math.inf, []
    if not smaller_is_best:  best_opt *= -1
    for model, opt in answers.with_optimization:
        # print('OPT, MODEL:', opt[0], model)
        if (opt[0] < best_opt) if smaller_is_best else (opt[0] > best_opt):
            best_opt, models = opt[0], []  # model will be given again as last model, so no need to include it twice
        else:
            models.append(model)
    return tuple(models)


def search_seeds(graph_data:str, start_seeds:iter=(), forbidden_seeds:iter=(),
                 targets:set=(), graph_filename:str=None, enum_mode:str=EnumMode.Enumeration,
                 explore_pareto:bool=False, greedy:bool=False, **kwargs) -> [{set}]:
    "Wrapper around all seeds search methods. The used methods depends of given parameters."
    enum_mode = EnumMode(enum_mode)
    if not targets:  # no target, just activate everything
        func = search_seeds_activate_all
    elif explore_pareto:  # explore the pareto front
        func = search_pareto_front
    elif greedy:  # non efficient search of targets
        func = search_seeds_activate_targets_greedy  # TODO: this is not tested…
    else:  # efficient search of targets
        func = search_seeds_activate_targets_iterative
    yield from func(graph_data, frozenset(start_seeds), frozenset(forbidden_seeds), frozenset(targets), graph_filename, enum_mode, **kwargs)


def search_pareto_front(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration) -> [set]:
    """Yield the set of seeds for each found solution on the pareto front."""
    if enum_mode in {EnumMode.Union, EnumMode.Intersection}:
        raise ValueError(f"Mode {str(enum_mode)} is not supported for this routine.")
    # treat parameters
    if start_seeds & forbidden_seeds:
        raise ValueError(f"start_seeds and forbidden_seeds shares some seeds: {start_seeds & forbidden_seeds}.")
    targets = frozenset(map(quoted, targets))
    start_seeds = frozenset(map(quoted, start_seeds))
    forbidden_seeds = frozenset(map(quoted, forbidden_seeds))
    # data representation
    start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
    targets_repr = ' '.join(f'target({quoted(t)}).' for t in targets)
    forb_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
    data_repr = graph_data + start_seeds_repr + targets_repr + forb_repr

    files = ASP_SRC_PARETO_OPTIMIZATIONS, ASP_SRC_GREEDY_TARGET_SEED_SOLVING
    models = solve(files, inline=data_repr, clingo_bin_path='asprin').with_optimality.discard_quotes.by_predicate
    for model, opts, isoptimum in models:
        if isoptimum:
            yield frozenset(args[0] for args in model.get('seed', ()))


def search_seeds_activate_targets_iterative(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                                            graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration,
                                            compute_optimal_solutions:bool=False, filter_included_solutions:bool=True) -> [set]:
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

    """
    if enum_mode in {EnumMode.Union, EnumMode.Intersection}:
        raise ValueError(f"Mode {str(enum_mode)} is not supported for this routine.")
    print(start_seeds, forbidden_seeds, start_seeds & forbidden_seeds)
    if start_seeds & forbidden_seeds:
        raise ValueError(f"start_seeds and forbidden_seeds shares some seeds: {start_seeds & forbidden_seeds}.")

    def find_aim(scc_name:str, hypothesis:tuple) -> [str]:
        "Yield targets or reactions associated with given scc for given hypothesis"
        scc_targets = targets & sccs[scc_name]
        seeds, scc_reactions, fullfilled = hypothesis
        for scc, reactions in scc_reactions.items():
            for reaction, external_reactants in reactions.items():
                scc_targets |= external_reactants
        print('\tFIND AIM:', scc_name, hypothesis, '->', scc_targets)
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

    sccs, scc_dag = compute_sccs(graph_data, graph_filename=graph_filename)
    terminals = frozenset(get_terminal_nodes(scc_dag))
    rev_scc_dag = dict(inverted_dag(scc_dag))
    print('  SCC DAG:', scc_dag)
    # print('TERMINALS:', terminals)
    # print('         :', frozenset(rev_scc_dag[None]))
    assert terminals == frozenset(rev_scc_dag[None])
    # hypothesis = defaultdict(list)  # {scc: hypothesis} with hypothesis == iterable of (seeds, reactions)
    # associate an empty hypothesis for each SCC having a target.
    targets = frozenset(map(quoted, targets))
    for scc, nodes in sccs.items():
        print('SEARCHING TARGETS:', nodes, targets, nodes & targets)
        if nodes & targets:  # the terminal SCC has an aim
            all_hypothesis.append(get_null_hypothesis(scc))
    # iteratively find hypothesis
    while len(scc_dag) > 1:  # last valid key is None
        for terminal in frozenset(get_terminal_nodes(scc_dag)):
            print('TERMINAL:', terminal)
            if terminal is None: continue
            self_hypothesis = tuple(get_hypothesis_of(terminal))
            print('ALL HYPS:', all_hypothesis, '\t(does not contains TERM HYP)')
            print('TERM HYP:', self_hypothesis)
            print('    DAG :', scc_dag)
            for current_hypothesis in self_hypothesis:
                aim = find_aim(terminal, current_hypothesis)
                if not aim:  # the hypothesis has nothing to produce
                    seeds, scc_reactions, fullfilled = current_hypothesis
                    for parent in rev_scc_dag[terminal] or [None]:
                        new_scc_reactions = dict(scc_reactions)  # don't share among parents
                        new_scc_reactions[parent] = scc_reactions[terminal]
                        del new_scc_reactions[terminal]
                        new_hypothesis = set(seeds), new_scc_reactions, fullfilled & targets
                        all_hypothesis.append(new_hypothesis)
                        print(f"\tTerminal {terminal} has nothing to do ; propagation of hypothesis '{new_hypothesis}' to parent {parent}.")
                        break  # don't propagate such trivial hypothesis to all parents. Only one will have the job done
                    continue  # the current hypothesis has been treated
                start_seeds_repr = ' '.join(f'seed({quoted(s)}).' for s in start_seeds)
                targets_repr = ' '.join(f'target({quoted(t)}).' for t in aim)
                forb_repr = ' '.join(f'forbidden({quoted(s)}).' for s in forbidden_seeds)
                scc_repr = ' '.join(f'scc({terminal},{node}).' for node in sccs[terminal])
                scc_data = graph_data + f'current_scc({terminal}). {scc_repr} {start_seeds_repr} {forb_repr} {targets_repr}'
                print('\tCURRENT HYP:', current_hypothesis)
                print('\tCURR.  AIM :', aim)
                # print('\t  SCC DATA:', scc_data)
                for new_seeds, new_targets, new_fullfill in _compute_hypothesis_from_scc(terminal, scc_data, sccs, rev_scc_dag):
                    new_seeds |= current_hypothesis[0]
                    new_targets = {**current_hypothesis[1], **new_targets}
                    del new_targets[terminal]  # remove self from the hypothesis
                    new_fullfill |= current_hypothesis[2]
                    new_hypothesis = new_seeds, new_targets, new_fullfill & targets
                    print('\t\tADDED HYPOTHESIS:', new_hypothesis)
                    all_hypothesis.append(new_hypothesis)
            # now, remove terminal from the dag
            remove_terminal(terminal, scc_dag, frozenset(rev_scc_dag.get(terminal, ())))
    print('OUTPUT HYPOTHESIS:', all_hypothesis, targets, 'compute_optimal_solutions=', compute_optimal_solutions)
    return _solutions_from_hypothesis(all_hypothesis, targets, compute_optimal_solutions, filter_included_solutions)

def _solutions_from_hypothesis(all_hypothesis:list, targets:set, compute_optimal_solutions:bool=False, filter_included_solutions:bool=True) -> frozenset:
    """Compute the solutions from hypothesis and targets.
    If compute_optimal_solutions is given, will filter out all non-optimal solutions.

    """
    if not all_hypothesis:  # no hypothesis
        return frozenset()
    # verify that all targets are reachable
    reachables = set()
    for _, __, reachable in all_hypothesis:
        reachables |= targets & reachable
    if reachables != targets:  # not all are reachables
        # print('WARNING: Not all targets are reachables:', ', '.join(targets - reachables) + '. An error will be raised.')
        raise ValueError(f"Not all targets are reachables: {targets - reachables}.")
        # targets = reachables  # TODO: should raise an error, instead ?
    # quick data integrity check
    unfullfilled_targets = set(targets)
    for seeds, scc_reactions, fullfilled in all_hypothesis:
        assert len(scc_reactions) == 1, scc_reactions
        assert None in scc_reactions, scc_reactions
        assert len(scc_reactions[None]) == 0, scc_reactions[None]
        assert seeds, seeds
        assert fullfilled, fullfilled
        seeds = frozenset(map(unquoted, seeds))
        unfullfilled_targets -= fullfilled
        # print(f'SOLUTION: seeds {set(seeds)} are fullfilling {set(map(unquoted, fullfilled))}')
    assert not unfullfilled_targets, unfullfilled_targets  # some targets where not reachable ; what is this magic ?
    # collect, for each target, the set of seeds reaching it
    seeds_for_target = defaultdict(list)  # target: {seeds}
    for seeds, scc_reactions, fullfilled in all_hypothesis:
        for target in fullfilled:
            seeds_for_target[target].append(seeds)
    # generate all possible solutions
    # print('SEEDS FOR TARGETS:', seeds_for_target)
    # print('                 :', tuple(itertools.product(*seeds_for_target.items())))
    solutions = set()  # a solution is a set of seeds activating all targets
    for seeds_sets in itertools.product(*seeds_for_target.values()):
        # print('SSEEDS SETS:', seeds_sets)
        seeds_set = set.union(*seeds_sets)
        solutions.add(frozenset(map(lambda x: x.strip('"'), seeds_set)))
        print(f'OPT SOLUTION: seeds {seeds_set} are fullfilling {targets}')
    if filter_included_solutions:
        filtered_solutions = set()  # a solution is a set of seeds activating all targets
        for seeds_sets in solutions:
            if not any(other_set < seeds_sets for other_set in solutions):
                filtered_solutions.add(seeds_sets)
        solutions = filtered_solutions
        # print('FILTERED:', solutions)
    if compute_optimal_solutions:
        # reduce solutions using seed number
        opt_seed_number = min(map(len, solutions))
        solutions = (s for s in solutions if len(s) == opt_seed_number)
    return frozenset(solutions)


def _compute_hypothesis_from_scc(scc_name:str, scc_encoding:set, sccs:dict, rev_scc_dag:dict) -> [(set, dict, set)]:
    """Yield hypothesis computed from given scc_name to consider for next SCCs"""
    # the following call will provide us a model for each hypothesis.
    models = solve(ASP_SRC_ITERATIVE_TARGET_SEED_SOLVING__AIM, inline=scc_encoding, options='--opt-mode=optN')
    models = opt_models_from_clyngor_answers(models.by_predicate)
    print('\tCOMPUTE ALL HYPOTHESIS:', scc_name)
    print('                NB MODEL:', len(models))
    print('             SCC PARENTS:', rev_scc_dag[scc_name])
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
        print('\t\tFound Hypothesis:   targets:', new_targets, '\tseeds:', new_seeds)
        # create hypothesis with each parent SCC that have a reactant in it (or None for roots)
        scc_reactions = {None: dict(new_targets)}  # default case: no parent
        alien_reactants = frozenset(itertools.chain.from_iterable(new_targets.values()))
        print('\t\tSCC_REACTION BUILDING…')
        print('\t\t\t', rev_scc_dag[scc_name], new_targets)
        print('\t\t\t', alien_reactants)
        if rev_scc_dag[scc_name] and new_targets:
            scc_reactions = {parent: dict(new_targets) for parent in rev_scc_dag[scc_name]
                             if any(reactant in sccs[parent] for reactant in alien_reactants)}
        print('\t\t\t', scc_reactions)
        yield new_seeds, scc_reactions, new_fullfilled


def search_seeds_activate_targets_greedy(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(),
                                         graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration, compute_optimal_solutions:bool=False,
                                         filter_included_solutions:bool=True) -> [{set}]:
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
                   inline=data_repr, options='--opt-mode=optN ' + enum_mode.clingo_option).discard_quotes
    models = opt_models_from_clyngor_answers(models.by_predicate)
    for model in models:
        seeds = frozenset(args[0] for args in model['seed'] if len(args) == 1)
        yield seeds


def search_seeds_activate_all(graph_data:str, start_seeds:iter=(), forbidden_seeds:set=(), targets:set=(), graph_filename:str=None, enum_mode:EnumMode=EnumMode.Enumeration) -> [{set}]:
    """Yield the set of seeds for each found solution.

    This implements the most simple solution: no targets, find the minimum sets
    of seeds that activate everything.

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


def compute_sccs(graph_data:str, graph_filename:str=None) -> [{str}]:
    """Return the DAG and nodes of Strongly Connected Components found in given graph.

    graph_data -- ASP encoding of the metabolic network.
    graph_filename -- filename of SBML encoded graph, used to speed up SCC mining.
    return -- (SCC, DAG), with DAG a SCC id -> {SCC id sucessor} mapping,
              and SCC a SCC id -> {node} mapping.

    If no graph_filename given, ASP will be used to infer them from graph_data.

    """
    if not graph_filename or True:  # the networkx method is not yet implemented
        models = solve(ASP_SRC_ENUM_CC, inline=graph_data)
        for model in models.by_predicate:
            print('SCC MODEL:', model)
            roots = set(args[0] for args in model.get('noinput', ()) if len(args) == 1)
            sccs = defaultdict(set)  # SCC identifier: nodes in SCC
            for scc_name, node in model.get('scc', ()):
                sccs[scc_name].add(node)
            scc_dag = defaultdict(set)  # SCC identifier: successor SCCs
            for scc_name, scc_succ in model.get('sccedge', ()):
                scc_dag[scc_name].add(scc_succ)
            scc_dag[None] = roots
            sccs = dict(sccs)
            scc_dag = dict(scc_dag)
        print('SCC FINAL:', sccs)
        print()
        print('         :', scc_dag)
        print()
        return sccs, scc_dag
    else:  # use graph_filename
        raise NotImplementedError("Networkx based SCC extraction is not yet implemented")
        digraph = sbml_module.read_SBML_network_as_simple_graph(graph_filename)
        sccs = {min(nodes): frozenset(nodes) for nodes in nx.strongly_connected_components(digraph)}
        return sccs, scc_dag
