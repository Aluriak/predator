% Search of the minimal set of seeds in one SCC. Other SCC may be provided,
%  but are not necessary.
% Nodes that are product of foreign reactions are considered activated.
%  This encodes the fact that previous SCC are completely activated.
% This is a standalone/proof of concept program. See pareto_utnu*
%  for implementations to ground with standard non-optimized search.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - scc(S,X): node X is in scc S
%   - current_scc(S): scc S is the scc to treat for this call
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
% OUTPUTS:  one model for each set of seed that activate all metabolites
%   - seed(S): node S is a seed

% A metabolite is a reactant or product belonging to the current SCC.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).

%target(M) :- metabolite(M).

% A metabolite is activated if produced by a reaction with all its reactants activated.
activated(M) :- metabolite(M) ; product(M,R) ; activated(T): reactant(T,R).

% A seed is activated
activated(S) :- seed(S).

% Seeds among targets can happen
seed_target(M) :- seed(M), target(M).

% Choose a set of seeds among metabolites of the current SCC.
{ seed(M): metabolite(M), not forbidden(M) }.

% Number of seeds
nb(seed,N) :- N={seed(_)}.
nb(activated,N) :- N={activated(_)}.
nb(target,N) :- N={target(_)}.
% #minimize{N:nb_seed(N)}.

:- target(M), not activated(M).


% Preferences
% maximize size of scope
#preference(scope,more(cardinality)){activated(N)}.
% minimize number of seeds
#preference(seeds,less(cardinality)){seed(N)}.
% minimize number of targets among seeds
#preference(targets_as_seeds,less(cardinality)){seed(N), target(M)}.
% Pareto
#preference(all,pareto){**scope; **seeds; **targets_as_seeds}.

#optimize(all).
#show seed/1.
#show nb/2.
