% Search of the minimal set of seeds or input reactions in given scc activating
%  all targets.
% It should be grounded with another encoding containing optimization criteria.
%  See iterative-target-seed-solving--aim--minimality-constraint.lp
%  for the minimality constraints to append to it for greedy search,
%  or pareto_iterative*.lp for pareto front search.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - target(T): node T must be activated
%   - scc(S,X): node X is in scc S
%   - current_scc(S): scc S is the scc to treat for this call
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
% OUTPUTS:  one model for each set of seed that activate all metabolites
%   - seed(S): node S is a seed
%   - new_target(T,R): node T is a reactant of reaction R, and a target that previous SCC will have to activate

% A metabolite is a reactant or product.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).

% A metabolite is activated if produced by a reaction with all its reactants activated.
activated(M) :- metabolite(M) ; scc(S,M) ; current_scc(S) ; product(M,R) ; activated(T): reactant(T,R).

% A seed is activated
activated(S) :- seed(S).

% A target in another SCC is activated.
activated(T) :- target(T) ; not scc(S,T) ; current_scc(S).

% Choose a set of seeds among metabolites of the current SCC.
{ seed(M): metabolite(M), not forbidden(M), scc(S,M), current_scc(S) }.

% Detect the set of input reactions that are feeding current SCC.
input_reaction(I):- reaction(I) ; current_scc(S) ; not scc(S,R) ; reactant(R,I) ; scc(S,M) ; product(M,I).

% Choose a set of input reactions.
{ input(I): input_reaction(I) }.

% Mark reactant of choosen input reactions as targets for future computations.
new_target(T,I) :- input(I) ; reactant(T,I) ; not scc(S,T) ; current_scc(S).
% Alien reactants of choosen input reactions are considered activated.
activated(M) :- input(I) ; reactant(M,I) ; not scc(S,M) ; current_scc(S).
% Mark targets that have been activated in that SCC.
activated_local_target(T) :- target(T) ; scc(S,T) ; current_scc(S).

% Discard any model that is not activating all targets.
:- target(M) ; not activated(M).
:- new_target(M) ; not activated(M).
