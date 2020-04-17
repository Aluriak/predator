% Search of the minimal set of seeds in all graph activating all targets.
%  This is a greedy method, that is not expected to scale.
% It should be grounded with another encoding containing optimization criteria,
%   such as greedy-target-seed-solving-seed-minimality-constraint.lp for simple search,
%   or pareto_utnu*.lp for pareto exploration with the asprin system.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - target(T): node T must be activated
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
% OUTPUTS:  one model for each set of seed that activate all metabolites
%   - seed(S): node S is a seed

% A metabolite is a reactant or product.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).

% A metabolite is activated if produced by a reaction with all its reactants activated.
activated(M) :- metabolite(M) ; product(M,R) ; activated(T): reactant(T,R).

% A seed is activated
activated(S) :- existing_seed(S).
activated(S) :- seed(S).

% try to cut the set of possible seeds
can_reach(M) :- target(M).
can_reach(M) :- reactant(M,R) ; product(P,R) ; can_reach(P).

% % Choose a set of seeds among metabolites that can reach the targets somehow.
% { seed(M): can_reach(M), not forbidden(M) }.
% located in the constraint/minimisation file

% % Discard any model that is not activating all targets.
% :- target(M) ; not activated(M).
% located in the constraint/minimisation file

% nb(seed,N) :- N={seed(_)}.
% nb(activated,N) :- N={activated(_)}.
% nb(target,N) :- N={target(_)}.
% nb(canreach,N) :- N={can_reach(_)}.

% #show seed/1.
% located in the constraint/minimisation file