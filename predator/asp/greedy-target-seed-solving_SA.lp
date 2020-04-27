% Search of the minimal set of seeds in all graph activating all targets.
%  This is a greedy method, that is not expected to scale.
% It should be grounded with another encoding containing optimization criteria,
%   such as greedy-target-seed-solving-seed-minimality-constraint.lp for simple search,
%   or pareto_utnu*.lp for pareto exploration with the asprin system.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - target(T): node T must be present
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
% OUTPUTS:  one model for each set of seed that activate all metabolites
%   - seed(S): node S is a seed

%%%%%%% Common in all versions %%%%%%%%%%
% there can be existing seeds
% A seed is present
any_seed(S) :- existing_seed(S).
any_seed(S) :- seed(S).

% A metabolite is a reactant or product.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).


% a reaction is feasible if all reactants are present
feasible(R) :- present(M) : reactant(M,R); reaction(R).

% a seed is present
present(M) :- any_seed(M).

% % every target must be present
% :- target(X), not present(X).
% located in the constraint/minimisation file

%%%%%%% version 1 %%%%%%%%

% % a metabolite is present or not
% { present(M) } :- metabolite(M).

% products of a feasible reaction are producible
producible(M) :- feasible(R), product(M,R).

% % every present metabolite must be producible
% :- present(M), not producible(M), not any_seed(M).

%%%%%%% version 2 %%%%%%%%


% try to cut the set of possibly present metabolites
reach(X) :- any_seed(X).
reach(X) :- product(X,R), reactant(Y,R), reach(Y).

% a metabolite might be present or not
{ present(X) } :- product(X,R), reach(X). 

:- present(X), not producible(X).

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% #show seed/1.
% located in the constraint/minimisation file
