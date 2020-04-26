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

%%%%%%% Common in all versions %%%%%%%%%%
% there can be existing seeds
% A seed is activated
any_seed(S) :- existing_seed(S).
any_seed(S) :- seed(S).

% A metabolite is a reactant or product.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).

% try to cut the set of possible seeds
can_reach(M) :- target(M).
can_reach(M) :- reactant(M,R), product(P,R), can_reach(P).

% % choose a set of seed from the possible seeds
% { seed(M): can_reach(M), not forbidden(M) }.
% located in the constraint/minimisation file

% a reaction is applied if all reactants are activated
feasible(R) :- activated(M) : reactant(M,R); reaction(R).

% a seed is activated
activated(M) :- any_seed(M).

% % every target must be activated
% :- target(X), not activated(X).
% located in the constraint/minimisation file

%%%%%%% version 1 %%%%%%%%

% % a metabolite is activated or not
% { activated(M) } :- metabolite(M).

% % products of a feasible reaction are producible
% producible(M) :- feasible(R), product(M,R).

% % every activated metabolite must be producible
% :- activated(M), not producible(M), not any_seed(M).

%%%%%%% version 2 %%%%%%%%

% a metabolite is activated if it is given as a seed or if it is the product of an applied reaction
activated(X) :- product(X,R), feasible(R).

reach(X) :- any_seed(X).
reach(X) :- product(X,R), reactant(Y,R), reach(Y).

% a metabolite might be activated or not
{ activated(X) } :- product(X,R), reach(X). 

% atleast one reaction must be producing X
is_produced(X) :- feasible(R), product(X,R).
:- not is_produced(X), activated(X), not any_seed(X).

% there must be a path from a seed to the activated metabolite
founded(X) :- any_seed(X).
founded(X) :- product(X,R), reactant(Y,R), founded(Y), feasible(R).

:- activated(X), not founded(X).

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% #show seed/1.
% located in the constraint/minimisation file
