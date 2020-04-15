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

% A metabolite is a reactant or product.
metabolite(M) :- reactant(M,_).
metabolite(M) :- product(M,_).

% try to cut the set of possible seeds
can_reach(M) :- target(M).
can_reach(M) :- reactant(M,R) ; product(P,R) ; can_reach(P).

% choose a set of seed from the possible seeds
{ seed(M): can_reach(M), not forbidden(M) }.

%%%%%%% version 1 %%%%%%%%

% a seed is present
present(M) :- seed(M).

% a metabolite is present or not
{ present(M) } :- metabolite(M).

% a reaction is feasible if all reactants are present
feasible(R) :- present(M) : reactant(M,R); reaction(R).

% products of a feasible reaction are producible
producible(M) :- feasible(R), product(M,R).

% every present metabolite must be producible
:- present(M), not producible(M), not seed(M).

% every target must be present
:- target(X), not present(X).

%%%%%%% version 2 %%%%%%%%

% % a reaction is applied if all reactants are available
% areaction(R) :- available(M) : reactant(M,R); reaction(R).
% % otherwise its inactive
% ireaction(R) :- reaction(R), not areaction(R).

% % a metabolite is available if it is given as a seed or if it is the product of an applied reaction
% available(X) :- seed(X).
% available(X) :- product(X,R), areaction(R).

% reach(X) :- seed(X).
% reach(X) :- product(X,R), reactant(Y,R), reach(Y).

% % a metabolite might be available or not
% { available(X) } :- product(X,R), reach(X). 

% % atleast one reaction must be producing X
% % if all are inactive it is not a solution 
% is_produced(X) :- areaction(R), product(X,R).
% :- not is_produced(X), available(X), not seed(X).

% % there must be a path from a seed to the available metabolite
% founded(X) :- seed(X).
% founded(X) :- product(X,R), reactant(Y,R), founded(Y), areaction(R).

% :- available(X), not founded(X).

% % every target must be available
% :- target(X), not available(X).

%%%%%%%%%%%%%%%%%%%%%%%%%%%


#show seed/1.
