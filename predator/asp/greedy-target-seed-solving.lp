% Search of the minimal set of seeds in all graph activating all targets.
%  This is a greedy method, that is not expected to scale.
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
activated(S) :- seed(S).

% Choose a set of seeds among metabolites of the current SCC.
{ seed(M): metabolite(M), not forbidden(M) }.

% Ensure its minimality
nb_seed(N) :- N={seed(_)}.
#minimize{N:nb_seed(N)}.

% Discard any model that is not activating all targets.
:- target(M) ; not activated(M).
