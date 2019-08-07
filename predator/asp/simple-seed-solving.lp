% Search of the minimal set of seeds in one SCC. Other SCC may be provided,
%  but are not necessary.
% Nodes that are product of foreign reactions are considered activated.
%  This encodes the fact that previous SCC are completely activated.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - scc(S,X): node X is in scc S
%   - current_scc(S): scc S is the scc to treat for this call
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
%   - metabolite(M): M is a metabolite of the graph.
% OUTPUTS:  one model for each set of seed that activate all metabolites
%   - seed(S): node S is a seed

% A metabolite is a reactant or product belonging to the current SCC.
scc_metabolite(M) :- reactant(M,_) ; current_scc(S) ; scc(S,M).
scc_metabolite(M) :- product(M,_) ; current_scc(S) ; scc(S,M).
scc_metabolite(M) :- metabolite(M) ; current_scc(S) ; scc(S,M).

% All foreign reactants of (partially) foreign reactions are activated.
activated(T) :- not scc_metabolite(T) ; reaction(R) ; reactant(T,R) ; scc_metabolite(M) ; product(M,R).

% A metabolite is activated if produced by a reaction with all its reactants activated.
activated(M) :- scc_metabolite(M) ; product(M,R) ; activated(T): reactant(T,R).

% A seed is activated
activated(S) :- seed(S).

% Choose a set of seeds among metabolites of the current SCC.
{ seed(M): scc_metabolite(M), not forbidden(M) }.

% Ensure its minimality
nb_seed(N) :- N={seed(_)}.
#minimize{N:nb_seed(N)}.

% Discard any model that is not activating all metabolites of the SCC.
:- scc_metabolite(M) ; not activated(M).
