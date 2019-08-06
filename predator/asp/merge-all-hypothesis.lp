% Tentative to implement a hypothesis merger.
% Instead of treating each hypothesis alone, it consider them in a single shot.
% This code needs to be grounded with predator/asp/iterative-target-seed-solving--aim.lp
%  (and its optimization method, either pareto or seed minimization).
% Each model is a set of hypothesis used to obtain the output, and therefore to populate seeds and new_targets,
%  along with the standard model output (seeds, new_target,…),
%  constituting an hypothesis for parent SCC.
% In the absence of hypothesis (typically the starting terminal SCCs), this code has no effect.
% If all hypothesis are exclusive in term of optimization, only models using exactly one of them will be yielded, preserving the standard behavior.
% INPUTS:
%   - seed(S): node S is a seed
%   - forbidden(S): node S cannot be a seed
%   - target(T): node T must be activated
%   - scc(S,X): node X is in scc S
%   - current_scc(S): scc S is the scc to treat for this call
%   - reaction(R): R is a reaction.
%   - reactant(T,R): T is a reactant of reaction R.
%   - product(P,R): P is a product of reaction R.
%   - hypothesis(H),hyp_seed(H,S),hyp_new_target(H,NT,R),hyp_fullfilled(H,F):
%       Hypothesis H fullfills targets F using seeds S and new_targets NT in reactions R.
% SUPPLEMENTARY OUTPUTS:
%   - choosen_hypothesis(H): H is an hypothesis to use. Of no use, except for debugging.


% Take any number of hypothesis.
{choosen_hypothesis(H): hypothesis(H)}.
% NOTE: we need to minimize the number of hypothesis with low priority,
%  in order to avoid equivalent models with one using a superset of hypothesis of another.
nb(hypothesis,X) :- X={choosen_hypothesis(_)}.
#minimize{X,5:nb(hypothesis,X)}.

% Consequences of the hypothesis.
seed(X) :- choosen_hypothesis(H) ; hyp_seed(H,X).
new_target(NT,R) :- choosen_hypothesis(H) ; hyp_new_target(H,NT,R).

% All fullfilled targets must be activated.
%  (don't choose a subset of hypothesis loosing anything)
activated(X) :- choosen_hypothesis(H) ; hyp_fullfilled(H,X).
:- hyp_fullfilled(_,F) ; not activated(F).
