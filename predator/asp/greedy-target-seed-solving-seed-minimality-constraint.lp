% Encoding of the minimality constraint for greedy target seed solving.

% choose a set of seed from the possible seeds
{ seed(M): can_reach(M), not forbidden(M) }.

% Discard any model that is not activating all targets.
:- target(M) ; not activated(M).

% Ensure its minimality
nb_seed(N) :- N={seed(_)}.
#minimize{N:nb_seed(N)}.

#show seed/1.
#show nb_seed/1.
