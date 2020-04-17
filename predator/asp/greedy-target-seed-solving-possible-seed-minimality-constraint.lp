% Encoding of the minimality constraint for greedy target seed solving when a set of possible seeds is given as input.

% choose a set of seed from the possible seeds
{ seed(M): possible_seed(M), not forbidden(M) }.

producible_target(M) :- target(M) ; activated(M).
unproducible_target(M) :- target(M) ; not activated(M).


% count producible targets and seed
nb_seed(N) :- N={seed(_)}.
nb_prod_target(N) :- N={producible_target(_)}.
nb_unprod_target(N) :- N={unproducible_target(_)}.

#minimize{N@2:nb_unprod_target(N)}.
#minimize{N@1:nb_seed(N)}.

#show seed/1.
#show nb_seed/1.
#show nb_prod_target/1.
#show nb_unprod_target/1.
#show producible_target/1.
#show unproducible_target/1.