% Encoding of the minimality constraint for iterative target seed solving.

% Ensure its minimality
nb_seed(N) :- N={seed(_)}.
nb_input(N) :- N={input(_)}.
#minimize{S+I:nb_seed(S), nb_input(I)}.
