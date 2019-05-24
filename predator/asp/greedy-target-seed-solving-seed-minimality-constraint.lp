% Encoding of the minimality constraint for greedy target seed solving.

% Ensure its minimality
nb_seed(N) :- N={seed(_)}.
#minimize{N:nb_seed(N)}.
