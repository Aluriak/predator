% Number of seeds
nb(seed,N) :- N={seed(_)}.
nb(activated,N) :- N={activated(_)}.
nb(target,N) :- N={target(_)}.
% #minimize{N:nb_seed(N)}.

:- target(M), not activated(M).

% Preferences
% maximize size of scope
#preference(scope,more(cardinality)){activated(N)}.
% minimize number of seeds
#preference(seeds,less(cardinality)){seed(N)}.
% minimize number of targets among seeds
#preference(targets_as_seeds,less(cardinality)){seed(N), target(M)}.
% Pareto
#preference(all,pareto){**scope; **seeds; **targets_as_seeds}.

#optimize(all).
#show seed/1.
#show nb/2.