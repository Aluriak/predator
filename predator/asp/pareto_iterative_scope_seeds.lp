% Asprin definitions to ground with iterative seeds search to explore the pareto
%  fronts of the two metrics: number of seeds and size of scope.

% Number of seeds
nb(seed,N) :- N={seed(_)}.
nb(new_target,N) :- N={new_target(_,_)}.
nb(seed_or_new_target,N+M) :- nb(seed,N) ; nb(new_target,M).
nb(activated,N) :- N={activated(_)}.
nb(target,N) :- N={target(_)}.

seed_or_new_target(X) :- seed(X).
seed_or_new_target(X) :- new_target(X,_).

% Preferences
% maximize size of scope
#preference(scope,more(cardinality)){activated(N)}.
% minimize number of seeds
#preference(seeds,less(cardinality)){seed_or_new_target(N)}.
% Pareto
#preference(all,pareto){**scope; **seeds}.

#optimize(all).
