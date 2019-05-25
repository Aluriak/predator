% Asprin definitions to ground with greedy seeds search to explore the pareto
%  fronts of the three metrics: number of seeds, size of scope,
%  number of targets that are also seeds.

% Number of seeds
nb(seed,N) :- N={seed(_)}.
nb(activated,N) :- N={activated(_)}.
nb(target,N) :- N={target(_)}.

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
