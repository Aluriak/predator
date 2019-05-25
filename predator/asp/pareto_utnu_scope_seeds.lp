% Asprin definitions to ground with greedy seeds search to explore the pareto
%  fronts of the two metrics: number of seeds and size of scope.

% Number of seeds
nb(seed,N) :- N={seed(_)}.
nb(activated,N) :- N={activated(_)}.

% Preferences
% maximize size of scope
#preference(scope,more(cardinality)){activated(N)}.
% minimize number of seeds
#preference(seeds,less(cardinality)){seed(N)}.
% Pareto
#preference(all,pareto){**scope; **seeds}.

#optimize(all).
#show seed/1.
#show nb/2.
