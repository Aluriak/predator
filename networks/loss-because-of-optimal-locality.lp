% See unit test in test/test_predator.py:test_loss_because_of_optimal_locality for a complete
%  testing of that file.
% SBML/XML version in networks/loss-because-of-optimal-locality.sbml
% NetSeed will find seeds A, B and C, but not the combinations of GIJ and HKL.

reactant(X,R):- reaction(X,R,_).
product(X,R) :- reaction(_,R,X).
reaction(R) :- reaction(_,R,_).

reaction((a;b;c),r1,(d;e;f)).
reaction((d;e;f),r1r,(a;b;c)).
reaction((d;e;f),r2,(g;h)).
reaction(g,r3,i).
reaction(i,r4,j).
reaction(j,r5,g).
reaction(h,r6,k).
reaction(k,r7,l).
reaction(l,r8,h).
reaction((j;k),r9,m).


% Variation upon the network:
% reaction(z,r0,(a;b;c)).

% With that variation, the results suddenly change in a weird way.
