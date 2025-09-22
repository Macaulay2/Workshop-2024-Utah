debug needsPackage "CoxCategories"
--needs "~/utah/DbCY/GlobalExt.m2"
debug Varieties

P2 = toricProjectiveSpace 2
X = toricBlowup({0,1}, P2)
S = ring X
p = X^[]

-- Orlov's formula for the blow up
L = { p^* OO_P2^{0}, p^* OO_P2^{1}, p^* OO_P2^{2}, sheaf(S^1/S_3) }
elapsedTime ExtTable(X, complex \ L) -- upper triangular means exceptional
-- now compare with complexes version
--assert(0 == ExtTable(X, L) - ExtTable'(X, apply(L, M -> freeResolution module M)))

--
-- (L1,L2) = (OO_X^1, sheaf(S^1/S_3))
-- (C,D) = freeResolution \ module \ (L1,L2)
-- (F,G) = module \ sheaf \ (C,D)

-- ExtTable(X, {F,G})
-- ExtTable(X, {L1,L2})

-- globalExt(0,C,D), Ext^0(L1, L2)
-- globalExt(1,C,D), Ext^1(L1, L2)
-- globalExt(2,C,D), Ext^2(L1, L2)

Y = toricProjectiveSpace 4
phi = map(Y, X, matrix {{1, 0}, {0, 1}, {1, 1}, {2, 1}})
R = quotient ideal phi
Z = Proj R

-- TODO: what are the pushforwards of L in P^4?
L

L' = {
    coker map(R^{{1},4:{0}}, ,{{0, 0, x_1*x_4, 0, x_0*x_4, x_0*x_4, 0, x_1^2, x_2*x_4, x_2*x_4, 0, 0, x_2*x_3, 0, x_2^2}, {0, x_3, 0, x_2, 0, 0, x_0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -x_4, -x_3, -x_3, -x_2, -x_2, -x_1, -x_0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -x_1, -x_1, x_3, x_2, 0, x_0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x_4, -x_3, -x_1, -x_1, -x_0}}),
    coker map(R^6,R^{14:{-1}, {0}},{{x_3, x_0, x_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-x_4, -x_1, 0, x_3, x_2, x_0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -x_4, -x_4, -x_3, -x_1, -x_4, -x_3, 0, 0, -x_3, -x_2, 0, 0, 0}, {0, 0, 0, 0, 0, 0, x_1, x_0, x_3, x_2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -x_4, -x_3, x_1, x_0, -x_3, -x_2, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x_1, x_0, 0}}),
    coker map(R^3,R^{6:{-1}},{{x_3, x_0, x_2, 0, 0, 0}, {-x_4, -x_1, 0, x_3, x_2, x_0}, {0, 0, -x_4, -x_4, -x_3, -x_1}}),
    coker map(R^3,R^{5:{-1}, {0}},{{x_3, x_2, x_0, 0, 0, 0}, {-x_4, -x_3, -x_1, -x_3, -x_2, 0}, {0, 0, 0, x_1, x_0, 0}}),
    coker map(R^3,R^{5:{-1}, {0}},{{x_3, x_2, x_0, 0, 0, 0}, {-x_4, -x_3, -x_1, -x_3, -x_2, 0}, {0, 0, 0, x_1, x_0, 0}}),
    coker map(R^5,R^{12:{-1}},{{x_3, x_0, x_2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-x_4, -x_1, 0, x_3, x_2, x_0, 0, 0, 0, 0, 0, 0}, {0, 0, -x_4, -x_4, -x_3, -x_1, -x_4, -x_3, 0, 0, -x_3, -x_2}, {0, 0, 0, 0, 0, 0, x_1, x_0, x_3, x_2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -x_4, -x_3, x_1, x_0}}),
    coker map(R^{{1},2:{0}}, ,{{0, 0, 0, x_1*x_4, 0, x_0*x_4, x_0*x_4, 0, x_1^2}, {0, 0, x_3, 0, x_2, 0, 0, x_0, 0}, {0, 0, -x_4, -x_3, -x_3, -x_2, -x_2, -x_1, -x_0}}),
    coker map(R^2,R^{3:{-1}},{{x_3, x_2, x_0}, {-x_4, -x_3, -x_1}})
    }

L
netList apply(subsets(5, 2), ell -> hilbertPolynomial sheaf pullback_phi flattenModule coker (vars R)_ell)
hilbertPolynomial L#3

apply({R^1, }, E -> prune sheaf pullback_phi flattenModule E)

netList apply(L', E -> hilbertPolynomial(pullback_phi flattenModule(E ** R^{-2}), Projective => false)), netList apply(L, hilbertPolynomial)
hilbertPolynomial(pullback_phi flattenModule(L'#0 ** R^{-2}), Projective => false), netList apply(L, hilbertPolynomial)

L' =
debug needsPackage "DirectSummands"
errorDepth=2
summands pushForward'(inducedMap phi, module L#0, options pushForward)

L' = dual \ reverse {
    R^{{-1}},
    image map(R^4, R^{2:{-1}}, {{x_5, x_4}, {x_4, x_3}, {x_3, x_2}, {x_1, x_0}}),
    image map(R^2, R^{4:{-1}}, {{x_5, x_4, x_3, x_1}, {x_4, x_3, x_2, x_0}}),
    R^1
    }

-- CoherentSheaf == CoherentSheaf := Boolean => (F, G) -> (
--     X := variety F;
--     hilbertPolynomial_X F === hilbertPolynomial_X G
--     and module prune F == module prune G)
assert(sheaf HH_0 OC == sheaf pullback_phi flattenModule OC'_0)
assert(L == apply(L', F -> prune sheaf pullback_phi flattenModule F))

----------
-- Orlov's blow up formula

psi^*(P112_0)

end--
restart
needs "hirzebruch.m2"

debugLevel=1
errorDepth=2
netList_2 L, netList_2 apply(L', sheaf)
elapsedTime ExtTable(X, complex \ L) -- ~6s
elapsedTime ExtTable(Z, complex \ L') -- ~6s

(C, D) = (sheaf OC', sheaf complex L'_0)
globalExt(C, D)
leftMutation(C, D)
ev = derivedEvaluationMap(C, D)
D , cone ev

netList_2 apply(L', E -> globalExt(sheaf OC', sheaf complex E))
L2' = apply(L', E -> leftMutation(sheaf OC', sheaf complex E))
-- elapsedTime ExtTable(X, L2')
netList L2'
netList apply(L2', E -> prune HH E)
T = schedule(() -> elapsedTime ExtTable'(X, L2')) -- verrry slow

cachedExtTable(X, L2')

-- L2' = apply(L', E -> rightMutation(sheaf complex E, sheaf OC'))
-- netList apply(L2', E -> prune HH E)
-- elapsedTime ExtTable'(X, L2) -- exceptional!

-- only works on my system for now ...
-- needs "GlobalExt.m2"
-- ExtTable = (X, L) -> (
--     T := (degreesRing 1)_0;
--     V := table(L, L, globalExt);
--     matrix applyTable(V, rankPolynomial_T))
-- elapsedTime ExtTable(X, complex \ L_{0,1,2,3})
-- | 1 0 0 0 |
-- | 2 1 0 0 |
-- | 4 2 1 0 |
-- | 6 4 2 1 |

L2 = apply(L, E -> leftMutation(OC, complex E))
netList apply(L2, E -> prune HH E)
ExtTable'(X, L2) -- exceptional!




---
restart
debug needsPackage "CoxCategories"
needs "GlobalExt.m2"
needs "BTtrunc.m2"
X = hirzebruchSurface 3
S = ring X
BT = complex \ zonotopeBundles X
degs = zonotopeDegrees S
sheaf freeResolution simple_S degs#2
sheaf freeResolution simple_S degs#3
sheaf (koszulComplex (vars S)_{1,3})

sheaf rightMutation(complex S^1, complex S^{{3,-2}})

globalExt = ToricExt
sheaf \ BT
sheaf \ mutate(0, mutate(1, BT))
sheaf \ mutate(1, mutate(2, mutate(0, mutate(1, BT))))



L = elapsedTime leftOrthogonal BT
netList apply(6, i -> {prune HH sheaf minimize L#i, prune HH sheaf minimize freeResolution simple_S degs#(-i-1)})
L#0
(prune minimize L#0)


E = BT; apply({4,3,2,1,0}, i -> E = mutate(i, E))
netList apply(6, i -> sheaf E#i)
