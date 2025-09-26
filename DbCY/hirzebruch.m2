debug needsPackage "CoxCategories"
--needs "~/utah/DbCY/GlobalExt.m2"
debug Varieties

X = hirzebruchSurface 2
S = ring X

OC = koszulComplex (vars S)_{1} -- TODO: why is this only spherical on H_2?
L = dual \ reverse { OO_X^{{-1,-1}}, OO_X^{{0,-1}}, OO_X^{{-1,0}}, OO_X^1 }

-- OC is spherical:
-- assert(hilbertPolynomial(sheaf OC ** OO toricDivisor X) == hilbertPolynomial(sheaf OC))
-- ExtTable(X, {OC}) -- should be 1+T^(dim X)

--
Y = toricProjectiveSpace 5
phi = map(Y, X, matrix {{1, 0}, {0, 1}, {1, 1}, {2, 1}, {3, 1}})
R = quotient ideal phi
Z = Proj R

OC' = complex coker (vars R)_{5,4,3,2}
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
C[2], source ev
D == target ev
D , cone ev
-- FIXME
ev' = derivedCoevaluationMap(C, D)
C, source ev'
D[-2], target ev'

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
X = toricProjectiveSpace 2
S = ring X
BT = complex \ zonotopeBundles X
degs = zonotopeDegrees S
sheaf freeResolution simple_S degs#2
sheaf freeResolution simple_S degs#3
sheaf (koszulComplex (vars S)_{1,3})

sheaf rightMutation(complex S^1, complex S^{{3,-2}})

globalExt = ToricExt
BT = sheaf \ dual \ reverse BT
L = leftOrthogonal BT
ExtTable_X mutate(1, BT)
ExtTable_X mutate(0, mutate(1, BT)) -- Note, this changes!
ExtTable_X mutate(1, mutate(0, mutate(1, BT)))
mutate(1, mutate(2, mutate(0, mutate(1, BT))))

-- FIXME: some terms are truncated too much
L = rightOrthogonal BT
L = mutate(BT, 0)
globalExt(L#1, L#2)
rightMutation(L#1, L#2)
mutate(mutate(mutate(BT, 0), 1), 0)
mutate(1, mutate(0, mutate(1, BT)))
OO toricDivisor X



L = elapsedTime leftOrthogonal BT
netList apply(6, i -> {prune HH sheaf minimize L#i, prune HH sheaf minimize freeResolution simple_S degs#(-i-1)})
L#0
(prune minimize L#0)


E = BT; apply({4,3,2,1,0}, i -> E = mutate(i, E))
netList apply(6, i -> sheaf E#i)
