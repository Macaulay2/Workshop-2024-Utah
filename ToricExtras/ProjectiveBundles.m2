--needs "helpers.m2"
needsPackage "NormalToricVarieties"

importFrom_Core "concatCols"
importFrom_NormalToricVarieties "cartierCoefficients"

-----------------------------------------------------------------------------
-- helpers

vecs = m -> entries transpose m

-- given a list return a permutation w such that L_w is increasingly ordered
-- e.g. perm {2,1,4} gives {1,0,2}
--perm = L -> last \ sort(reverse \ toList pairs L)

-- given a permutation and a list, permute elements of L according to w
-- e.g. move_{1,0,2} {2,0,1} gives {2,1,0}
move = (w, L) -> apply(L, i -> position(w, j -> i == j))

-- given the Cox ring of a projective bundle, give the indices of
-- rays corresponding to the base and the fiber of the bundle
-- TODO: this is not very robust and depends on a particular basis of Pic X
-- it should instead be determined by selecting the rays in the fan
baseVars  = S -> select(numgens S, i -> last degree S_i == 0)
fiberVars = S -> select(numgens S, i -> last degree S_i == 1)

-- ad-hoc check for two toric varieties being the same up to permutation of the rays
-- compare X and Y, possibly with a permutation of the coordinates given by perm
sameToricVariety  = (X, Y) -> X === Y or fan X == fan Y or any(permutations(d := dim X), perm -> sameToricVariety'(X, Y, id_(ZZ^d)_perm))
sameToricVariety' = (X, Y, m) -> try (f := inducedMap map(X, Y, m)) then f ideal X == ideal Y else false

-----------------------------------------------------------------------------
-- Projectivization of a sum of line bundles on a toric variety

protect Fiber
protect ProjectiveBundle

-- see CLS 7.3
PP = method()
PP CoherentSheaf := E -> E.cache.ProjectiveBundle ??= (
    X := variety E;
    r := rank E - 1; -- a rank r+1 sheaf gives a PP^r bundle on X
    if not isFreeModule module E then error "projectivization of arbitrary sheaves if not yet implemented";
    if r == 0 then return X;
    -- see CLS pp. 337
    p := inverse fromWDivToCl X * fromPicToCl X;
    -- TODO: should the degrees be sorted?
    a := apply(-degrees E, deg -> first vecs(p * transpose matrix {deg}));
    --
    d := dim X;
    N0 := ZZ^d; -- lattice of the base X
    N1 := ZZ^r; -- lattice of the fiber PP^r
    i0 := (id_N0 ++ 0 * id_N1)_{0 ..< d};
    i1 := (0 * id_N0 ++ id_N1)_{d ..< d+r};
    B := transpose matrix rays X;
    -- see CLS (7.3.2)
    Facets := affineImage_i1 \ facesAsCones(0, fan toricProjectiveSpace r);
    -- see CLS (7.3.3)
    Sigmas := apply(#rays X,
	rho -> i0 * B_{rho} + i1 * sum toList apply(1 .. r, k -> (a_k_rho - a_0_rho) * N1_{k-1}));
    -- This one is kind of out of order
    Y0 := normalToricVariety fan flatten table(max X, r + 1,
	(ell, i) -> coneFromVData(concatCols Sigmas_ell) + Facets_i);
    -- TODO: this part doesn't _always_ work
    eff0 := inverse nefGenerators Y0 * effGenerators Y0;
    -- FIXME: the columns are not sorted the way one might expect
    perm := sortColumns eff0;
    Y := normalToricVariety((rays Y0)_perm, move_perm \ max Y0,
	CoefficientRing => X.cache.CoefficientRing,
	Variable        => X.cache.Variable,
	WeilToClass     => eff0_perm);
    Y.cache.Base = X;
    Y.cache.Fiber = E;
    Y)


-- TODO: this should project the fan in some way
-- TODO: if X is not a fiber bundle, should this be an error instead?
base NormalToricVariety  := X -> if X.cache.?Base  then X.cache.Base  else X
-- TODO: technically not always correct
fiber = method()
fiber NormalToricVariety := X -> if X.cache.?Fiber then X.cache.Fiber else X.cache.CoefficientRing^(dim X + 1)

-- TODO:
-- give the toric bundle map p^*E^* --> OO_X(1), or similar

-- TODO:
-- give the exact sequence from CLS Thm 8.1.6

--- for personal tests only, to be removed or modified
kk = ZZ/101
PP ZZ   := NormalToricVariety => memoize(n -> toricProjectiveSpace(n, CoefficientRing => kk))
PP List := NormalToricVariety => memoize(w -> weightedProjectiveSpace(w, CoefficientRing => kk))

end--
restart
needs "ProjectiveBundles.m2"

-- PP^1
X = PP 1
assert isWellDefined X; nefGenerators X, effGenerators X

-- Hirzebruch of type 3
Y = PP(OO_X^1 ++ OO_X^{3})
assert isWellDefined Y; nefGenerators Y, effGenerators Y

-- A PP^1 bundle over Y
Z = PP OO_Y^2
assert isWellDefined Z; nefGenerators Z, effGenerators Z

-- Z has the same fan as HH_3 ** PP^1
assert sameToricVariety(Z, Y ** X)
assert sameToricVariety(Z, X ** Y)
-- finding the permutation of rays that induces an isomorphism
select(permutations 3, perm -> sameToricVariety'(Z, X ** Y, id_(ZZ^3)_perm))

-- A PP^1 bundle over Z
W = PP(OO_Z^1 ++ OO_Z^{{3,0,0}})
assert isWellDefined W; nefGenerators W, effGenerators W

-- W has the same fan as HH_3 ** HH_3
-- FIXME: should this be true?
sameToricVariety(W, Y ** Y)

--
restart
needs "ProjectiveBundles.m2"

-- PP^1
X = PP 1
Y = PP(OO_X^1 ++ OO_X^{3})
Z = PP(OO_Y^1 ++ OO_Y^{{2,4}})
assert isWellDefined Z; nefGenerators Z, effGenerators Z

---
restart
needs "ProjectiveBundles.m2"

F0 = PP 1
F1 = PP(OO_F0^1 ++ OO_F0^1{1})
assert isWellDefined F1; nefGenerators F1, effGenerators F1

F2 = PP(OO_F1^1 ++ OO_F1^1{0,1})
assert isWellDefined F2; nefGenerators F2, effGenerators F2

F2' = fano(3, 6)
assert isWellDefined F2'; nefGenerators F2', effGenerators F2'

sameToricVariety'(F2, F2', matrix{{1,0,0},{0,-1,0},{0,0,1}})

primaryDecomposition ideal F2
primaryDecomposition ideal F2'

---
restart
needs "ProjectiveBundles.m2"

A = PP{1,2,3}
nefGenerators A, effGenerators A
B = PP(OO_A^3)
nefGenerators B, effGenerators B
-- TODO: is A**A a projective bundle?
