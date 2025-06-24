needsPackage "NormalToricVarieties"

ToricMap^* := f -> pullback_f

-- TODO: move to Core, c.f. https://github.com/Macaulay2/M2/issues/3844
-- assume 'test' is a monotonic function, i.e. false for all i < n then true for i >= n
binarySearch = method()
-- return the first index of an element in L such that test(L#i) is true
binarySearch(List,        Function) := (L,          test) -> binarySearch(0, #L-1, i -> test(L#i))
-- shorthand for search in [0, n)
binarySearch(         ZZ, Function) := (      high, test) -> binarySearch(0, high-1, test)
-- shorthand for when a lower bound isn't known
binarySearch(Nothing, ZZ, Function) := (null, high, test) -> (
    dist := 1;
    while true do if test(high - dist)
    then (high, dist) = (high - dist, dist * 2)
    else break binarySearch(high - dist + 1, high, test));
-- shorthand for when an upper bound isn't known
binarySearch(ZZ, Nothing, Function) := (low, null, test) -> (
    dist := 1;
    while true do if test(low + dist)
    then break binarySearch(low, low + dist, test)
    else (low, dist) = (low + dist + 1, dist * 2))
-- standard binary search
binarySearch(ZZ, ZZ, Function) := (low, high, test) -> (
    -- TODO: are the first two lines standard?
    if     test(low)  then return low;
    --if not test(high) then return high + 1;
    while high - low > 1 do (
	mid := (high + low) // 2;
	if test(mid) then high = mid else low = mid);
    high)

debug needsPackage "Truncations"
toricDivisor(List, Ring) := opts -> (d, S) -> (
    if not instance(X := variety S, NormalToricVariety)
    then error "expected a degree and the Cox ring of a normal toric variety";
    P := basisPolyhedron(effGenerators S, transpose matrix {d});
    toricDivisor(first entries transpose interiorPoint P, X, opts))

ampleDegree = X -> X.cache#"AmpleDegree" ??= (
    nef := nefGenerators X;
    L := apply(entries reducedRowEchelonForm(nef ** QQ),
	row -> position(row, not zero));
    entries sum(L, i -> nef_i))
ampleDivisor = X -> toricDivisor(ampleDegree X, ring X)

globalExt = method()
globalExt(ZZ, CoherentSheaf, CoherentSheaf) := Module => (m, F, G) -> F.cache#(symbol globalExt, m, F, G) ??= (
    -- computing global Ext^m(M, N(v))
    (M, N) := (module F, module G);
    X := variety F;
    d := dim X;
    u := ampleDegree X;
    v := 0 * u;
    nef := coneFromVData nefGenerators X;
    -- find e that satisfies Greg's conditions
    C := freeResolution(M, LengthLimit => m);
    D := freeResolution(N, LengthLimit => d-m);
    e := binarySearch(sum min degrees sum C, sum max degrees sum D, e -> (
	    -- TODO: the paper asks for S_{e*u} M, is truncation the same?
	    C = freeResolution(truncate(e * u, M, MinimalGenerators => false), LengthLimit => m);
	    -- TODO: in the single graded case we can just take the maximum degree
	    -- but to make this work in the multigraded case, that may not work!
	    all((0,0) .. (m,d-m),
		(k,i) -> all(unique degrees C_(m-k) ** unique degrees D_i, -- TODO: why m-k and not just k?
		    (aM, aN) -> contains(nef, transpose matrix {v + aM - aN})))
	    ));
    if debugLevel > 0 then printerr("using truncation limit ", toString(e * u));
    -- TODO: recover the Yoneda sheaf extension
    E := part(v, Ext^m(truncate(e * u, M, MinimalGenerators => false), N, MinimalGenerators => false));
    E.cache.Ext = (m, F, G);
    E)

globalHom = (C, D) -> C.cache#("globalHom", C, D) ??= (
    z := degree 1_(ring C);
    minimize part_z Hom(C, D, DegreeLimit => z,
	MinimalGenerators => false))

globalExt(ZZ, Complex, Complex) := Module => (m, F, G) -> F.cache#(symbol globalExt, m, F, G, 1) ??= (
    -- TODO: also need to push forward to ambient projective space
    (C, D) := (F, G); -- (module F, module G);
    X := variety ring F;
    d := dim X; -- should be embedding dimension
    u := ampleDegree X;
    v := 0 * u;
    nef := coneFromVData nefGenerators X;
    s := first concentration D;
    if 0 != s then (C, D) = (C[s], D[s]);
    -- find r that satisfies inequality in Theorem 2.14
    if #u == 1 then (
	r := max for j to last concentration D
	list max for i to pdim D_j -- TODO: what are n and l in the paper?
	-- TODO: translate to a containment of cones for the toric case
	list max apply(keys betti freeResolution D_j, (k, aa, s) -> aa) - d * u);
    -- just for fun, we compute the bound a different way and compare
    C' := freeResolution(C, LengthLimit => m);
    D' := freeResolution(D, LengthLimit => d-m);
    e := binarySearch(sum min degrees sum C, sum max degrees sum D', e -> (
	    -- TODO: the paper asks for S_{e*u} M, is truncation the same?
	    C' = freeResolution(truncate(e * u, C, MinimalGenerators => false), LengthLimit => m);
	    -- TODO: in the single graded case we can just take the maximum degree
	    -- but to make this work in the multigraded case, that may not work!
	    all((0,0) .. (m,d-m),
		(k,i) -> all(unique degrees C'_(m-k) ** unique degrees D'_i, -- TODO: why m-k and not just k?
		    (aC, aD) -> contains(nef, transpose matrix {v + aC - aD})))
	    ));
    if debugLevel > 0 then printerr("using truncation limit ", toString r, " vs ", toString(e * u));
    -- use e for the multigraded case
    --if #u > 1 then
    r = e * u;
    C' = freeResolution(truncate(r, C, MinimalGenerators => false),
	-- TODO: this +2 seems extra, but some examples fail without it
	LengthLimit => m - min(0, first concentration C) + 2);
    E := globalHom(C', D); -- ~70% of the computation
    E.cache.Ext = (F, G);
    -- TODO: why is it -m here?!
    E_(-m))

ExtTable = (X, L) -> (
    T := (degreesRing 1)_0;
    matrix table(L, L,
	(F, G) -> sum(dim X + 1,
	    m -> T^m * rank globalExt(m, F, G))))

------------
end
restart
needs "GlobalExt.m2"
debugLevel=1

n = 2
X = toricProjectiveSpace n
S = ring X

assert(globalExt(n, OO_X^{n+1}, OO_X^{0}) === QQ^1)

-- Beilinson's collection of O's
L = apply(n+1, i -> OO_X^{i})
ExtTable(X, L) -- upper triangular means exceptional
elapsedTime assert(0 == ExtTable(X, L) - ExtTable(X, complex \ module \ L))

-- Beilinson's collection of Omega's
L = apply(n+1, i -> cotangentSheaf(i, X) ** OO_X^{i})
ExtTable(X, L) -- upper triangular means exceptional
-- uuhhhh did nobody notice that this list is backwards??
elapsedTime assert(0 == ExtTable(X, L) - ExtTable(X, freeResolution \ module \ L))

N = module cotangentSheaf(1, X) ** S^{1}
M = module cotangentSheaf(2, X) ** S^{2}
0 == Ext^1(sheaf N, sheaf M)
0 == globalExt(1, sheaf N, sheaf M)
0 == globalExt(1, freeResolution N, freeResolution M)

--
K = koszulComplex vars S
L1 = apply(n+1, i -> (naiveTruncation(K, (i+1, n+1)))[i+1] ** S^{i})
L2 = apply(n+1, i -> (naiveTruncation(K, (0,     i)))[i]   ** S^{i})
assert all(n+1, i -> sheaf HH_0 L1#i == L#i and sheaf HH_0 L2#i == L#i)

-- show that L1 (or L2) form exceptional collections (up to appropriate shifts)
ExtTable(X, L) -- upper triangular means exceptional
assert(0 == ExtTable(X, L) - ExtTable(X, L1))
assert(0 == ExtTable(X, L) - ExtTable(X, L2))

------------
end
restart
needs "GlobalExt.m2"

Y = toricProjectiveSpace 2
X = toricBlowup({0,1}, Y)
S = ring X
p = X^[]

-- Orlov's formula for the blow up
L = { p^* OO_Y^{0}, p^* OO_Y^{1}, p^* OO_Y^{2}, sheaf(S^1/S_3) }
ExtTable(X, L) -- upper triangular means exceptional
-- now compare with complexes version
assert(0 == ExtTable(X, L) - ExtTable(X, apply(L, M -> freeResolution module M)))



--
(L1,L2) = (OO_X^1, sheaf(S^1/S_2))
--ExtTable(X, {L1,L2})
(C,D) = freeResolution \ module \ (L1,L2)
(F,G) = sheaf \ (C,D)

globalExt(0,F,G)
globalExt(1,G,F)

ExtTable(X, {F,G})
ExtTable(X, {L1,L2})

part_{0,0} Ext^1(module L2, module L2)
part_{0,0} Hom(truncate_{3,3} D, D[1])
globalExt(1, L2, L2)



globalExt(0,C,D), Ext^0(L1, L2)
globalExt(1,D,C), Ext^1(L2, L1)

minimize part_{0,0} Hom(D, C)



end--
restart
needs "GlobalExt.m2"

-- FIXME
X = toricProjectiveSpace 2
S = ring X
F = cotangentSheaf(2, X)
C = freeResolution module cotangentSheaf(1, X) ** S^{1}
minimize part_0 Hom(truncate_1 C, freeResolution module F ** S^{2})
minimize part_0 Hom(truncate_1 C, freeResolution module prune F ** S^{2})
minimize part_0 Hom(truncate_1 C, freeResolution module(F ** OO_X^{2}))
minimize part_0 Hom(truncate_1 C, module((sheaf freeResolution module F) ** OO_X^{2}))


minimize part_0 Hom(truncate_1 C, module(sheaf module freeResolution module M ** OO_X^{2}))

restart
needs "GlobalExt.m2"
n = 2
X = toricProjectiveSpace 2
S = ring X

N = module cotangentSheaf(1, X) ** S^{1}
M = module cotangentSheaf(2, X) ** S^{2}
Ext^1(sheaf N, sheaf M)
globalExt(1, sheaf N, sheaf M)
debugLevel = 1
globalExt(1, freeResolution N, freeResolution M)
globalExt(1, freeResolution N, freeResolution module prune sheaf M)


Ext^2(OO_X^{n+2}, OO_X^{0})
globalExt(2, OO_X^{n+2}, OO_X^{0})
globalExt(2, complex module OO_X^{n+2}, complex module OO_X^{0})
globalExt(2, complex S^{n+2}, complex prune truncate(1, S^1))
globalExt(2, complex S^{n+2}, complex module prune sheaf truncate(1, S^1))
