needsPackage "NormalToricVarieties"

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

ampleDegree = X -> (
    nef := nefGenerators X;
    L := apply(entries reducedRowEchelonForm(nef ** QQ),
	row -> position(row, not zero));
    entries sum(L, i -> nef_i))
ampleDivisor = X -> toricDivisor(ampleDegree X, ring X)

globalExt = (m, F, G) -> (
    -- computing global Ext^m(M, N(v))
    (M, N) := (module F, module G);
    X := variety F;
    d := dim X;
    u := ampleDegree X;
    v := 0 * u;
    nef := coneFromVData nefGenerators X;
    -- find e that satisfies Greg's conditions
    e := binarySearch(0, , e -> (
	    -- TODO: the paper asks for S_{e*u} M, is truncation the same?
	    C := freeResolution(truncate(e * u, M), LengthLimit => m);
	    D := freeResolution(N, LengthLimit => d-m);
	    all((0,0) .. (m,d-m),
		(k,i) -> all(unique degrees C_(m-k) ** unique degrees D_i, -- TODO: why m-k and not just k?
		    (aM, aN) -> contains(nef, transpose matrix {v + aM - aN})))
	    ));
    -- TODO: recover the Yoneda sheaf extension
    E := part(v, Ext^m(truncate(e * u, M), N));
    E.cache.Ext = (m, F, G);
    E)

ExtTable = (X, L) -> (
    T = (degreesRing 1)_0;
    matrix table(L, L,
	(F, G) -> sum(dim X + 1,
	    m -> T^m * rank globalExt(m, F, G))))

------------
end
restart
needs "GlobalExt.m2"

n = 2
X = toricProjectiveSpace n
S = ring X

assert(globalExt(n, OO_X^{n+1}, OO_X^{0}) === QQ^1)

-- Beilinson's collection of O's
L = apply(n+1, i -> OO_X^{i})
ExtTable(X, L)

-- Beilinson's collection of Omega's
-- uuhhhh did nobody notice that this list is backwards??
L = apply(n+1, i -> cotangentSheaf(i, X) ** OO_X^{i})
ExtTable(X, L)

------------
end
restart
needs "GlobalExt.m2"

n = 2
X = toricProjectiveSpace n
S = ring X

K = koszulComplex vars S
F = cotangentSheaf X

L  = apply(n+1, i -> (naiveTruncation(K, (i+1, n+1)))[i+1])
L' = apply(n+1, i -> (naiveTruncation(K, (0,     i)))[i])
assert all(n+1, i -> sheaf HH_0 L#i == sheaf HH_0 L'#i)
-- FIXME: why does Omega^2 have coefficients?
apply(n+1, i -> sheaf HH_0 L#i == exteriorPower(i, F))

-- TODO: show that L (or L') form exceptional collections (up to appropriate shifts)

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
ExtTable(X, L)

-- TODO: now try their resolutions
apply(L, M -> sheaf freeResolution module M)

