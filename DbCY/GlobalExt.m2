needsPackage "NormalToricVarieties"
needsPackage "Complexes"
importFrom("Varieties", {"flattenModule", "flattenComplex"})
needs "binarySearch.m2"

-- my system defines this already, but it's not in 1.25.06
if not isMember((module, Complex), methods(module, Complex)) then module Complex := identity

debug needsPackage "Truncations"
toricDivisor(List, Ring) := opts -> (d, S) -> (
    if not instance(X := variety S, NormalToricVariety)
    then error "expected a degree and the Cox ring of a normal toric variety";
    P := basisPolyhedron(effGenerators S, transpose matrix {d});
    toricDivisor(first entries transpose interiorPoint P, X, opts))

RHom' = method()
RHom'(ZZ, CoherentSheaf, CoherentSheaf) := (m, F, G) -> RHom'(m, complex module F, complex module G)
RHom'(ZZ, Complex, Complex) := Module => (m, F, G) -> F.cache#(symbol RHom, m, F, G, 1) ??= (
    -- TODO: also need to push forward to ambient projective space
    (C, D) := (module F, module G);
    X := variety ring F;
    d := dim X; -- should be embedding dimension
    w := sum degrees ring F;
    z := {0};
    u := {1};
    s := first concentration D;
    if 0 != s then (C, D) = (C[s], D[s]);
    -- find r that satisfies inequality in Theorem 2.14
    if #u == 1 then (
	r := max for j to max D
	list max for i to pdim flattenModule D_j
	list max apply(keys betti freeResolution flattenModule D_j, (k, aa, s) -> aa) - w + {1});
    C' = freeResolution(truncate(r, C, MinimalGenerators => false),
	-- TODO: this +2 seems extra, but some examples fail without it
	LengthLimit => m - min(0, first concentration C) + 2);
    E := minimize prune part_z Hom(C', D, DegreeLimit => z,
	MinimalGenerators => false);
    E.cache.RHom = (F, G);
    E^m)

RHom'(CoherentSheaf, CoherentSheaf) := (F, G) -> RHom'(complex module F, complex module G)
RHom'(Complex, Complex) := Complex => (F, G) -> (
    -- TODO: need to push forward to ambient projective space
    (C, D) := (module F[min F], module G[min G]);
    --
    Y := youngest(C, D);
    if Y.cache#?(symbol RHom, C, D) then return (
	Y.cache#(symbol RHom, C, D)[min F - min G]);
    --
    R := ring F;
    K := coefficientRing R;
    X := variety R;
    d := dim X; -- should be embedding dimension
    z := {0};
    u := {1};
    -- find r that satisfies inequality in Theorem 2.14
    if #u == 1 then (
	r := max for j to max D
	list max for i to pdim flattenModule D_j
	list max apply(keys betti freeResolution(D_j, LengthLimit => d), (k, aa, s) -> aa) - d * u);
    C' := freeResolution(truncate(r, C, MinimalGenerators => false),
	LengthLimit => max(0, d + max C - min C + 1));
    H := prune homology Hom(C', D); -- ~70% of the computation
    M := for i in (min C - max D .. max C + dim X - min D) list K^(numcols basis_z H_(-i));
    E := complex(M, Base => -max H);
    --
    E.cache.RHom = (C, D);
    Y.cache#(symbol RHom, C, D) = E;
    E[min F - min G])

globalExt = RHom'

ExtTable = (X, L) -> (
    T := (degreesRing 1)_0;
    matrix table(L, L,
	(F, G) -> sum(dim X + 1,
	    m -> T^m * rank RHom'(m, F, G))))

rankPolynomial = (T, C) -> sum(pairs C.module, (i, E) -> rank E * T^(-i))
ExtTable = (X, L) -> (
    T := (degreesRing 1)_0;
    V := table(L, L, RHom');
    matrix applyTable(V, rankPolynomial_T))

------------
end
restart
notify = true
needs "./GlobalExt.m2"
debugLevel=1

n = 2
X = toricProjectiveSpace n
S = ring X

assert(RHom'(OO_X^{n+1}, OO_X^{0}) == QQ^1[3])
assert(RHom'(n, OO_X^{n+1}, OO_X^{0}) === QQ^1)

-- Beilinson's collection of O's
L = apply(n+1, i -> OO_X^{i})
ExtTable(X, L) -- upper triangular means exceptional
elapsedTime assert(0 == ExtTable(X, L) - ExtTable(X, complex \ module \ L))

-- Beilinson's collection of Omega's
L = apply(n+1, i -> prune cotangentSheaf(i, X) ** OO_X^{i})
ExtTable(X, L) -- upper triangular means exceptional
-- uuhhhh did nobody notice that this list is backwards??
elapsedTime assert(0 == ExtTable(X, L) - ExtTable(X, freeResolution \ module \ L))

N = module cotangentSheaf(1, X) ** S^{1}
M = module cotangentSheaf(2, X) ** S^{2}
0 == RHom'^1(sheaf N, sheaf M)
0 == RHom'(1, sheaf N, sheaf M)
0 == RHom'(1, freeResolution N, freeResolution M)

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

ToricMap^* := f -> pullback_f

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

RHom'(0,F,G)
RHom'(1,G,F)

ExtTable(X, {F,G})
ExtTable(X, {L1,L2})

part_{0,0} RHom'^1(module L2, module L2)
part_{0,0} Hom(truncate_{3,3} D, D[1])
RHom'(1, L2, L2)



RHom'(0,C,D), RHom'^0(L1, L2)
RHom'(1,D,C), RHom'^1(L2, L1)

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
RHom'^1(sheaf N, sheaf M)
RHom'(1, sheaf N, sheaf M)
debugLevel = 1
RHom'(1, freeResolution N, freeResolution M)
RHom'(1, freeResolution N, freeResolution module prune sheaf M)


RHom'^2(OO_X^{n+2}, OO_X^{0})
RHom'(2, OO_X^{n+2}, OO_X^{0})
RHom'(2, complex module OO_X^{n+2}, complex module OO_X^{0})
RHom'(2, complex S^{n+2}, complex prune truncate(1, S^1))
RHom'(2, complex S^{n+2}, complex module prune sheaf truncate(1, S^1))
