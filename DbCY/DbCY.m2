needsPackage "Depth"

path = prepend("../packages", path)
loadPackage("Truncations", Reload => true)
loadPackage("Complexes", Reload => true)
loadPackage("Varieties", Reload => true)

orlovTruncateLess = method()
-- Input: a complex F of graded free modules and an integer i.
-- Output: the subcomplex of F given by summands of the form R(j) with j > -i (so R(j) is generated in degree < i)
orlovTruncateLess(Complex,    ZZ) := (F,   i) -> complex applyValues(F.dd.map, f -> submatrixByDegrees(f, (, i-1), (, i-1)))
-- Input: a map psi of complexes and an integer i.
-- Output: the induced map on subcomplexes as above.
orlovTruncateLess(ComplexMap, ZZ) := (psi, i) -> map(
    orlovTruncateLess(target psi, i), -- target
    orlovTruncateLess(source psi, i), -- source
    applyValues(psi.map, f -> submatrixByDegrees(f, (, i-1), (, i-1))))

orlovTruncateGeq = method()
-- Input: a complex F of graded free modules and an integer i.
-- Output: the quotient of F given by summands of the form R(j) with j <= -i (so R(j) is generated in degree >= i).
orlovTruncateGeq(Complex,    ZZ) := (F,   i) -> complex applyValues(F.dd.map, f -> submatrixByDegrees(f, (i, ), (i, )))
-- Input: a map psi of complexes and an integer i.
-- Output: the induced map on quotient complexes as above.
orlovTruncateGeq(ComplexMap, ZZ) := (psi, i) -> map(
    orlovTruncateGeq(target psi, i), -- target
    orlovTruncateGeq(source psi, i), -- source
    applyValues(psi.map, f -> submatrixByDegrees(f, (i, ), (i, ))))


supTruncate = method();
-- Input: a finitely generated module M over a graded Gorenstein ring with nonnegative Gorenstein parameter, and an integer i.
--        We recall that the Gorenstein parameter is the integer a such that Ext^d_R(k, R) = k(-a) (up to a homological
--        shift), where R is the ring of M, d is the dimension of R, and k is the residue field of R. 
-- Output: an integer, call it N, satisfying the following: if F is the minimal free resolution of M, then the
--         homology of the dual of orlovTruncateGeq(F, i) is concentrated in homological degrees -N, ..., 0.
supTruncate(Module, ZZ) := (M, i) -> (
    R := ring M;
    d := dim R;
    t := min flatten degrees M;-- this is the minimum generating degree of M
    if i >= t then d+i-t else d
)
-- Input: a Complex C given by a finitely generated module concentrated in a single homological degree, and an integer i. The ring
--	  of the module should be as in the input of supTruncate(Module, ZZ). 
-- Output: an integer, call it N, satisfying the following: if F is the minimal free resolution of C, and
--	   C is concentrated in degree m, then the homology of the dual of orlovTruncateGeq(F, i) is concentrated in homological
--	   degrees -N, ..., -m.
supTruncate(Complex, ZZ) := (C, i) -> (
    R := ring C;
    m := min C;
    d := dim R;
    t := min flatten degrees C_m;-- this is the minimum generating degree of C_m
    if i >= t then d+i-t + m else d + m
)

orlovTruncateGeqDualize = method();
-- Input: a graded module M and an integer i
-- Output: a smart truncation of the dual of orlovTruncationGeq(F, i) that is quasi-isomorphic to
--         the complex orlovTruncationGeq(F, i), where F is the (typically infinite) minimal free resolution of M.
orlovTruncateGeqDualize(Module, ZZ) := (M, i) -> (
    F := freeResolution(M, LengthLimit => supTruncate(M, i) + 2);
    Fi := orlovTruncateGeq(F, i);
    Fidual := dual Fi;
    canonicalTruncation(Fidual, -supTruncate(M, i) -1, )
)
orlovTruncateGeqDualize(Complex, ZZ) := (C, i) -> (
    F := resolution(C, LengthLimit => supTruncate(C, i) + 2);
    Fi := orlovTruncateGeq(F, i);
    Fidual := dual Fi;
    canonicalTruncation(Fidual, -supTruncate(C, i) -1, )
)
-- THIS FUNCTION DOESN'T WORK YET! We need the canonicalTruncation function for maps of complexes. See comment in code.
-- Input: a morphism f of graded modules and an integer i.
-- Output: the induced map on truncateGeqDualize applied to the source and target of f (and i).
orlovTruncateGeqDualize(Matrix, ZZ) := (f, i) -> (
    M := source f;
    N := target f;
    s := max{supTruncate(M) + 2, supTruncate(N) + 2};
    g := freeResolution(f, LengthLimit => s);
    gi := orlovTruncateGeq(ftilde, i);
    gidual := dual gi;
    canonicalTruncation(gidual, -s - 1)--this function doesn't exist for ComplexMaps yet.
)

sup = method();
sup(Complex) := (C) -> (
    for i from -max C to -min C -1 do (
	if prune HH_(-i)(C) != 0 then return -i
	);
    )

singularityToDerived = method();
--Input: a finitely generated module M over a graded Gorenstein ring with nonnegative Gorenstein parameter, and integers i and j.
--       We recall that the Gorenstein parameter is the integer a such that Ext^d_R(k, R) = k(-a) (up to a homological
--       shift), where R is the ring of M, d is the dimension of R, and k is the residue field of R.
--Output: Let D^{sing}(R) denote the singularity category of R, i.e. the quotient of the bounded derived category
--	  of graded R-modules by the subcategory perfect complexes. As in Orlov's paper "Derived categories
--	  of coherent sheaves and triangulated categories of singularities", we denote by \Phi_i the fully
--	  faithful functor D^{sing}(R) --> D^b(Proj(R)) constructed in that paper (see Theorem 2.5).
--	  View M as an object in D^b(R) concentrated in homological degree 0, and hence also an object in D^{sing}(R).
--	  This method outputs \Phi_i(M) (thought of as a complex of graded modules, rather than sheaves). This complex
--        is unbounded in negative homological degrees; we brutally truncate its tail so that it has length j.
--CAVEAT: Any object in D^{sing}(R) is isomorphic to a (maximal Cohen-Macaulay) module, but concentrated in some
--	  (possibly nonzero) homological degree. This method assumes the module is concentrated in homological degree
--	  zero. Should allow for more generality.

singularityToDerived(Module, ZZ, ZZ) := (M, i, j) -> (
    R := ring M;
    d := dim R;
    kk := coker vars R;
    if (flatten degrees prune Ext^d(kk, R^1))_0 > 0 then error "The Gorenstein parameter is negative.";
    G := resolution(orlovTruncateGeqDualize(M, i), LengthLimit => j);
    Gdual := dual G;
    orlovTruncateLess(Gdual, i)
)

--Input: a bounded Complex C of finitely generated modules as in the above function, and i and j as in the above function.
--Output: \Phi_i(C), where \Phi_i(C) is as described in the previous function, with tail brutally truncated so that
--	  it has length j.
--WARNING: This function needs more testing. It may not work.
--CAVEAT: this function is not functorial. Modeling maps between objects in the singularity category
--	  is only feasible when the objects are MCM modules concentrated in the same homological degree.
singularityToDerived(Complex, ZZ, ZZ) := (C, i, j) -> (
    R := ring C;
    d := dim R;
    kk := coker vars R;
    if (flatten degrees prune Ext^d(kk, R^1))_0 > 0 then error "The Gorenstein parameter is negative.";
    s := sup C;
    F := (resolution C, LengthLimit => s + 1 - max{0,min C});
    trunc := naiveTruncation(F, s, );
    N := HH_(max C) trunc;
    M := complex({map(N, R^0, 0)}, Base => max C);
    G := resolution(orlovTruncateGeqDualize(M, i), LengthLimit => j);
    Gdual := dual G;
    orlovTruncateLess(Gdual, i)
)

-- THIS FUNCTION DOESN'T WORK YET! We need two things:
--    (1) need truncateGeqDualize to work for a matrix.
--    (2) given a map of complexes, need to be able to compute the induced map on minimal free resolutions of the complexes.
-- Input: a morphism f of graded maximal Cohen-Macaulay modules and integers i and j.
-- Output: the induced map on singularityToDerived applied to the source and target of f (and also i and j). Note:
--         the space of morphisms between MCM modules in the singularity category is given by "stable" R-linear maps;
--         see Proposition 1.11 in Orlov's paper  "Derived categories of coherent sheaves and triangulated categories
--	   of singularities". In particular, every morphism between these objects in the singularity category can
--	   be represented by an honest map of modules.
singularityToDerived(Matrix, ZZ, ZZ) := (f, i, j) -> (
    R := ring f;
    d := dim R;
    kk := coker vars R;
    if (flatten degrees prune Ext^d(kk, R^1))_0 > 0 then error "The Gorenstein parameter is negative.";
    if depth(target f) < d or depth(source f) < d then error "Not a map of MCM modules.";
    g := resolution(orlovTruncateGeqDualize(f, i), LengthLimit => j);
    gdual := dual g;
    orlovTruncateLess(gdual, i)
)

end;


restart
load "DbCY.m2"
R = ZZ/101[x_0] / ideal(x_0^3)
M = coker vars R
assert try ( singularityToDerived(M, 1, 1); false ) else true

restart
load "DbCY.m2"
R = ZZ/101[x_0..x_4] / ideal(x_0*x_1, x_2*x_3*x_4)
X = Proj R
M = coker matrix{{x_0*x_2}}
D = minimize singularityToDerived(M, 3, 7)
D = sheaf D
errorDepth=2
debugLevel=2
RHom(OO_X^1, D)

prune HH D
D.dd
elapsedTime prune HH RHom(D, OO_X^1) -- 22s -> 16s, ranks 113 103 69 7
elapsedTime RHom(OO_X^1, D)
debug Varieties


restart
load "DbCY.m2"
S = ZZ/101[x_0..x_4]
f = sum for i from 0 to 4 list x_i^5
R = S / ideal(f)
M = coker vars R
i = 0;
j = 7;
F = singularityToDerived(M, i, j)
dual((res M)[4]) -- this is RHom(k, R)[4], which is isomorphic to k(a) = k by the Gorenstein property
--This implies that \widetilde{F} = O_X[-3].
G = res M
K = coker G.dd_5
singularityToDerived(K, 0, 9)
prune ker oo.dd_(-1)
F = singularityToDerived(M, i, 6)
G = singularityToDerived(K, 0, 10)
G[-4] == F
F.dd_3
(G[-4]).dd_3
for i from -1 to 3 do print HH_i(o8) == 0

prune HH_0(o8) == 0
prune HH_(-1)(o8) == 0
prune HH_1(o8) == 0
prune HH_2(o8) == 0
prune HH_3(o8) == 0


--DEMO
restart
load "DbCY.m2"
S = ZZ/101[x_0..x_4]
f = sum for i from 0 to 4 list x_i^5
R = S / ideal(f) --affine cone of the quintic
kk = coker vars R --the residue field of R
F = singularityToDerived(kk, 0, 7)
--This computes \Phi(kk) (up to sheafification).
--We will eventually interface with "Varieties" package to implement the associated complex of sheaves.
--The complex \Phi(kk) has infinite length; it is unbounded on the right. We take a length 7 approximation here.
--The tail is a matrix factorization.

--What is the homology of this complex?
homologies = for i from 0 to 3 list prune HH_i F
--Let's sheafify these:
sheaf HH_0 F == 0
sheaf prune HH_3 F
--Thus, \Phi_0(kk) = O_X[-3].

--Future work: implement the functor going the opposite direction.

