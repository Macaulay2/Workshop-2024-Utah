-*
TODO:
- add documentation (Souvik, Geoff, Tim)
- add tests from the paper (Guanyu, Geoff)
- clean up naming and input conventions (Michael)
- RHom and liftAlongQuasiIsomorphism changes (Michael, Mahrud)
- turn examples from section 2 to code in section 4 (Guanyu, Tim, Geoff)
*-

newPackage(
    "OrlovFunctors",
    Headline => "Orlov's functors between singularity and derived categories for Calabi-Yau varieties",
    Version => "0.1",
    Date => "June 3, 2025",
    Authors => {
	{ Name => "Michael K. Brown", Email => "mkb0096@auburn.edu", HomePage => "http://webhome.auburn.edu/~mkb0096/" },
    { Name => "Souvik Dey", Email => "souvikd@uark.edu", HomePage => "https://sites.google.com/view/souvikdey"}, 
	{ Name => "Mahrud Sayrafi", Email => "mahrud@mcmaster.ca", HomePage => "https://mahrud.github.io" },
	{ Name => "Guanyu Li", Email => "gl479@cornell.edu", HomePage => "https://sites.google.com/view/guanyu-li-math/home" },
	{ Name => "Geoffrey Fatin", Email => "glf55@cornell.edu", HomePage => "https://physics.cornell.edu/geoffrey-fatin"},
	{ Name => "Tim Tribone", Email => "tim.tribone@utah.edu", HomePage => "https://timtribone.com/"}
    },
    PackageImports => { "Depth" },
    PackageExports => { "Complexes", "Varieties" },
    AuxiliaryFiles => true,
    DebuggingMode => true,
)

export {
    "orlovTruncateGeq",
    "orlovTruncateGeqDualize",
    "orlovTruncateLess",
    "singularityToDerived",
    "supTruncate",
    }

-* Code section *-

------------------------------------------------------------------------------

orlovTruncateLess = method()
-- Input: a complex F of graded free modules and an integer i.
-- Output: the subcomplex of F given by summands of the form R(j) with j > -i (so R(j) is generated in degree < i)
orlovTruncateLess(ZZ, Complex)    := Complex    => (i, F) -> (
    complex applyValues(F.dd.map, f -> submatrixByDegrees(f, (, i-1), (, i-1))))
-- Input: a map psi of complexes and an integer i.
-- Output: the induced map on subcomplexes as above.
orlovTruncateLess(ZZ, ComplexMap) := ComplexMap => (i, psi) -> map(
    orlovTruncateLess(i, target psi), -- target
    orlovTruncateLess(i, source psi), -- source
    applyValues(psi.map, f -> submatrixByDegrees(f, (, i-1), (, i-1))))


orlovTruncateGeq = method()
-- Input: a complex F of graded free modules and an integer i.
-- Output: the quotient of F given by summands of the form R(j) with j <= -i (so R(j) is generated in degree >= i).
orlovTruncateGeq(ZZ, Complex)    := Complex    => (i, F) -> (
    complex applyValues(F.dd.map, f -> submatrixByDegrees(f, (i, ), (i, ))))
-- Input: a map psi of complexes and an integer i.
-- Output: the induced map on quotient complexes as above.
orlovTruncateGeq(ZZ, ComplexMap) := ComplexMap => (i, psi) -> map(
    orlovTruncateGeq(i, target psi), -- target
    orlovTruncateGeq(i, source psi), -- source
    applyValues(psi.map, f -> submatrixByDegrees(f, (i, ), (i, ))))

------------------------------------------------------------------------------

supTruncate = method()
-- Input: a finitely generated module M over a graded Gorenstein ring with nonnegative Gorenstein parameter, and an integer i.
--        We recall that the Gorenstein parameter is the integer a such that Ext^d_R(k, R) = k(-a) (up to a homological
--        shift), where R is the ring of M, d is the dimension of R, and k is the residue field of R. 
-- Output: an integer, call it N, satisfying the following: if F is the minimal free resolution of M, then the
--         homology of the dual of orlovTruncateGeq(F, i) is concentrated in homological degrees -N, ..., 0.
supTruncate(ZZ, Module) := (i, M) -> (
    R := ring M;
    d := dim R;
    t := min flatten degrees M;-- this is the minimum generating degree of M
    if i >= t then d+i-t else d)
-- Input: a Complex C given by a finitely generated module concentrated in a single homological degree, and an integer i. The ring
--	  of the module should be as in the input of supTruncate(Module, ZZ). 
-- Output: an integer, call it N, satisfying the following: if F is the minimal free resolution of C, and
--	   C is concentrated in degree m, then the homology of the dual of orlovTruncateGeq(F, i) is concentrated in homological
--	   degrees -N, ..., -m.
supTruncate(ZZ, Complex) := (i, C) -> (
    R := ring C;
    m := min C;
    d := dim R;
    t := min flatten degrees C_m;-- this is the minimum generating degree of C_m
    if i >= t then d+i-t + m else d + m)

------------------------------------------------------------------------------

orlovTruncateGeqDualize = method()
-- Input: a graded module M and an integer i
-- Output: a smart truncation of the dual of orlovTruncationGeq(F, i) that is quasi-isomorphic to
--         the complex orlovTruncationGeq(F, i), where F is the (typically infinite) minimal free resolution of M.
orlovTruncateGeqDualize(ZZ, Module) := (i, M) -> (
    F := freeResolution(M, LengthLimit => supTruncate(i, M) + 2);
    Fi := orlovTruncateGeq(i, F);
    Fidual := dual Fi;
    canonicalTruncation(Fidual, -supTruncate(i, M) -1,  ))
-- Input: a Complex C given by a finitely generated module concentrated in a single homological degree, and an integer i
-- Output: a smart truncation of the dual of orlovTruncationGeq(F, i) that is quasi-isomorphic to
--         the complex orlovTruncationGeq(F, i), where F is the (typically infinite) minimal free resolution of M.
orlovTruncateGeqDualize(ZZ, Complex) := (i, C) -> (
    F := freeResolution(C, LengthLimit => supTruncate(i, C) + 2);
    Fi := orlovTruncateGeq(i, F);
    Fidual := dual Fi;
    canonicalTruncation(Fidual, -supTruncate(i, C) -1,  ))
-- THIS FUNCTION DOESN'T WORK YET! We need the canonicalTruncation function for maps of complexes. See comment in code.
-- Input: a morphism f of graded modules and an integer i.
-- Output: the induced map on truncateGeqDualize applied to the source and target of f (and i).
orlovTruncateGeqDualize(ZZ, Matrix) := (i, f) -> (
    M := source f;
    N := target f;
    s := max{supTruncate(i, M) + 2, supTruncate(i, N) + 2};
    g := freeResolution(f, LengthLimit => s);
    gi := orlovTruncateGeq(i, g);
    gidual := dual gi;
    canonicalTruncation(gidual, -s - 1, ))

-- f = map(N, M, 1)
-- assert isHomogeneous f
-- F = res(f, LengthLimit => 3)
-- F1 = orlovTruncateGeq(1, F)
-- F1' = dual F1
-- TODO: orlovTruncateGeqDualize(i, f) needs to match phi:
-- phi = canonicalTruncation(F1', -2, 0)

-- TODO: this function doesn't exist for ComplexMaps yet.

------------------------------------------------------------------------------

isMCM = M -> depth M >= dim ring M

sup = method()
sup(Complex) := (C) -> (
    for i from -max C to -min C -1 do (
	if prune HH_(-i)(C) != 0 then return -i
	);
    )

singularityToDerived = method(Options => { LengthLimit => null })
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
singularityToDerived(ZZ, Module) := Complex => opts -> (i, M) -> (
    -- TODO: check if M is MCM?
    -- TODO: what does LengthLimit for a resolution of a complex do?
    G := freeResolution(orlovTruncateGeqDualize(i, M) -*, opts *-);
    orlovTruncateLess(i, dual G))

--Input: a bounded Complex C of finitely generated modules as in the above function, and i and j as in the above function.
--Output: \Phi_i(C), where \Phi_i(C) is as described in the previous function, with tail brutally truncated so that
--	  it has length j.
--WARNING: This function needs more testing. It may not work.
--CAVEAT: this function is not functorial. Modeling maps between objects in the singularity category
--	  is only feasible when the objects are MCM modules concentrated in the same homological degree.
singularityToDerived(ZZ, Complex) := Complex => opts -> (i, C) -> (
    -- TODO: check if C is MCM?
    s := sup C;
    F := freeResolution(C, LengthLimit => s + 1 - max{0,min C});
    M := complex(HH_(max C) naiveTruncation(F, s, ), Base => max C);
    G := freeResolution(orlovTruncateGeqDualize(i, M), opts);
    orlovTruncateLess(i, dual G)
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
singularityToDerived(ZZ, Matrix) := ComplexMap => opts -> (i, f) -> (
    if not isMCM source f or not isMCM target f
    then error "expected Maximally Cohen-Macaulay source and target";
    -- g := freeResolution(orlovTruncateGeqDualize(i, f), opts);
    -- orlovTruncateLess(i, dual g)
    F := freeResolution(f, opts);
    F1' := dual orlovTruncateGeq(i, F);
    Gs := singularityToDerived(i, source f);
    Gt := singularityToDerived(i, target f);
    -- FIXME: figure out the bounds here
    phi := canonicalTruncation(F1', -2, 0);
    -- TODO: we need the target to be === to G0t
    -- but canonicalTruncation isn't ===-functorial!
    -- c.f. https://github.com/Macaulay2/M2/issues/3865
    G0s := canonicalTruncation(source F1', -2, 0);
    G0t := canonicalTruncation(target F1', -2, 0);
    phi = map(G0t, source phi, phi);
    f = phi * resolutionMap G0s;
    g := resolutionMap G0t;
    assert(target f === target g);
    h := liftMapAlongQuasiIsomorphism(f, g); -- or f // g
    -- homotopyMap h -- is this useful for anything?
    orlovTruncateLess(i, dual h))

-- TODO: this requires functorial RHom first
-- singularityToDerived(ZZ, ComplexMap) := ComplexMap => opts -> (i, f) -> ()

------------------------------------------------------------------------------

-* Documentation section *-
beginDocumentation()

load "./OrlovFunctors/docs.m2"

-* Test section *-
load "./OrlovFunctors/tests.m2"

end--

-* Development section *-
restart
needsPackage "OrlovFunctors"
check "OrlovFunctors"
