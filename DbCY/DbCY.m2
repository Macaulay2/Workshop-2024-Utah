needsPackage "Complexes"
needsPackage "TorAlgebra"


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

-- Input: a finitely generated module M over a graded Gorenstein ring with nonnegative Gorenstein parameter, and an integer i.
--        We recall that the Gorenstein parameter is the integer a such that Ext^d_R(k, R) = k(-a) (up to a homological
--        shift), where R is the ring of M, d is the dimension of R, and k is the residue field of R. 
-- Output: an integer, call it N, satisfying the following: if F is the minimal free resolution of M, then the
--         homology of the dual of orlovTruncateGeq(F, i) is concentrated in homological degrees -N, ..., 0.
supTruncate = method();
supTruncate(Module, ZZ) := (M, i) -> (
    R := ring M;
    d := dim R;
    t := min flatten degrees M;
    -- get the min degree of all gens of the module
    if i >= t then d+i-t else d
)


truncateGeqDualize = method();
truncateGeqDualize(Module, ZZ) := (M, i) -> (
    F := freeResolution(M, LengthLimit => supTruncate(M, i));
    Fi := orlovTruncateGeq(F, i);
    Fidual := dual Fi;
    canonicalTruncation(Fidual, -supTruncate(M, i) + 1, )
)

singularityToSheaves = method();
singularityToSheaves(Module, ZZ, ZZ) := (M, i, j) -> (
    R := ring M;
    d := dim R;
    kk := coker vars R;
    if (flatten degrees prune Ext^d(kk, R^1))_0 > 0 then error "The Gorenstein parameter is negative.";
    G := resolution(truncateGeqDualize(M, i), LengthLimit => j);
    D := dual G;
    orlovTruncateLess(D, i)
)
end;

restart
load "DbCY.m2"
R = ZZ/101[x_0] / ideal(x_0^3)
M = coker vars R
singularityToSheaves(M, 1, 1)
restart
load "DbCY.m2"
R = ZZ/101[x_0..x_4] / ideal(x_0*x_1, x_2*x_3*x_4)
d = dim(R)
kk = coker vars R
degrees Ext^d(kk, R^1)

M = coker matrix{{x_0*x_2}}
singularityToSheaves(M, 3, 7)
F = freeResolution (M, LengthLimit => 3)
i = 3
f = id_F
orlovTruncateGeq(id_F, i)
dual id_F
flatten degrees Ext^d(kk, R^1)
