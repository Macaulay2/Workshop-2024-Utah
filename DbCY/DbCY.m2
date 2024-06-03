needsPackage "Complexes"

orlovTruncateLess = method()
orlovTruncateLess(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list submatrixByDegrees(F.dd_j, (, i-1), (, i-1));
    complex(mapList, Base => m)
    )

--orlovTruncateLess = method()
orlovTruncateLess(ComplexMap, ZZ) := (f, i) -> (
    m := min{min source f, min target f};
    M := max{max source f, max target f};
    mapList := for j from m to M list submatrixByDegrees(f_j, (, i-1), (, i-1));
    map(orlovTruncateLess(target f, i), orlovTruncateLess(source f, i), mapList)
    )

orlovTruncateGeq = method()
orlovTruncateGeq(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list submatrixByDegrees(F.dd_j, (i,), (i,));
    complex(mapList, Base => m)
    )

orlovTruncateGeq(ComplexMap, ZZ) := (f, i) -> (
    m := min{min source f, min target f};
    M := max{max source f, max target f};
    mapList := for j from m to M list submatrixByDegrees(f_j, (i, ), (i,));
    map(orlovTruncateGeq(target f, i), orlovTruncateGeq(source f, i), mapList)
    )

-- Input: a complex F which is already minimal resolution, a number i
-- Output: a number, which is the upper bound
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
    canonicalTruncation(Fidual, -supTruncate(M, i) + 1, 10000)
)

singularityToSheaves = method();
singularityToSheaves(Module, ZZ, ZZ) := (M, i, j) -> (
    G := resolution(truncateGeqDualize(M, i), LengthLimit => j);
    D := dual G;
    orlovTruncateLess(D, i)
)



end;

restart
load "DbCY.m2"
R = ZZ/101[x_0..x_4] / ideal(x_0*x_1, x_2*x_3*x_4)
M = coker matrix{{x_0*x_2}}
F = freeResolution (M, LengthLimit => 3)
i = 3
f = id_F
orlovTruncateGeq(id_F, i)
dual id_F
