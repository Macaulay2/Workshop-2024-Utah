needsPackage "Complexes"

orlovTruncateLess = method()
orlovTruncateLess(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list {submatrixByDegrees(F.dd_j, (-10000, i-1), (-10000, i-1)), j};
    --for i from 0 to #mapList - 1 do (
	--if mapList_i = 0 then j = j +1; 
	--);
    complex(mapList)[m]
    )

orlovTruncateGeq = method()
orlovTruncateGeq(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list submatrixByDegrees(F.dd_j, (i, 10000), (i, 10000));
    complex(mapList)[m]
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
i = 3
j = supTruncate(M, i) + dim(R)
singularityToSheaves(M, i, supTruncate(M, i) + dim(R))

F
