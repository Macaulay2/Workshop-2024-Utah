
orlovTruncateLess = method()
orlovTruncateLess(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list submatrixByDegrees(F.dd_j, (0, i-1), (0, i-1));
    complex(mapList)
    )

orlovTruncateGeq = method()
orlovTruncateGeq(Complex, ZZ) := (F, i) -> (
    m := min F;
    M := max F;
    mapList := for j from m+1 to M list submatrixByDegrees(F.dd_j, (i, 10000), (i, 10000));
    complex(mapList)
    )

end;

restart
needsPackage "Complexes"
load "Singularity_to_derived.m2"
R = ZZ/101[x_1..x_4, Degrees => {1, 1, 3, 4}]
F = koszulComplex vars R
i = 3
orlovTruncateLess(F, i)
orlovTruncateGeq(F, i)
F
