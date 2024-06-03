debug loadPackage("Truncations", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Truncations.m2", Reload => true);
debug loadPackage("Complexes", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Complexes.m2", Reload => true);
debug installPackage("Varieties", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Varieties.m2");


R = QQ[x,y];

Sleuler(Complex) := C -> (
    d = length(C);
    c = 0;
    for i from -d to d do (
        c = c + (-1)^i * (euler(HH_i()))
    );
    return c
)

koszulComplex(Ideal) := I -> (
    return koszulComplex(gens I);
);
koszulComplex(Variety, Variety) := (X, Y) -> (
    return sheaf(koszulComplex(ideal Y));
);
koszulComplex(SheafMap) := phi -> (
    C := koszulComplex((phi)#map);
    return sheaf C;
);
