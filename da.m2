debug loadPackage("Truncations", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Truncations.m2", Reload => true);
debug loadPackage("Complexes", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Complexes.m2", Reload => true);
debug installPackage("Varieties", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Varieties.m2");


R = QQ[x,y];

closedSSSequence := method();




Sleuler(Complex) := C -> (
    d = length(C);
    c = 0;
    for i from -d to d do (
        c = c + (-1)^i * (euler(HH_i()))
    );
    return c
)