debug loadPackage("Truncations", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Truncations.m2", Reload => true);
debug loadPackage("Complexes", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Complexes.m2", Reload => true);
debug installPackage("Varieties", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Varieties.m2");

koszulComplex(Ideal) := I -> (
    return koszulComplex(gens I);
);
koszulComplex(SheafMap) := phi -> (
    C := koszulComplex((phi).map);
    return sheaf C;
);
