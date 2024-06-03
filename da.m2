debug loadPackage("Truncations", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Truncations.m2", Reload => true);
debug loadPackage("Complexes", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Complexes.m2", Reload => true);
debug installPackage("Varieties", FileName => "/Users/daniel/Documents/GitHub/Workshop-2024-Utah/packages/Varieties.m2");


koszulComplex(Ideal) := Complex => {Concentration => null} >> opts -> I -> (
    koszulComplex(gens I)
);
koszulComplex(SheafMap) := Complex => {Concentration => null} >> opts -> phi -> (
    sheaf koszulComplex(matrix phi)
);
koszulComplex(ProjectiveVariety) := Complex => {Concentration => null} >> opts -> X -> (
    sheaf koszulComplex(vars(ring X))
);
koszulComplex(AffineVariety) := Complex => X -> {Concentration => null} >> opts -> (
    sheaf koszulComplex(vars(ring X))
);
