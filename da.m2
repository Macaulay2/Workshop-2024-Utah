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

-- ComplexMap * ZZ, RingElement works once we add this to SheafMaps.m2

SheafMap * ZZ := SheafMap * RingElement := (f, r) -> r * f

-- SheafMap * RingElement is odd, similarly for complexes.

X = Proj QQ[x,y,z]
F = OO_X
id_F * x_1

-- Attempt at twisting maps

ZZ * SheafMap := RingElement * SheafMap := (r, f) -> (
    if isHomogeneous(r) == false then (
        error "expected homogeneous element"
    );
    (map(target f(first degree r), (source f), r * id_(module target f) * matrix f))
);

-- SheafMap + RingElement => SheafMap + Complex not defined


