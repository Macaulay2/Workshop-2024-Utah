ProjectiveBundle = new Type of HashTable

projectiveBundle = method()
projectiveBundle(ProjectiveVariety, CoherentSheaf) := (X, E) -> (
--checks on X and E
    if variety E =!= X then return "error: expected coherent sheaf on X";
    if rank E == 0 then return "error: expected a nonzero sheaf";
    if not isSmooth X then print "warning: X is not smooth and results may be incorrect";
--TODO: add option to bypass local freeness check and pruning
    if not isLocallyFree E then return "error: expected a locally free sheaf";
    new ProjectiveBundle from{symbol Variety => X, symbol CoherentSheaf => E, symbol rank => rank E}
)
