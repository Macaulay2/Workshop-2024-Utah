newPackage(
    "RuledSurfaces",
    Date     => "2 Jun 2024",
    Version  => "0.1",
    Keywords => { "TBD" },
    Headline => "TBD",
    Authors  => {
	{   Name => "Devlin Mallory,",
	    Email => "malloryd@math.utah.edu",
	    HomePage => "https://www.math.utah.edu/~malloryd/" }
        },
    Authors  => {
	{   Name => "Fei Xiang,",
	    Email => "fxiang2@uci.edu",
	    HomePage => "" }
        },
    {   Name => "Swaraj Pande,",
	    Email => "swarajsp@umich.edu",
	    HomePage => "" }
        },
    AuxiliaryFiles => true,
    Reload => true
)
export {
    -- Types
    "ProjectiveBundle",
    -- methods
    "projectiveBundle",
    "findGlobalGeneratorsOfTwist",
}

needsPackage "Varieties"

ProjectiveBundle = new Type of HashTable

projectiveBundle = method()
projectiveBundle(ProjectiveVariety, CoherentSheaf) := (X, E) -> (
--checks on X and E
    if variety E =!= X then error "expected a coherent sheaf on X";
    if rank E == 0 then error "expected a nonzero sheaf";
    if not isSmooth X then print "warning: X is not smooth and results may be incorrect";
--TODO: add option to bypass local freeness check and pruning
    if not isLocallyFree E then error "expected a locally free sheaf";
    new ProjectiveBundle from{symbol Variety => X, symbol CoherentSheaf => E, symbol rank => rank E}
)


findGlobalGeneratorsOfTwist = method()
findGlobalGeneratorsOfTwist(CoherentSheaf) := E -> (
    i := regularity module E;
    while sheaf coker basis(i, module E) == 0 do i = i - 1;
    (i + 1, sheaf basis(i + 1, module E))
)



--load "./RuledSurfaces/RS-tests.m2"





