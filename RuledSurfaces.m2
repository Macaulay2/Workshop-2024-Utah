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
    AuxiliaryFiles => true
)
export {
    -- Types
    "ProjectiveBundle",
    -- methods
    "projectiveBundle",
}

needsPackage "Varieties"


ProjectiveBundle = new Type of HashTable

projectiveBundle = method()
projectiveBundle(ProjectiveVariety, CoherentSheaf) := (X, E) -> (
--checks on X and E
    if variety E =!= X then error "expected coherent sheaf on X";
    if rank E == 0 then error "expected a nonzero sheaf";
    if not isSmooth X then print "warning: X is not smooth and results may be incorrect";
--TODO: add option to bypass local freeness check and pruning
    if not isLocallyFree E then return "error: expected a locally free sheaf";
    new ProjectiveBundle from{symbol Variety => X, symbol CoherentSheaf => E, symbol rank => rank E}
)

load "./RuledSurfaces/RS-tests.m2"

