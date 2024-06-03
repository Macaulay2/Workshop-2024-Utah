newPackage(
    "RuledSurfaces",
    Date     => "2 Jun 2024",
    Version  => "0.1",
    Keywords => { "TBD" },
    Headline => "TBD",
    Authors  => {
	{   Name => "Devlin Mallory,",
	    Email => "malloryd@math.utah.edu",
	    HomePage => "https://www.math.utah.edu/~malloryd/" },
        {   Name => "Swaraj Pande,",
	    Email => "swarajsp@umich.edu",
	    HomePage => "" },
	{   Name => "Fei Xiang,",
	    Email => "fxiang2@uci.edu",
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
    "multiProjEmbedding",
}

needsPackage "Varieties"
needsPackage "Divisor"

ProjectiveBundle = new Type of HashTable

projectiveBundle = method()
projectiveBundle(CoherentSheaf) := E -> projectiveBundle(variety E, E)

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

multiProjEmbedding = method()
multiProjEmbedding(CoherentSheaf) := E -> (
    q := matrix last findGlobalGeneratorsOfTwist E;
    ker symmetricAlgebra q
)

--given a projective bundle p: P(E) -> X, a divisor D on X, and an integer m, finds the closure image of P(E) under the rational map corresponding to p^* O_X(D) ** OO_P(e)(m)
imageOfLinearSeries = method()
imageOfLinearSeries(ProjectiveBundle, Divisor, ZZ) := (PE, D, m) -> (
    --first do the map phi_D x id_P^r:
    --if OO_X(-a) -> E is our surjection, global sections of E(a) live in degree {1,a}; i.e., linear series p^*OO_X(n) ** OO_PE(m) should live in degree {m, n}.
    B:=basis({m,n},T);
    Proj quotient ker map(T, (coefficientRing T)[rank B], B)
)


--load "./RuledSurfaces/RS-tests.m2"





