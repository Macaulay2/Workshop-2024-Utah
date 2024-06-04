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
    "generatingDegree",
    "findGlobalGeneratorsOfTwist",
    "multiProjEmbedding",
    "lineBundleToDivisor",
    "divisorToLineBundle",
    "imageOfLinearSeries",
}

needsPackage "Varieties"
needsPackage "Divisor"

ProjectiveBundle = new Type of MutableHashTable
LineBundleOnProjectiveBundle = new Type of HashTable

projectiveBundle = method()
projectiveBundle(CoherentSheaf) := E -> projectiveBundle(variety E, E)

projectiveBundle(ProjectiveVariety, CoherentSheaf) := (X, E) -> (
--checks on X and E
    if variety E =!= X then error "expected a coherent sheaf on X";
    if rank E == 0 then error "expected a nonzero sheaf";
    if not isSmooth X then print "warning: X is not smooth and results may be incorrect";
--TODO: add option to bypass local freeness check and pruning
    if not isLocallyFree E then error "expected a locally free sheaf";
    new ProjectiveBundle from{
        symbol variety => X, 
        symbol CoherentSheaf => E, 
        symbol rank => rank E,
    }
)

variety ProjectiveBundle := PE -> PE.variety
sheaf   ProjectiveBundle := PE -> PE.CoherentSheaf
generatingDegree = method()
generatingDegree ProjectiveBundle := PE -> first PE.generators

findGlobalGeneratorsOfTwist = method()
findGlobalGeneratorsOfTwist(ProjectiveBundle) := PE -> if not PE.?generators then PE.generators = findGlobalGeneratorsOfTwist(sheaf PE)
findGlobalGeneratorsOfTwist(CoherentSheaf) := E -> (
    i := regularity module E;
    while sheaf coker basis(i, module E) == 0 do i = i - 1;
    (i + 1, sheaf basis(i + 1, module E))
)

multiProjEmbedding = method()
multiProjEmbedding(ProjectiveBundle) := PE -> multiProjEmbedding(sheaf PE)
multiProjEmbedding(CoherentSheaf) := E -> (
    q := matrix last findGlobalGeneratorsOfTwist E;
    quotient ker symmetricAlgebra q
)


lineBundleToDivisor = method()
lineBundleToDivisor(CoherentSheaf) := L -> (
    --TODO: add checks that L really is a line bundle
    if not (isLocallyFree L and rank L == 1) then error "expected a line bundle";
    X := variety L;
    H := Hom(OO_X^1,L);
    --if we have OO_X -> L, then get L^-1 -> OO_X, and the image is an ideal representing the divisor
    if H != 0 then divisor ideal image matrix dual homomorphism H_0 else(
        a := first min degrees L;
        Ha := Hom(OO_X^1,L(a));
        divisor(ideal image matrix dual homomorphism Ha_0) - a * divisor( (ring X)_0)
    )
)

divisorToLineBundle = method()
divisorToLineBundle(WeilDivisor) := D -> prune sheaf OO(D)

--given a projective bundle p: P(E) -> X, a divisor D on X, and an integer m, finds the closure image of P(E) under the rational map corresponding to p^* O_X(D) ** OO_P(e)(m)
imageOfLinearSeries = method()
imageOfLinearSeries(ProjectiveBundle, CoherentSheaf, ZZ) := (PE, L, m) -> imageOfLinearSeries(PE, lineBundleToDivisor L, m)

imageOfLinearSeries(ProjectiveBundle, WeilDivisor, ZZ) := (PE, D, m) -> (
    --TODO: check that D lives on variety PE
    if ring D =!= ring variety PE then error "expected a divisor on base variety of bundle";
    --first do the map phi_D x id_P^r:
    --if OO_X(-a) -> E is our surjection, global sections of E(a) live in degree {1,a}; i.e., linear series p^*OO_X(n) ** OO_PE(m) should live in degree {m, ?}.
    findGlobalGeneratorsOfTwist PE;
    a := generatingDegree PE;
    Da := a*divisor( (ring variety PE)_0);
    E := sheaf PE;
    T := multiProjEmbedding(E);
    phiD := mapToProjectiveSpace(D);
    SD := source phiD;
    --phi is phiD x id_P^r
    phi := map(T, SD monoid T, vars T | matrix phiD);
    --preimDa is the preimage of the twisting divisor OO_X(a) under phiD x id_P^r
    preimDa := preimage_phi sub(ideal Da, T);
    T' := prune quotient ker phi;
    preimDa' := sub(preimDa, T');
    --the line below corresponds to "basis({m,1/a + m},T')
    Proj quotient ker mapToProjectiveSpace(-divisor preimDa' + m * divisor(T'_0) + divisor last gens T') 
)


--load "./RuledSurfaces/RS-tests.m2"





