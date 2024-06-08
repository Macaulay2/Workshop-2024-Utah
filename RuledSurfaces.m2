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
    Reload => true,
    DebuggingMode => true
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
    "sectionFromLineBundleQuotient",
    "isMultipleOf",
    "minimalEmbedding",
    "checkSurjectivityOnH0",
    "planeTangoCurve",
    "frobeniusSheafMap",
    "minVeryAmpleTwist",
}

protect symbol Bound

needsPackage "Varieties"
needsPackage "Divisor"
needsPackage "DirectSummands"

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
generatingSurjection = method()
generatingSurjection ProjectiveBundle := PE -> last PE.generators

isLineBundle = method()
isLineBundle(CoherentSheaf) := L -> isLocallyFree L and rank L == 1

findGlobalGeneratorsOfTwist = method()
findGlobalGeneratorsOfTwist(ProjectiveBundle) := PE -> if not PE.?generators then PE.generators = findGlobalGeneratorsOfTwist(sheaf PE)
findGlobalGeneratorsOfTwist(CoherentSheaf) := E -> (
    i := regularity module E;
    --TODO: should we also check that this map induces a surjection on H^0?
    while sheaf coker basis(i, module E) == 0 do i = i - 1;
    (i + 1, sheaf basis(i + 1, module E))
)

checkSurjectivityOnH0 = method()
checkSurjectivityOnH0(ProjectiveBundle, WeilDivisor, ZZ) := (PE, D, m) -> checkSurjectivityOnH0(PE, divisorToLineBundle D, m)
checkSurjectivityOnH0(ProjectiveBundle, CoherentSheaf, ZZ) := (PE, L, m) -> (
    q := last findGlobalGeneratorsOfTwist sheaf PE;
    --TODO: erase the first line below once symmetric powers work better
    if m == 1 then isSurjective HH^0(q ** L) else
    isSurjective HH^0(symmetricPower(m, q) ** L)
)


multiProjEmbedding = method()
multiProjEmbedding(ProjectiveBundle) := PE -> multiProjEmbedding(sheaf PE)
multiProjEmbedding(CoherentSheaf) := (cacheValue symbol ideal) (E -> (
    q := matrix last findGlobalGeneratorsOfTwist E;
    quotient ker symmetricAlgebra q)
)


lineBundleToDivisor = method()
lineBundleToDivisor(CoherentSheaf) := L -> (
    --TODO: add checks that L really is a line bundle
    if not isLineBundle L then error "expected a line bundle";
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
    if ring D =!= ring variety PE then error "expected a divisor on base variety of bundle";
    --first do the map phi_D x id_P^r:
    --if OO_X(-a) -> E is our surjection, global sections of E(a) live in degree {1,a}; i.e., the linear series p^*OO_X(n) ** OO_PE(m) should live in degree {m, n - m * a }.
    findGlobalGeneratorsOfTwist PE;
    a := generatingDegree PE;
    X := variety PE;
    xx := (ring X)_0;
    E := sheaf PE;
    T := multiProjEmbedding(E);
    T' := prune T;
    --if the divisor D is a multiple of OO_X(1) then just use basis for T':
    d0 := isMultipleOf(divisorToLineBundle D, OO_X(1));
    if instance(d0, ZZ) then(
        B := basis({m, d0}, T');
        b := rank source B;
        t := symbol t;
        return Proj quotient ker map(T', (coefficientRing ring X)[t_0..t_(b-1)], gens image B));
    --if not, map T using D + d H with d minimal
    (d, LDH) := minVeryAmpleTwist divisorToLineBundle D;
    DH := lineBundleToDivisor LDH;
    phiD := mapToProjectiveSpace DH;
    SD := source phiD;
    --phi is phiDH x id_P^r
    phi := map(T, SD monoid T, vars T | matrix phiD);
    T' = prune quotient ker phi;
    --if D + dH divides OO_X(ma) (this is rare), returns OO_X(ma)/(D+dH):
    multipleFlag := isMultipleOf(OO_X(m * a), divisorToLineBundle(DH));
    if false and instance(multipleFlag, ZZ) then (
        B = basis({m, 1 + multipleFlag}, T');
        b = rank source B;
        T = symbol T;
        Proj quotient ker map(T', (coefficientRing ring X)[T_0..T_(b-1)], gens image B))
    else (
        --pH is the preimage of the twisting divisor H = OO_X(1) under phiD x id_P^r
        --pDa is the preimage of the twisting divisor Da = OO_X(a) under phiD x id_P^r
        --pD is the preimage of D+d H
        pH := divisor sub(preimage_phi sub(ideal xx, T), T');
        pD := divisor ideal last gens T';
        --the mapping divisor below corresponds to "basis({m,m/deg D *a + m},T')
        mappingDivisor := - (m * a + d) * pH + pD + m * divisor(T'_0);
        Proj quotient ker mapToProjectiveSpace mappingDivisor)
)

--only "minimal" relative to the choice of OO_X(1)
minimalEmbedding = method()
minimalEmbedding(CoherentSheaf)    := E -> minimalEmbedding projectiveBundle E
minimalEmbedding(ProjectiveBundle) := PE -> (
    findGlobalGeneratorsOfTwist PE;
    a := generatingDegree PE;
    imageOfLinearSeries(PE, OO_(variety PE)(a + 1), 1)
)

--checks whether L2 divides L1, i.e., if L1 is L2^n for some n
isMultipleOf = method(Options => {Bound => 5})
isMultipleOf(CoherentSheaf,CoherentSheaf) := opts -> (L1, L2) -> (
    if not (isLineBundle L1 and isLineBundle L2) then error "expected two line bundles";
    (M1, M2) := (module prune L1, module prune L2);
    if not isFreeModule M1 and isFreeModule M2 then return false;
    if isFreeModule M1 and isFreeModule M2 then(
        if regularity M1 % regularity M2 == 0 then return regularity M1//regularity M2 else return false
    );
    (D1, D2) := (lineBundleToDivisor L1, lineBundleToDivisor L2);
    b := opts.Bound;
    for i from -b to b do(
        if isLinearEquivalent(D1, i*D2, IsGraded => true) then return i
    )
)


--alternative to the method below, but seems to be slower
--minVeryAmpleTwist = method()
--minVeryAmpleTwist(CoherentSheaf) := (cacheValue symbol min) (D -> (
--    i := 0;
--    while not isVeryAmple lineBundleToDivisor D(i) do i = i + 1;
--    (i, D(i))
--)
--)

minVeryAmpleTwist = method()
minVeryAmpleTwist(CoherentSheaf) := (cacheValue symbol min) (L -> (
    D := lineBundleToDivisor L;
    H := divisor (ring D)_0;
    i := 0;
    while not isVeryAmple(D + i * H)  do i = i + 1;
    (i, divisorToLineBundle(D + i * H))
)
)


sectionFromLineBundleQuotient = method()
sectionFromLineBundleQuotient(SheafMap) := beta -> sectionFromLineBundleQuotient(projectiveBundle source beta, beta)
sectionFromLineBundleQuotient(ProjectiveBundle, SheafMap) := (PE, beta) -> (
    if not isSurjective beta then error "expected a surjection";
    E := source beta;
    L := target beta;
    if not E == sheaf PE then error "expected source to be bundle on PE";
    if not isLineBundle L then error "expected a line bundle";
    findGlobalGeneratorsOfTwist PE;
    qtilde := beta * generatingSurjection PE;
    ker map(quotient ker symmetricAlgebra matrix qtilde, multiProjEmbedding E)
)


--load "./RuledSurfaces/RS-tests.m2"
load "./RuledSurfaces/TangoCurves.m2"





