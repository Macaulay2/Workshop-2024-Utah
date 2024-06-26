newPackage(
       "NumberFields",
    Version => "0.0.0", 
        Date => "June 3, 2024",
        Authors => {
            {Name=>"Jack J Garzella", Email=>"jgarzell@ucsd.edu", HomePage=>"https://mathweb.ucsd.edu/~jjgarzel"},
            {Name=>"Nicholas Gaubatz", Email=>"nmg0029@auburn.edu", HomePage=>"https://nicholasgaubatz.github.io/"},
            {Name=>"Ethan Toshihiro Mowery"},
            {Name => "Karl Schwede", Email=>"schwede@math.utah.edu", HomePage=>"http://www.math.utah.edu/~schwede"}
        },
        Headline => "number fields",
        Keywords => {"field extension"},
    PackageImports => {},
    PackageExports => {"PushForward", "MinimalPrimes"},
    Reload => true,
    DebuggingMode => true
    )

export{
   "NumberField", 
   "numberField",
   "NumberFieldExtension",
   "numberFieldExtension",
   "TempNumberField",
   "tempNumberField",
   "isGalois",
   "splittingField",
   "compositums",
   "simpleExt",
   "getRoots",
   "ringElFromMatrix",
   "matrixFromRingEl",
   "asExtensionOfBase",--this probably shouldn't be exposed to the user long term
   "remakeField",--this probably shouldn't be exposed to the user long term
   "minimalPolynomial",
   "vectorSpace",
   --"ringMapFromMatrix",
   "isFieldAutomorphism",
   --"matrixFromRingMap"
};

NumberField = new Type of QuotientRing

-*
TempNumberField = new Type of QuotientRing

tempNumberField = method(Options => {Variable=>null})
tempNumberField(RingElement) := opts -> f1 -> (
    R1 := ring f1;
    if not isField coefficientRing R1 then error("Expected a polynomial over a field.");
    if #(gens R1) != 1 then error("Expected a polynomial in one variable.");
    if char R1 != 0 then error("Expected characteristic 0.");

    -- Verifies that the resulting quotient is a field.
    if f1 == 0 then error("Expected nonzero polynomial.");
    --if not isPrime ideal(f1) then error("Expected an irreducible polynomial.");

    tempNumberField(R1/ideal(f1), opts)
)



tempNumberField(Ring) := opts -> R1 -> (
    if R1===QQ then return new NumberField from {
            ring => R1, 
            pushFwd => pushFwd(map(QQ[],QQ)),
            cache => new CacheTable from {},
            String => "QQ, rational numbers"
        };
    
    if not isPrime (ideal 0_R1) then error("Expected a field.");
    if not dim R1 == 0 then error("Expected a field.");
    
    if char R1 != 0 then error("Expected characteristic 0.");
    outputRing := remakeField(R1, Variable=>opts.Variable);
    iota := map(outputRing,QQ);
    local myPushFwd;
    try myPushFwd = pushFwd(iota) else error("Not finite dimensional over QQ");
    genMinPolys := apply(gens outputRing, h->minimalPolynomial(h));

    genList := gens outputRing;
    genCt := #genList;
    i := 0;
    deg := degree outputRing;
    
    myStr := "Number Field, degree ";
    myStr = myStr | toString(deg) | " / Q, generated by:";
    if genCt == 0 then myStr = myStr + " nothing";
    while (i < genCt) do (
        myStr = myStr | " " | toString(genList#i) | " (" | toString (genMinPolys#i) | ")";
        i = i+1;
        if (i < genCt) then myStr = myStr | ", ";
    );

    tempNumField := new TempNumberField from outputRing;        

    tempNumField#pushFwd = myPushFwd;
    tempNumField#String = myStr;
    tempNumField#minimalPolynomial = genMinPolys;
    tempNumField#cache#degree = deg;

    tempNumField
)
*-

--
--The entries should be:
--ring (the ring it points to)
--cache (cached data)
--

--****************************
--NumberField constructors
--****************************

--The following function takes a ring and renames the variables
--in the future, it might also try to clean up relations (prune, trim, our own combination of prune/trim)
remakeField = method(Options => {Variable=>null, Degree => 1, NoPrune => false});
remakeField(Ring) := opts -> R1 -> (
    local var;
    local amb;
    local myIdeal;
    local finalRing2;
    a := local a;
    R2 := (flattenRing R1)#0;
    if instance(R2, QuotientRing) then (amb = ambient R2; myIdeal = ideal R2) else (amb = R2; myIdeal = ideal(sub(0,R2)));
    numVars := #(gens amb);

    if opts.Variable === null then (var = a) else (var = opts.Variable);
    varList := {var_1..var_numVars};
    varDegrees := apply(numVars, i->opts.Degree);
    newRing1 := QQ(monoid[varList, Degrees=>varDegrees]);--we may want to think about the monomial order, we also might want to consider messing around with choosing a better presentation.
    phi := map(newRing1, amb, gens newRing1);
    finalRing := newRing1/phi(myIdeal);    
    psi := map(finalRing, R1, matrix phi);
    if not opts.NoPrune then (
        finalRing2 = prune finalRing;
        psi = (finalRing.minimalPresentationMap)*psi;
    )
    else(
        finalRing2 = finalRing;
    );

    (finalRing2, psi)
)

numberField = method(Options => {Variable=>null})
numberField(RingElement) := opts -> f1 -> (
    R1 := ring f1;
    if not isField coefficientRing R1 then error("Expected a polynomial over a field.");
    if #(gens R1) != 1 then error("Expected a polynomial in one variable.");
    if char R1 != 0 then error("Expected characteristic 0.");

    -- Verifies that the resulting quotient is a field.
    if f1 == 0 then error("Expected nonzero polynomial.");
    --if not isPrime ideal(f1) then error("Expected an irreducible polynomial.");

    numberField(R1/ideal(f1), opts)
)



numberField(Ring) := opts -> R1 -> (
    if R1===QQ then return new NumberField from {
            ring => R1, 
            pushFwd => pushFwd(map(QQ[],QQ)),
            cache => new CacheTable from {},
            String => "QQ, rational numbers"
        };
    
    if not isPrime (ideal 0_R1) then error("Expected a field.");
    if not dim R1 == 0 then error("Expected a field.");
    
    if char R1 != 0 then error("Expected characteristic 0.");
    outputRing := (remakeField(R1, Variable=>opts.Variable))#0;
    iota := map(outputRing,QQ);
    local myPushFwd;
    try myPushFwd = pushFwd(iota) else error("Not finite dimensional over QQ");
    if not isFreeModule(myPushFwd#0) then error "numberField: something went wrong, this should be a free module over QQ";
    genMinPolys := apply(gens outputRing, h->minimalPolynomial(h));

    genList := gens outputRing;
    genCt := #genList;
    i := 0;
    deg := rank(myPushFwd#0);
    
    myStr := "Number Field, degree ";
    myStr = myStr | toString(deg) | " / Q, generated by:";
    if genCt == 0 then myStr = myStr | " nothing";
    while (i < genCt) do (
        myStr = myStr | " " | toString(genList#i) | " (" | toString (genMinPolys#i) | ")";
        i = i+1;
        if (i < genCt) then myStr = myStr | ", ";
    );

    -*
    new NumberField from {
        ring => toField outputRing, 
        pushFwd => myPushFwd,
        String => myStr,
        minimalPolynomial => genMinPolys, --a list of min polys of gens, for display purposes
        cache => new CacheTable from {degree => deg}
    }
    *-
    tempNumField := new NumberField from outputRing;        

    tempNumField#pushFwd = myPushFwd;
    tempNumField#String = myStr;
    tempNumField#minimalPolynomial = genMinPolys;
    tempNumField#cache#degree = deg;

    tempNumField
)

internalNumberFieldConstructor := R1 -> (
    
);

--*****************************
--NumberField display
--*****************************


net NumberField := nf -> (
    -*
    myStr := "Number Field, degree ";
    myStr = myStr | toString(degree nf) | " / Q, generated by:";
    i := 0;
    genCt := #(gens ring nf);
    genList := gens ring nf;
    if genCt == 0 then myStr = myStr + " nothing";
    while (i < genCt) do (
        myStr = myStr | " " | toString(genList#i) | " (" | toString (nf#minimalPolynomial#i) | ")";
        i = i+1;
        if (i < genCt) then myStr = myStr + ", ";
    );
    --do some code
    *-

    nf#String
)

--****************************
--NumberField basic operations
--****************************

-*ring(NumberField) := R1 -> (
    if (class (R1#ring) === PolynomialRing) then (
        return coefficientRing(R1#ring);
    );
    R1#ring
    varList := gens R1;
    myMon := monoid[varList];
    amb := QQ[varList];
    amb/sub(ideal R1, amb)
)*-




--degreeNF = method(Options => {})
degree(NumberField) := nf -> (
    if (nf#cache#?degree) then return nf#cache#degree;
    -*iota := map(ring(nf),QQ);
    rk := rank((pushFwd(iota))#0);
    nf#cache#degree = rk;
    rk*-
    --Karl:  something is wrong with pushFwd in this context, I rewrote this function for now.  The old version is above.
    rank(nf#pushFwd#0)
)

--this gives the basis for the numberField over QQ
basis(NumberField) := opts -> nf -> (
    first entries (nf#pushFwd#1)
);

vectorSpace = method(Options=>{})
vectorSpace(NumberField) := opts -> nf -> (
    nf#pushFwd#0
)

vector(RingElement, NumberField) := (f1, nf) -> (
    if not (ring f1 === nf) then error "Expected an element of the NumberField";
    (nf#pushFwd#2)(f1)
);



NumberFieldExtension = new Type of RingMap

numberFieldExtension = method(Options => {})
numberFieldExtension(RingMap) := opts -> phi1 -> (
    answer := new NumberFieldExtension from phi1;
    answer#cache#String = ("Number field extension, degree " | degree answer);
    answer

    -*new NumberFieldExtension from {
        source=>numberField source phi1, 
        target=>numberField target phi1, 
        "map"=>phi1, 
        cache => new CacheTable from {}
    }*-
);

net NumberFieldExtension := nfe -> (nfe#cache#String)

numberFieldExtension(RingElement) := opts -> f1 -> (
    if not (gens ring f1 == 1) then error "Expected a polynomial in a single variable";
    baseField := numberField coefficientRing ring f1;
    

);

--source(NumberFieldExtension) := phi1 -> (source phi1);
--target(NumberFieldExtension) := phi1 -> (target phi1);
map(NumberFieldExtension) := opts -> phi1 -> (phi1);

degree(NumberFieldExtension) := nfe -> (
    if (nfe#cache#?degree) then return nfe#cache#degree;
    rk := rank((pushFwd(map(target nfe, source nfe, matrix nfe)))#0);
    nfe#cache#degree = rk;
    rk
)

norm(RingElement) := (elt) ->(
    S := ring elt;
    return det pushFwd(map(S^1, S^1, {{elt}}));
);
trace(RingElement) := (elt) -> (
    S := ring elt;
    return trace pushFwd(map(S^1, S^1, {{elt}}));
);
--*************************
--Methods
--*************************

isGalois = method(Options =>{})
isGalois(NumberField) := opts -> K -> {
    mapList := compositums(K,K);
    degs := apply(mapList, x -> x#3);
    L := all(degs, d -> d == degs#0);
    L
}

isGalois(NumberFieldExtension) := opts -> iota -> (
     myMapList := {}; --replace with Jack's function when ready
    --assuming iota : K -> L, myMapList is a list of maps L -> L_i where 
    --L_i is one of the components of L **_K L.
    
    --check if all degrees are the same, and equal to 1.
)

isGalois(RingMap) := opts -> iota -> (
   isGalois(numberFieldExtension iota)
)

-- splittingField method
--****KARL:  THIS IS CURRENTLY BROKEN, I TRIED TO MAKE IT FASTER...*****
splittingField = method(Options => {Variable=>null, Verbose=>false})
splittingField(RingElement) := opts -> f1 -> (
    --R1 := QQ[x];
    R1 := ring f1;
    curf1 := f1;
    curf1old := f1;
    varName := gens R1;
    if not (#varName == 1) then error "Expected a polynomial ring in a single variable";
    S1 := R1;
    Svar := (gens R1)#0;
    SvarOld := Svar;
    K1 := coefficientRing R1;
    --K2 := (remakeField( coefficientRing R1, Degree=>0))#0;
    psi := map(K1, K1);
    totalPsi := psi;
    local unMadeField;
    local phi1;
    local linTerm;
    local finalAnswer;
    local flatTargetRing;
    local flatPsi;
    local newPsi;
    local kappa;
    variableIndex := 1;
    finished := false;
    i := 1;
    idealList := {ideal curf1};
    L1 := flatten apply(idealList, z->decompose z);
    local var;
    a := local a;
    if opts.Variable === null then (var = a) else (var = opts.Variable);
    while not finished do (
        print ("Starting a loop : " | toString(idealList));        
        idealList = select(idealList, z->not isLinear z);        --let's only keep the good ones.
        --print i;
        --print curf1;
        if opts.Verbose then print idealList;
        i += 1;
        finished = true;
        --executeForLoop := true;
        if (#idealList > 0) then (
            curIdeal := idealList#0;
            currentEntry := (entries gens curIdeal)#0;   --grab a polynomial to work with
            finished = false;
            newTargetRing := K1[var_variableIndex];
            (flatTargetRing,flatPsi) = flattenRing(newTargetRing);
            newPsi = flatPsi*map(newTargetRing, S1, gens newTargetRing);
            kappa = map(S1, K1);
            K1 = flatTargetRing/newPsi(curIdeal);
            psi = newPsi * kappa;
            1/0;
            -*unMadeField = R1/(idealList#0);
            totalPsi = (map(unMadeField, target totalPsi))  * totalPsi;
            (K1, psi) = flattenRing(K1[local a_variableIndex]);*-
            
            --(K1, psi) = remakeField (unMadeField, Degree=>0, NoPrune=>true);                    
            totalPsi = psi*totalPsi;
            --S1 = K1[local a_variableIndex];                    
            S1 = K1[varName];
            SvarOld = Svar;
            Svar = sub(varName#0, S1);
            linTerm = Svar - newPsi(SvarOld);
            phi1 = map(S1, R1, {Svar});    --this is behaving badly, let me try sub
            print phi1;               
            --curf1old = phi1(curf1);
            --curf1 = curf1old;-- // linTerm;      --is this working? --it is not.
            --assert(linTerm*curf1 == curf1old);
            
            idealList = drop(apply(idealList, z->sub(z, S1)), 1);
            idealList = {saturate(sub(curIdeal, S1), linTerm)} | idealList;
            R1 = S1;
            variableIndex += 1;
        )

        
        -*for i from 0 to #idealList-1 do (
            if executeForLoop then (
                currentEntry := (entries gens L1#i)#0;         
                if opts.Verbose then print (toString(currentEntry) | " : " | toString(length(currentEntry)) | "," | toString(degree(currentEntry#0)) );
                if not (length(currentEntry)==1 and max(degree(currentEntry#0))==1) then (
                    finished = false;
                    unMadeField = R1/(L1#i);
                    totalPsi = (map(unMadeField, target totalPsi))  * totalPsi;
                    (K1, psi) = remakeField (unMadeField, Degree=>0);                    
                    totalPsi = psi*totalPsi;
                    --S1 = K1[local a_variableIndex];                    
                    S1 = K1[varName];
                    SvarOld = Svar;
                    Svar = sub(varName#0, S1);
                    linTerm = Svar - psi(SvarOld);
                    phi1 = map(S1, R1, {Svar});                    
                    curf1old = phi1(curf1);
                    curf1 = curf1old;-- // linTerm;      --is this working? --it is not.
                    --assert(linTerm*curf1 == curf1old);
                    R1 = S1;
                    variableIndex += 1;
                    executeForLoop = false;
                );                                
            );            
        );*-
    );
    --K1
    --numberField K1
    --numberFieldExtension map((flattenRing K1)#0[local y], R1)
    --numberFieldExtension (map(K1, K2))    
    (finalAnswer, psi) = remakeField(K1, Degree=>1, Variable=>opts.Variable);

    (numberField finalAnswer, numberFieldExtension(psi*totalPsi))
)

isLinear = method(Options=>{})
isLinear(Ideal) := opts -> (J1) -> (
    idealGens := (entries gens J1)#0;
    length(idealGens)<=1 and max(degree(idealGens#0))<=1
)

syntheticDivision = method(Options=>{})
syntheticDivision(RingElement, RingElement) := (f1, g1) -> ( --compute f1 / g1, where g1 = x-a, and where g1 divides f1

)

isFieldAutomorphism = method(Options=>{})

isFieldAutomorphism(NumberField, Matrix) := opts -> (NF1, sigma1) -> (
    R1 := ring NF1;
    C1 := coefficientRing R1;
    P1 := (pushFwd(map(R1, C1)))#2;
    phi1 := ringMapFromMatrix(NF1, sigma1);
    if not isWellDefined phi1 then (return false;);
    if not isInjective phi1 then (return false;);
    newBasis1 := apply(basis NF1, i -> phi1(i));
    newGensAsBasis1 := apply(newBasis1, P1);
    A1 := newGensAsBasis1#0;
    for i from 1 to #newGensAsBasis1-1 do (
        A1 |= newGensAsBasis1#i;
    );
    return (A1-sigma1)==0 -- Test matrix equality only on entry level; A==sigma1 would return false if either sources or target differ
)

-- Helper method for isFieldAutomorphism
-*ringMapFromMatrix = method(Options=>{})

ringMapFromMatrix(NumberField, Matrix) := opts -> (NF1, sigma1) -> (
*-
ringMapFromMatrix = (NF1, sigma1) -> (
    R1 := ring NF1;
    C1 := coefficientRing R1;
    P1 := (pushFwd(map(R1, C1)))#2;
    gensAsBasis1 := apply(gens R1, P1); -- Expresses each generator of R1 as a vector w.r.t. the basis of NF1
    newGensAsBasis1 := apply(gensAsBasis1, i -> sigma1*i);
    newGensAsMatrices1 := apply(newGensAsBasis1, i -> matrix({basis NF1})*i);
    newGens1 := apply(newGensAsMatrices1, i -> (entries i)#0#0);
    map(R1, R1, newGens1)
)
-*
matrixFromRingMap = method(Options=>{})

matrixFromRingMap(NumberField, RingMap) := opts -> (n1, psi) ->(
    print "Hello";
    matrixFromRingMap(n1, n1, psi)
)
matrixFromRingMap(ZZ,ZZ) := opts -> (i,j) -> (
    i+j
)

matrixFromRingMap(NumberField, NumberField, Thing) := opts -> (nf1, nf2, psi) -> (
*-
matrixFromRingMap = (nf1, nf2, psi) -> (
    if not ((target psi === ring nf1) and (source psi === ring nf2)) then error "expected the map to go from the the second argument to the first";
    
    B1 := basis nf2;
    matrixOut := vector(psi(B1#0), nf1);
    i := 1;
    while (i < #B1) do (
        matrixOut = matrixOut | vector(psi(B1#i), nf1);
        i = i+1;
    );
    matrixOut
)
matrixFromRingEl = method(Options => {});
matrixFromRingEl(NumberField, RingElement) := opts -> (nF, rEl) -> (
    R := ring nF;
    return pushFwd(map(R^1, R^1, matrix{{rEl}}));
)
ringElFromMatrix = method(Options => {});
ringElFromMatrix(NumberField, Matrix) :=opts -> (nF, mat) -> (
    --We basically turn the natural linear algebra basis of our number field into a matrix, then row reduce it to turn mat into an element in our number field.
    R0 := ring nF;
    R1 := coefficientRing R0;
    M0 := (pushFwd(map(R0,R1)))_1;
    vList := {};
    for i from 0 to ((numgens source M0)-1) do(

        Mi := matrixFromRingEl(nF, (M0_i)_0);
        vi := vector reshape(R1^((numgens target Mi)*(numgens source Mi)), R1^1, Mi);
        vList = append(vList, vi);
    );
    vList = append(vList, vector reshape(R1^((numgens target mat)*(numgens source mat)), R1^1, mat));
    M := matrix(vList);
    RRM := reducedRowEchelonForm M;
    lastCol := RRM_{numColumns RRM-1};
    el := 0_R0;
    for i from 0 to ((numgens source M0)-1) do(
        el = el + (lastCol_0)_i * (M0_i)_0;
    );
    return el;
)
getRoots = method(Options =>{});
getRoots(RingElement) := opts -> (f) -> (
    R := ring f;
    (S,M, MInv) := (flattenRing (R,Result=>3));
    primeFactors := decompose ideal M(f);
    linearTerms := {};
    for i from 0 to ((length primeFactors)-1) do(
        if (degree primeFactors#i_0)#0 == 1 then (
            linearTerms = append(linearTerms, MInv(primeFactors#i_0));
        );
    );
    return linearTerms;
)

simpleExt = method(Options => {});
simpleExt(NumberField) := opts -> nf ->(
    --We first get the degree of K as a field extension over Q and store it as D. 
    K := ring nf;
    D := degree K;
    --We find an element that produces a degree D field extension.
    d := 0;
    c := 0;
    primitiveElement := 0; -- Uncomment along with below chunk to get a simpler primitive element
    while d < D do 
    (
        r := random(1, K); -- Get a random homogeneous RingElement from K1 of degree 1

        -- Uncomment the following along with primitiveElement := 0 above to get a simpler primitive element
        (if primitiveElement==0 then (
            primitiveElement = sum gens K;
         )
        else (
            primitiveElement += (random(gens K))#0; -- Randomly shuffles the list of generators of K and then takes the first element
         )
        );
        r = primitiveElement;
        --

        xx := local xx;
        R := QQ[xx];
        phi := map( K, R, {r});
        if  isPrime (kernel phi) then (
            I := kernel phi *sub (( 1/(((coefficients (first entries gens kernel phi)_0)_1)_0)_0), R);
            simpleExt := numberField(R / I);
            d = degree simpleExt;
        );
    );
    return simpleExt;
)

minimalPolynomial = method(Options => {})
minimalPolynomial(RingElement) := opts -> f1 -> (
    R1 := ring f1;
    S1 := (coefficientRing(R1))[local y];
    P1 := pushFwd(map(R1, coefficientRing(R1)));
    A1 := (P1#2)(1_R1);
    pow1 := 1;
    while gens(kernel(A1))==0 do (
        A1 = A1|(P1#2)(f1^pow1);
        pow1 += 1;
    );
    M1 := matrix({{1_S1}});
    for i1 from 1 to (pow1-1) do (
        M1 |= y^i1;
    );
    (entries (M1*(gens(kernel(A1)))))#0#0
)
minimalPolynomial(List) := opts -> L1 -> (
    apply(L1, i -> minimalPolynomial(i))
)

--********************************
--******Compositums
--********************************

asExtensionOfBase = method(Options => {})
asExtensionOfBase(NumberFieldExtension) := opts -> iota -> (
--
--    -- get source and target
--    s := ring(source iota);
--    t := ring(target iota);
--    -- get ideal from target
--    --I := ideal target; 
--    -- calculate numgens of ideal
--    n := numgens t;
--    -- create polynomial ring in numbgens variables
--    polyring := s[b_1..b_n];
--    -- create map by assigning generators to other generators
--    extensionMap := map(iota);
--    images := (entries  (matrix extensionMap))#0;
--    m := map(t,polyring,join(images, gens(t)));
--    -- take the kernel
--    k := kernel m;
--    -- return quotient by kernel
--    polyring / k
)

--loadPackage ("NumberFields", Reload=>true)

compositums = method(Options => {})
compositums(NumberField,NumberField) := opts -> (K1,K2) -> (
    T := ring(K1) ** ring(K2);
    -- compositums correspond to prime ideals
    II := decompose (ideal 0_T);
    -- quotient rings
    QRs := apply(II, I -> T / I);
    -- make them number field objects?
    NFs := apply(QRs, qr -> numberField(qr));
    -- sort by degree
    sorted := sort(NFs, degree);
    -- get maps from K1 & K2 
    K1maps := apply(sorted, nf -> map(ring(nf),ring(K1)));
    K2maps := apply(sorted, nf -> map(ring(nf),ring(K2)));

    degs := apply(sorted, degree);
    -- a slight hack to package the data
    inds := toList(0..(length(sorted)-1));
    infoList := apply(inds, i -> (sorted#i,K1maps#i,K2maps#i,degs#i));

    infoList
)


--*****************************
--Documentation
--*****************************
beginDocumentation()

doc ///
    Key
        NumberFields
    Headline
        number fields and Galois theory
    Description
        Text
            {\em NumberFields} is a package that allows for clean definition and use of number fields. A number field is a extension of $\QQ$ of finite degree, e.g., $\QQ(\sqrt(2))$, which is isomorphic to $\QQ[x]/(x^2-2)$.
///

doc ///
    Key
        splittingField
        (splittingField,RingElement)
    Headline
        creates the splitting field of a polynomial
    Usage
        splittingField f
    Inputs
        f: RingElement
            a polynomial whose splitting field the method returns
    Outputs
        : NumberFieldExtension
    Description
        Text
            This method creates a @TO NumberFieldExtension@ whose source is $\mathbb{Q}$, target is a NumberField consisting of the splitting field of $f$, and map is inclusion. Recall that the splitting field of a polynomial $f$ with coefficients in a field $\FF$ is the smallest field extension of $\FF$ in which $f$ splits, or factors into linear factors.
        Example
            R = QQ[x];
            f = x^2-2;
            splittingField f
            f = x^2+x+1;
            splittingField f
///

doc ///
    Key
        minimalPolynomial
        (minimalPolynomial, RingElement)
        (minimalPolynomial, List)
    Headline
        computes the minimalPolynomial of a field element or list of field elements
    Usage
        g = minimalPolynomial f
        L1 = minimalPolynomial L
    Inputs
        f: RingElement
            the polynomial whose minimal polynomial the method computes
        L: List
            a list of polynomials whose minimal polynomials the method computes
    Outputs
        g: RingElement
            minimal polynomial
        L1: List
            list of minimal polynomials
    Description
        Text
            This method computes minimal polynomials. Recall that if $\mathbb{K}$ is a field extension of $\mathbb{F}$, the minimal polynomial of a field element $\alpha$ in $\mathbb{K}$ is, if it exists, the unique monic polynomial $f$ of lowest degree with coefficients in $\mathbb{F}$ such that $\alpha$ is a root of $f$.
        Example
            R = QQ[x]/(x^3-2);
            f = x;
            g = minimalPolynomial f
            L = {x, x^2, x^2+x};
            L1 = minimalPolynomial L
///

--*****************************
--Tests
--*****************************

TEST /// --Test #0
    K = QQ[x]/ideal(x^3-2)
    L = K[y]/ideal(y^2 + y + 1)
    assert(degree numberField K == 3)
    assert(degree numberField L == 6)
///

-*TEST /// --Test #1
    K = QQ[x]
    f = x^2-2
    assert( isIsomorphic(target splittingField(f), numberField(QQ[y]/(y^2-2))) )
///

--compositums(NumberFieldExtenison,NumberFieldExtension) := opts -> (iota,kappa) -> (
--    -- check for common base
--    -- write iota and kappa as entensions of common base
--    -- take tensor product wrt common base
--    -- do all the same things as compositums over QQ
--
--)

-*
This is a comment block.
*-
end
