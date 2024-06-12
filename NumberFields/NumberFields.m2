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
   "isGalois",
   "splittingField",
   "compositums",
   "simpleExt",
   "asExtensionOfBase",--this probably shouldn't be exposed to the user long term
   "remakeField",--this probably shouldn't be exposed to the user long term
   "minimalPolynomial",
   "vectorSpace",
   "ringMapFromMatrix",
   "isFieldAutomorphism",
   "matrixFromRingMap"
};

NumberField = new Type of HashTable

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
remakeField = method(Options => {Variable=>null});
remakeField(Ring) := opts -> R1 -> (
    local var;
    local amb;
    local myIdeal;
    a := local a;
    R2 := (flattenRing R1)#0;
    if instance(R2, QuotientRing) then (amb = ambient R2; myIdeal = ideal R2) else (amb = R2; myIdeal = ideal(sub(0,R2)));
    numVars := #(gens amb);

    if opts.Variable === null then (var = a) else (var = opts.Variable);
    varList := {var_1..var_numVars};
    newRing1 := QQ(monoid[varList]);--we may want to think about the monomial order, we also might want to consider messing around with choosing a better presentation.
    phi := map(newRing1, amb, gens newRing1);

    newRing1/phi(myIdeal)
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

    new NumberField from {
        ring => toField outputRing, 
        pushFwd => myPushFwd,
        String => myStr,
        minimalPolynomial => genMinPolys, --a list of min polys of gens, for display purposes
        cache => new CacheTable from {degree => deg}
    }
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

ring(NumberField) := R1 -> (
    if (class (R1#ring) === PolynomialRing) then (
        return coefficientRing(R1#ring);
    );
    R1#ring
)




--degreeNF = method(Options => {})
degree(NumberField) := nf -> (
    if (nf#cache#?degree) then return nf#cache#degree;
    -*iota := map(ring(nf),QQ);
    rk := rank((pushFwd(iota))#0);
    nf#cache#degree = rk;
    rk*-
    --Karl:  something is wrong with pushFwd in this context, I rewrote this function for now.  The old version is above.
    degree ((ring nf)^1)
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
    if not (ring f1 === ring nf) then error "Expected an element of the NumberField";
    (nf#pushFwd#2)(f1)
);

gens(NumberField) := opts -> nf -> (
    gens ring nf
)


NumberFieldExtension = new Type of HashTable

numberFieldExtension = method(Options => {})
numberFieldExtension(RingMap) := opts -> phi1 -> (
    new NumberFieldExtension from {
        source=>numberField source phi1, 
        target=>numberField target phi1, 
        "map"=>phi1, 
        cache => new CacheTable from {}
    }
);

numberFieldExtension(RingElement) := opts -> f1 -> (
    if not (gens ring f1 == 1) then error "Expected a polynomial in a single variable";
    baseField := numberField coefficientRing ring f1;
    

);

source(NumberFieldExtension) := phi1 -> (phi1#source);
target(NumberFieldExtension) := phi1 -> (phi1#target);
map(NumberFieldExtension) := opts -> phi1 -> (phi1#"map");

degree(NumberFieldExtension) := nfe -> (
    if (nfe#cache#?degree) then return nfe#cache#degree;
    rk := rank((pushFwd(map nfe))#0);
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
splittingField = method(Options => {})
splittingField(RingElement) := opts -> f1 -> (
    --R1 := QQ[x];
    R1 := ring f1;
    S1 := R1;
    K1 := coefficientRing R1;
    K2 := coefficientRing R1;
    variableIndex := 1;
    finished := false;
    i := 1;
    while not finished do (
        L1 := decompose ideal f1;
        print i;
        print L1;
        i += 1;
        finished = true;
        executeForLoop := true;
        for i from 0 to #L1-1 do (
            if executeForLoop then (
                currentEntry := (entries gens L1#i)#0;
                if not (length(currentEntry)==1 and max(degree(currentEntry#0))==1) then (
                    finished = false;
                    K1 = R1/(L1#i);
                    S1 = K1[local a_variableIndex];
                    phi1 := map(S1, R1, {a_variableIndex});
                    f1 = phi1(f1);
                    R1 = S1;
                    variableIndex += 1;
                    executeForLoop = false;
                );
            );
        );
    );
    --K1
    --numberField K1
    --numberFieldExtension map((flattenRing K1)#0[local y], R1)
    numberFieldExtension (map(K1, K2))
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
    Node
        Key
            NumberFields
        Headline
            an example Macaulay2 package
        Description
            Text
                {\em FirstPackage} is a basic package to be used as an example.
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
