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

numberField = method(Options => {})
numberField(RingElement) := opts -> f1 -> (
    R1 := ring f1;
    if not isField coefficientRing R1 then error("Expected a polynomial over a field.");
    if #(gens R1) != 1 then error("Expected a polynomial in one variable.");
    if char R1 != 0 then error("Expected characteristic 0.");

    -- Verifies that the resulting quotient is a field.
    if f1 == 0 then error("Expected nonzero polynomial.");
    print(f1);
    if not isPrime ideal(f1) then error("Expected an irreducible polynomial.");

    new NumberField from { 
        ring => toField (R1/ideal(f1)), 
        cache => new CacheTable from {}
    }
)

numberField(Ring) := opts -> R1 -> (
    if R1===QQ then return new NumberField from {ring => R1};
    
    if not isPrime (ideal 0_R1) then error("Expected a field.");
    if not dim R1 == 0 then error("Expected a field.");
    
    if char R1 != 0 then error("Expected characteristic 0.");

    flattenedR1 := (flattenRing(R1))#0;
    iota := map(flattenedR1,QQ);
    try pushFwd(iota) else error("Not finite dimensional over QQ");

    new NumberField from {
        ring => toField (flattenedR1), 
        cache => new CacheTable from {}
    }
)

internalNumberFieldConstructor := R1 -> (
    
);

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
    iota := map(ring(nf),QQ);
    rk := rank((pushFwd(iota))#0);
    nf#cache#degree = rk;
    rk
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
    while not finished do (
        L1 := decompose ideal f1;
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

-*
This is a comment block.
*-
end
