newPackage(
       "NumberFields",
    Version => "0.0.0", 
        Date => "June 3, 2024",
        Authors => {
            {Name=>"Jack Garzella"},
            {Name=>"Nicholas Gaubatz", Email=>"nmg0029@auburn.edu", HomePage=>"https://nicholasgaubatz.github.io/"},
            {Name=>"Ethan Toshihiro Mowery"},
            {Name => "Karl Schwede", Email=>"schwede@math.utah.edu", HomePage=>"http://www.math.utah.edu/~schwede"}
        },
        Headline => "number fields",
        Keywords => {"field extension"},
    PackageImports => {}, 
    PackageExports => {"PushForward"},
    Reload => true,
    DebuggingMode => true
    )

export{
   "NumberField", 
   "numberField",
   "NumberFieldExtension",
   "numberFieldExtension",
   "isGalois",
--   "degreeNF",
   "splittingField"
};

NumberField = new Type of HashTable

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

    new NumberField from { ring => toField (R1/ideal(f1))}
)

numberField(Ring) := opts -> R1 -> (
    if R1===QQ then return new NumberField from {ring => R1};
    
    if not isPrime (ideal 0_R1) then error("Expected a field.");
    if not dim R1 == 0 then error("Expected a field.");
    
    if char R1 != 0 then error("Expected characteristic 0.");

    flattenedR1 := (flattenRing(R1))#0;
    iota := map(flattenedR1,QQ);
    try pushFwd(iota) else error("Not finite dimensional over QQ");

    new NumberField from {ring => toField (flattenedR1)}
)

internalNumberFieldConstructor := R1 -> (
    
);

--****************************
--NumberField basic operations
--****************************

ring(NumberField) := R1 -> (
    R1#ring
)




--degreeNF = method(Options => {})
--degreeNF(NumberField) := opts -> nf -> (
--
--    iota := map(nf,QQ);
--
--    degree((pushFwd(iota))#0)
--)

NumberFieldExtension = new Type of HashTable

numberFieldExtension = method(Options => {})
numberFieldExtension(RingMap) := opts -> phi1 -> (
    new NumberFieldExtension from {
        source=>numberField source phi1, 
        target=>numberField target phi1, 
        map=>phi1
    }
);

source(NumberFieldExtension) := phi1 -> (source phi1);
target(NumberFieldExtension) := phi1 -> (target phi1);
map(NumberFieldExtension) := phi1 -> (map phi1);


--*************************
--Methods
--*************************

isGalois = method(Options =>{})
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
                    S1 = K1[local x_variableIndex];
                    phi1 := map(S1, R1, {x_variableIndex});
                    f1 = phi1(f1);
                    R1 = S1;
                    variableIndex += 1;
                    executeForLoop = false;
                );
            );
        );
    );
    K1
    --R1
)

end

loadPackage ("NumberFields", Reload=>true)

--compositums = method(Options => {})
--compusitums(NumberField,NumberField) := opts -> (K1,K2) -> (
--    T := K1 ** K2;
--
--    -- compositums correspond to prime ideals
--    II := primaryDecomposition (ideal 0_T);
--
--    -- quotient rings
--    QRs := apply(II, I -> T / I);
--
--    -- make them number field objects?
--    NFs := apply(QRs, qr -> numberField(qr));
--
--    -- calculate degrees
--
--    -- sort by degree
--
--    -- get maps from K1 & K2
--
--)

-*
This is a comment block.
*-