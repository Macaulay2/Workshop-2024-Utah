newPackage(
       "NumberFields",
    Version => "0.0.0", 
        Date => "June 3, 2024",
        Authors => {
            {Name=>"Jack Garzella"},
            {Name=>"Nicholas Gaubatz"},
            {Name=>"Ethan Toshihiro Mowery"},
            {Name => "Karl Schwede", Email=>"schwede@math.utah.edu", HomePage=>"http://www.math.utah.edu/~schwede"}
        },
        Headline => "number fields",
        Keywords => {"field extension"},
    PackageImports => {}, 
    PackageExports => {"PushForward"}
    )

export{
   "NumberField", 
   "numberField",
   "NumberFieldExtension",
   "numberFieldExtension",
   "isGalois",
--   "degreeNF",
};

NumberField = new Type of HashTable

--****************************
--NumberField constructors
--****************************

numberField = method(Options => {})
numberField(RingElement) := opts -> f1 -> (
    print("1");
    R1 := ring f1;
    print("2");
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

    
    if not isPrime (ideal 0_R1) then error("Expected a field.");
    if not dim R1 == 0 then error("Expected a field.");
    
    if char R1 != 0 then error("Expected characteristic 0.");

    flattenedR1 := (flattenRing(R1))#0;
    iota := map(flattenedR1,QQ);
    try pushFwd(iota) else error("Not finite dimensional over QQ");

    new NumberField from {ring => toField (flattenedR1)}
)

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