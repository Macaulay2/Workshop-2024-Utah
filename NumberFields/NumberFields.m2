newPackage(
       "NumberFields",
	Version => "0.0.0", 
        Date => "June 1, 2024",
    	Authors => {{}},
    	Headline => "number fields",
	PackageImports => {}, 
	PackageExports => {"PushForward"}
	)

export{
   "NumberField", 
   "numberField",
   "NumberFieldExtension",
   "numberFieldExtension",
   "isGalois",
};

NumberField = new Type of Ring

numberField = method(Options => {})
numberField(RingElement) := opts -> f1 -> (
	R1 := ring f1;

	if not isField coefficientRing R1 then error("Expected a polynomial over a field.");
	if #(gens R1) != 1 then error("Expected a polynomial in one variable.");
	if char R1 != 0 then error("Expected characteristic 0.");

	-- Verifies that the resulting quotient is a field.
	if f1 == 0 then error("Expected nonzero polynomial.");
	if not isPrime ideal(f1) then error("Expected an irreducible polynomial.");

	new NumberField from toField (R1/ideal(f1))
)

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
