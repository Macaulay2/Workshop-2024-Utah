--Loading packages
restart
loadPackage ("NumberFields", Reload=>true)
-- loadPackage ("Permutations", Reload=>true)
loadPackage ("TestNCAlgebra", Reload=>true)
loadPackage ("AssociativeAlgebras", Reload=>true)
stack apply(prefixPath, p -> p | Layout#1#"core" | "raw.m2")

----CORE STUFF
-- Core.Dictionary #? "rawQuotientRing"
-- importFrom_Core{"rawQuotientRing"}
-- importFrom(Core, "raw")
-- importFrom(Core, "sin")
-- importFrom(Core, "rawQuotientRing")
-- Play around with creating my own ring. What does this even mean?
A = ZZ{x,y,z}
f = y*z + z*y - x^2
-----------------------------------------------------------------------------------------------------------------------
-- NumberField = new Type of EngineRing

-- debug Core
-- presentation NumberField := (R) -> presentation R
-- generators NumberField := opts -> (R) -> generators R
-- isCommutative NumberField := (R) -> isCommutative R
-- degreeLength NumberField := (R) -> degreeLength R
-- numgens NumberField := (R) -> numgens R

R = QQ[x] / ideal(x^2-2)
NF = new NumberField from R 
-- new NumberField from QuotientRing := (QuotientRing, inits) -> new QuotientRing of RingElement from new HashTable from init
new NumberField from List := (QuotientRing, inits) -> new QuotientRing of RingElement from new HashTable from inits

NF = new NumberField from R 
--------------------------------------------------------
loadPackage ("testLocalRings", Reload=>true)
S = testlocalRing(QQ[x], ideal(x^2-2))
use S  
x 
loadPackage ("LocalRings", Reload=>true)
S = localRing(QQ[x], ideal(x^2-2))
use S 
x 
S = QQ[x]/(x^2-2)
S.generatorSymbols;

loadPackage ("testNumberField", Reload=>true)
R = QQ[x] / ideal(x^2-2)
S = testnumberField(R)

T0 = QQ[x]
S = testnumberField(T0)

T1 = QQ[x]
x --This is from T0
use R 
x  -- This is from R
use T0 
x -- This is from T0
S = testnumberField(R)
A0=QQ[x]
x --This is from S
use A0
x -- This is STILL from S
use R 
x -- This is from R as desired
use A0
x -- This is from R!?
use S 
x -- This is from S as expected


A0 = QQ[x]
A1 = A0/ideal(x^2)
A2 = QQ[x]
x

symbol x
R.generators 
R.?generatorSymbols
R.?generatorExpressions
R.?index        
R.?indexStrings 
R.?indexSymbols 
expression S := lookup(expression,R);
NumberField= new Type of EngineRing;
NumberField.synonym = "Number Field"

-------------------------------------------------------------
TestRing = new Type of Ring;



TestQuotientRing = new Type of TestRing;
TestPolynomialRing = new Type of TestRing;
TestRingElement = new Type of HashTable
TestMatrix = new Type of MutableHashTable
TestMonomial = new Type of HashTable
globalAssignment TestRing;

coefficientRing TestRing := A -> A.CoefficientRing

generators TestRing := opts -> A -> (
    if A.?generators then A.generators else {}
)

numgens TestRing := A -> #(A.generators)
isHomogeneous TestRing := A -> (
   if instance(A,TestPolynomialRing) then true
   else isHomogeneous ideal A
)

testMonomial = method()
testMonomial (List,TestRing) := (monL,B) -> (
   newMon := new TestMonomial from {(symbol monList) => monL,
                                  (symbol ring) => B};
   newMon
)
putInRing = method()
putInRing (TestMonomial, ZZ) := 
putInRing (TestMonomial, QQ) :=
putInRing (TestMonomial, RingElement) := (mon,coeff) ->
      new (mon.ring) from {(symbol ring) => mon.ring,
        (symbol cache) => new CacheTable from {("isReduced",false)},
                           (symbol terms) => new HashTable from {(mon,promote(coeff,coefficientRing (mon.ring)))}}
putInRing (List, TestRing, ZZ) := 
putInRing (List, TestRing, QQ) :=
putInRing (List, TestRing, RingElement) := (monList,A,coeff) -> (
    mon := testMonomial(monList,A);
    new A from {(symbol ring) => A,
                (symbol cache) => new CacheTable from {("isReduced",false)},
                (symbol terms) => new HashTable from {(mon,promote(coeff,coefficientRing A))}}
)

use TestRing := A -> (scan(A.generatorSymbols, A.generators, (sym,val) -> sym <- val); A)

--Polynomial testring
new TestPolynomialRing from List := (TestPolynomialRing, inits) -> new TestPolynomialRing of TestRingElement from new HashTable from inits
Ring List := (R, varList) -> (
   -- get the symbols associated to the list that is passed in, in case the variables have been used earlier.
   print("BOX 1");
   if #varList == 0 then error "Expected at least one variable.";
   if #varList == 1 and class varList#0 === Sequence then varList = toList first varList;
   varList = varList / baseName;
   print("BOX 2");
   A := new TestPolynomialRing from {(symbol generators) => {},
                                   (symbol generatorSymbols) => varList,
                                   (symbol degreesRing) => degreesRing 1,
				   (symbol CoefficientRing) => R,
                                   (symbol cache) => new CacheTable from {},
				   (symbol baseRings) => {ZZ},
                                   BergmanRing => false};
   newGens := apply(varList, v -> v <- putInRing({v},A,1));
    print(newGens);
--    if R === QQ or R === ZZ/(char R) then A#BergmanRing = true;

   A#(symbol generators) = newGens;
   
--    setWeights(A, toList (#(gens A):1));
      
--    --- all these promotes will need to be written between this ring and all base rings.
--    promote (ZZ,A) := (n,A) -> putInRing({},A,promote(n,R));


--    promote (QQ,A) := (n,A) ->  putInRing({},A,promote(n,R));
        
--    promote (R,A) := (n,A) -> putInRing({},A,n);
      
--    promote (TestMatrix,A) := (M,A) -> (
--        if M.source == {} or M.target == {} then
--           TestMatrix(A,M.target,M.source)
--        else (
--           prom := TestMatrix apply(M.matrix, row -> apply(row, entry -> promote(entry,A)));
--           if isHomogeneous M then
--              assignDegrees(prom,M.target,M.source);
--           prom
--        )
--    );

--    promote (A,A) := (f,A) -> f;
   
--    addVals := (c,d) -> (
--       e := c+d;
--       if e == 0 then continue else e
--    );

--    multVals := (c,d) -> c*d;
      
--    multKeys := (m,n) -> m | n;

--    A + A := (f,g) -> (
--       -- new way
--       newHash := removeZeroes merge(f.terms,g.terms,addVals);
--       if newHash === hashTable {} then newHash = (promote(0,f.ring)).terms;
--       new A from hashTable {(symbol ring, f.ring),
--                             (symbol cache, new CacheTable from {("isReduced",false)}),
--                             (symbol terms, newHash)}   
--    );

--    A ? A := (f,g) -> (
--       m := first pairs (leadMonomial f).terms;
--       n := first pairs (leadMonomial g).terms;
--       m ? n
--    );

--    A * A := (f,g) -> (
--       newHash := removeZeroes combine(f.terms,g.terms,multKeys,multVals,addVals);
--       if newHash === hashTable {} then newHash = (promote(0,f.ring)).terms;
--       new A from hashTable {(symbol ring, f.ring),
--                             (symbol cache, new CacheTable from {("isReduced",false)}),
--                             (symbol terms, newHash)}
--    );

--    A ^ ZZ := (f,n) -> product toList (n:f);

--    R * A := (r,f) -> promote(r,A)*f;
--    A * R := (f,r) -> r*f;
--    QQ * A := (r,f) -> promote(r,A)*f;
--    A * QQ := (f,r) -> r*f;
--    ZZ * A := (r,f) -> promote(r,A)*f;

--    A * ZZ := (f,r) -> r*f;
--    A - A := (f,g) -> f + (-1)*g;
--    - A := f -> (-1)*f;
--    A + ZZ := (f,r) -> f + promote(r,A);
--    ZZ + A := (r,f) -> f + r;
--    A + QQ := (f,r) -> f + promote(r,A);
--    QQ + A := (r,f) -> f + r;
--    A + R  := (f,r) -> f + promote(r,A);
--    R + A  := (r,f) -> f + r;
   
--    A ? A := (f,g) -> (leadTestMonomial f) ? (leadTestMonomial g);

--    A == A := (f,g) -> #(f.terms) == #(g.terms) and (sort pairs f.terms) == (sort pairs g.terms);
--    A == ZZ := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
--    ZZ == A := (n,f) -> f == n;
--    A == QQ := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
--    QQ == A := (n,f) -> f == n;
--    A == R := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
--    R == A := (n,f) -> f == n;

   A
)
new TestPolynomialRing from List := (TestPolynomialRing, inits) -> new TestPolynomialRing of TestRingElement from new HashTable from inits

Ring List := (R, varList) -> (
   -- get the symbols associated to the list that is passed in, in case the variables have been used earlier.
   print("BOX 1");
   if #varList == 0 then error "Expected at least one variable.";
   if #varList == 1 and class varList#0 === Sequence then varList = toList first varList;
   varList = varList / baseName;
   print("BOX 2");
   A := new TestPolynomialRing from {(symbol generators) => {},
                                   (symbol generatorSymbols) => varList,
                                   (symbol degreesRing) => degreesRing 1,
				   (symbol CoefficientRing) => R,
                                   (symbol cache) => new CacheTable from {},
				   (symbol baseRings) => {ZZ},
                                   BergmanRing => false
                };
    print("BOX 3");
    print (varList);
    newGens := apply(varList, v -> v);

    -- newGens := apply(varList, v -> v <- putInRing({v},A,1));
    print(newGens);

   A#(symbol generators) = newGens;
   A
)
V:= QQ{x}
use V
x
ZZ[x]
use TestRing
x
R

numgens R
generators R
coefficientRing R
isHomogeneous R
x

-- setWeights = method()
-- setWeights (NCRing,List) := (A,weightList) -> (
--    gensA := A.generatorSymbols;
--    A#(symbol weights) = new HashTable from apply(#gensA, i -> (gensA#i,weightList#i));
--    A
-- )



TestRing (Ring,ZZ,List) := (R,skewElt,varList) -> (
    R
) 
ZZ[x]
-- R = new TestRing from PolynomialRing := (TestRing, origPolyRing) -> new TestRing from PolynomialRing;
new TestRing from PolynomialRing := (TestRing, origPolyRing) -> new TestRing from PolynomialRing;

new TestRing from ZZ[x]

ZZ[x]

-- new NCPolynomialRing from List := (NCPolynomialRing, inits) -> new NCPolynomialRing of NCRingElement from new HashTable from inits

-- Create a test NC Ring.



----------------------------------TESTING PARTS OF NC ALGEBRA

