--Loading packages
restart
loadPackage ("NumberFields", Reload=>true)
-- loadPackage ("Permutations", Reload=>true)
loadPackage ("TestNCAlgebra", Reload=>true)
-- loadPackage ("AssociativeAlgebras", Reload=>true)


-- Play around with creating my own ring. What does this even mean?
A = ZZ{x,y,z}
f = y*z + z*y - x^2

TestRing = new Type of Ring
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