TEST /// 
-*
needsPackage "ToricExtras"
*-
rayList = {{1} , {-1}}
coneList = {{1}, {0}}
PP1 = normalToricVariety (rayList, coneList)
D0 = toricDivisor ( { 0 , 0}, PP1)
D1 = toricDivisor ( {0 , 7} , PP1)
testH7 = projectivizationOfBundle({D0, D1})

assert(isWellDefined testH7)
assert(isSmooth testH7)
assert(isProjective testH7)
-- picardGroup will always be plus 1
--assert(rank picardGroup V == 3)
-- ask for dimension
-- check hirze is the same
///
