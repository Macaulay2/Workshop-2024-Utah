needsPackage "NormalToricVarieties"

divisorsToVariety := method();

--check that underlying varieties of given divisors are same

divisorsToVariety List := (listOfDivisors) -> (

inputVariety := variety(listOfDivisors#0);

-- 

rayList := {};
coneList := {};
outputVariety := normalToricVariety(rayList , coneList);

);

end ------------------------------

rayList = {{-1} , {1}}
coneList = {{1}, {0}}
X = normalToricVariety (rayList, coneList)
dim X
rays X
X1 = toricProjectiveSpace 1
X == X1
rays X1 == rays X
max X1
rayListH = { {1 , 0} , { 0 , 1} , {-1, 3}, {0 , -1}}
coneListH = { { 0 , 1} , { 1, 2} , { 2 , 3} , {3 , 0} }
Y = normalToricVariety (rayListH , coneListH)
Y1 = hirzebruchSurface 3
max Y1

D1 = toricDivisor ( {0 , 7} , X)
D0 = toricDivisor ( { 0 , 0}, X)
