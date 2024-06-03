needsPackage "NormalToricVarieties"

divisorsToVariety := method();

--check that underlying varieties of given divisors are same

divisorsToVariety List := (listOfDivisors) -> (

inputVariety := variety(listOfDivisors#0);


-- listOfDivisors = { {0 , 7, 9}, {8 , 0 , 1} }
-- extract coefficients
--  get r = numberOfDivisors
-- build e_i
-- set e_0 = -e_1 - ... - e_r
-- Construct Fi's: Fi = cone(e_0, ... , exclude e_i , ... , e_r )
-- 

numberOfDivisors := #listOfDivisors;

inputConeList := max(inputVariety);
inputRayList := rays(inputVariety);

-- for each cone sigma in inputConeList, for each ray rho in sigma
-- construct cone(inputRayList#rho +  )
-- 

rayList := {};
coneList := {};
outputVariety := normalToricVariety(rayList , coneList);

);

end ------------------------------

restart
needsPackage "NormalToricVarieties"
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
variety(D1)

-----
PP2RayList = {
    {1,0}, {0,1}, {-1,-1}
}
PP2ConeList = {
    {0,1},
    {1,2},
    {2,0}
}
PP2 = normalToricVariety(PP2RayList, PP2ConeList)
D1 = toricDivisor ( {0 , 7, 9} , PP2)
D0 = toricDivisor ( { 8 , 0, 1}, PP2)

listOfDivisors = { D0, D1 }
inputVariety = variety listOfDivisors#0
r = #listOfDivisors - 1
oldDim = #entries(listOfDivisors#0)
-- build e_i
-- get original rays and then add zeros to the end
origRays = rays inputVariety
-- whichever side we add 0s to for the original rays, we need
-- to add to the complementary side for the e_is
origRaysInBigSpace = rays | { }