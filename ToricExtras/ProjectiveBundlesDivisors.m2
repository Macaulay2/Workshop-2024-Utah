needsPackage "NormalToricVarieties"

projectivizationOfBundle = method()

--check that underlying varieties of given divisors are same

projectivizationOfBundle List := (listOfDivisors) -> (

inputVariety := variety(listOfDivisors#0);


-- listOfDivisors = { {0 , 7, 9}, {8 , 0 , 1} }
-- extract coefficients
--  get r = numberOfDivisors
-- build e_i
-- set e_0 = -e_1 - ... - e_r
-- Construct Fi's: Fi = cone(e_0, ... , exclude e_i , ... , e_r )
-- 

-- make list of coefficients of divisors 
listOfDivisorCoefficients := apply(listOfDivisors, i -> entries(i));

extendingDimension := #listOfDivisors - 1;

inputConeList := max(inputVariety);
inputRayList := rays(inputVariety);

numberOfInputRays := #inputRayList;

-- get dimension of N_R
inputFanDimension := #(inputRayList#0);

-- padd with zeros to include the rays of the original fan in N_R x R^r
extendedInputRayList := {};

for oldRay in inputRayList do(

    extendedInputRayList = append(extendedInputRayList, (oldRay | toList (extendingDimension : 0)));

);


--build rayliftingCoeffs

rayLiftingCoefficients := {};

for i from 1 to extendingDimension do(

rayLiftingCoefficients = append(rayLiftingCoefficients , listOfDivisorCoefficients#i - listOfDivisorCoefficients#0);

);

extendedStandardBasisMatrix := map(ZZ^inputFanDimension, ZZ^extendingDimension, 0) || id_(ZZ^extendingDimension);

--extendedRayLiftingMatrix := map(ZZ^inputFanDimension, ZZ^extendingDimension, 0) || matrix(rayLiftingCoefficients);


-- add lifted original rays to new ray list
newRayList := {};

for j from 0 to (#extendedInputRayList - 1) do(

    oldRay := extendedInputRayList#j;

    for i from 0 to (extendingDimension - 1) do(
        oldRay = oldRay + entries((rayLiftingCoefficients#i)#j*extendedStandardBasisMatrix_i);

    );
    
    newRayList = append(newRayList , oldRay);
    
    
    );

-- make e_0 and add to new ray list

eZero := entries((-extendedStandardBasisMatrix)*vector(toList (extendingDimension : 1)));

newRayList = append(newRayList , eZero);

--add padded standard basis vectors to new ray list

for i from 0 to (extendingDimension - 1) do(
        newRayList = append(newRayList , entries(extendedStandardBasisMatrix_i));

    );

 

newConeList := {};

shiftedEis := toList(0..extendingDimension) + toList ((extendingDimension + 1) : numberOfInputRays);

for i from 0 to (#inputConeList - 1) do(

    for j from 0 to extendingDimension do(

        newConeList = append(newConeList , inputConeList#i | drop(shiftedEis, {j , j}) );

    );




);

--newVariety := normalToricVariety(newRayList , newConeList);

normalToricVariety(newRayList , newConeList)
);

end ------------------------------

restart
needsPackage "NormalToricVarieties"
load("ProjectiveBundlesDivisors.m2")
rayList = {{1} , {-1}}
coneList = {{1}, {0}}
X = normalToricVariety (rayList, coneList)
D0 = toricDivisor ( { 0 , 0}, X)
D1 = toricDivisor ( {0 , 7} , X)
L = projectivizationOfBundle({D0, D1})

rayListY = {{1 , 0}, {0 , 1}, {-1, -1}}
coneListY = {{0, 1}, {1, 2}, {2 , 0}}
Y = normalToricVariety (rayListY, coneListY)
D0Y = toricDivisor ( { 9, 3 , 2}, Y)
D1Y = toricDivisor ( {1 , 4, 7} , Y)
L = projectivizationOfBundle({D0Y, D1Y})

dim X
rays X
X1 = toricProjectiveSpace 1
rays X1 == rays X
max X1
rayListH = { {1 , 0} , { 0 , 1} , {-1, 3}, {0 , -1}}
coneListH = { { 0 , 1} , { 1, 2} , { 2 , 3} , {3 , 0} }
Y = normalToricVariety (rayListH , coneListH)
Y1 = hirzebruchSurface 3
max Y1


variety(D1)

divisorsToVariety()

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